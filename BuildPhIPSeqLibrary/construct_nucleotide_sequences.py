import logging
import math
import time
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd

cnt_tries = 0

from BuildPhIPSeqLibrary.config import AMINO_INFO, BARCODE_IN_5_PRIME_END, BARCODE_NUC_LENGTHS, RESTRICTED_SEQUENCES, \
    UNCONVERTED_SEQUENCES_FILE, BARCODED_NUC_FILE, OLIGO_AA_LENGTH, NUM_MAPPING_THREADS, OUTPUT_DIR
from BuildPhIPSeqLibrary.read_pipeline_files import read_barcoded_nucleotide_files, read_unconverted_sequences


def has_no_restricted_sequences(nuc_seq):
    ret = all(list(map(lambda restricted: restricted not in nuc_seq, RESTRICTED_SEQUENCES)))
    return ret


def code_one_aa_sequence_to_nuc(aa_seq, num_tries=100, by_codon_frequencies=True, log=True):
    amino_info = AMINO_INFO.set_index('amino_acid')
    for _ in range(num_tries):
        ret = ''
        for aa in aa_seq:
            amino_acid_df = amino_info.loc[aa]
            if isinstance(amino_acid_df, pd.DataFrame):
                ret += np.random.choice(amino_acid_df['codon'].values,
                                        p=amino_acid_df[
                                            'corrected_relative_frequency'].values if by_codon_frequencies else None)
            else:
                ret += amino_acid_df['codon']
        # Ensure no restricted sequences are in the nucleotide sequence
        if has_no_restricted_sequences(ret):
            return ret
    if log:
        logging.debug(
            f"Failed to convert sequence from amino acids to nucleotides without adding restricted sequences "
            f"{aa_seq}")
    return None


def iterative_barcode_construction(aa_to_recode, nuc_prefix, nuc_suffix, existing_barcodes):
    """
    Iteratively creates a subset of the oligo.
    :param aa_to_recode:
    :param nuc_prefix:
    :param nuc_suffix:
    :param existing_barcodes:
    :return:
    """
    global cnt_tries

    if not BARCODE_IN_5_PRIME_END:
        start_pos = 3 * len(aa_to_recode) + len(nuc_prefix) + len(nuc_suffix) - sum(BARCODE_NUC_LENGTHS)
    if len(aa_to_recode) == 0:
        ret = nuc_prefix + nuc_suffix
        cnt_tries += 1
        if has_no_restricted_sequences(ret):
            return ret
        return None
    try:
        codon_opts = AMINO_INFO.set_index('amino_acid').loc[aa_to_recode[0]]['codon'].values
    except:
        codon_opts = [AMINO_INFO.set_index('amino_acid').loc[aa_to_recode[0]]['codon']]
    for codon in codon_opts:
        bad_opt = False
        if BARCODE_IN_5_PRIME_END:
            new_nuc_prefix = nuc_prefix + codon
            for i in range(len(BARCODE_NUC_LENGTHS)):
                if (len(nuc_prefix) < sum(BARCODE_NUC_LENGTHS[:i + 1])) and \
                        (len(new_nuc_prefix) >= sum(BARCODE_NUC_LENGTHS[:i + 1])):
                    barcode_i = new_nuc_prefix[:sum(BARCODE_NUC_LENGTHS[:i + 1])][-BARCODE_NUC_LENGTHS[i]:]
                    if barcode_i in existing_barcodes[f"barcode_{i}"].values:
                        bad_opt = True
                        break
            if not bad_opt:
                ret = iterative_barcode_construction(aa_to_recode[1:], new_nuc_prefix, nuc_suffix, existing_barcodes)
                if ret is not None:
                    return ret
        else:
            new_nuc_prefix = nuc_prefix + codon
            for i in range(len(BARCODE_NUC_LENGTHS))[::-1]:
                if ((len(nuc_prefix) - start_pos) < sum(BARCODE_NUC_LENGTHS[i:])) and \
                        ((len(new_nuc_prefix) - start_pos) >= sum(BARCODE_NUC_LENGTHS[i:])):
                    barcode_i = new_nuc_prefix[start_pos + sum(BARCODE_NUC_LENGTHS[i + 1:]):][:BARCODE_NUC_LENGTHS[i]]
                    if barcode_i in existing_barcodes[f"barcode_{i}"].values:
                        bad_opt = True
                        break
            if not bad_opt:
                ret = iterative_barcode_construction(aa_to_recode[1:], new_nuc_prefix, '', existing_barcodes)
                if ret is not None:
                    return ret
    return None


def iterative_correction_of_single_barcode(possible_nuc_sequences, oligo_row, existing_barcodes):
    """
    Given an oligo with a single Barcode that needs correcting, fixes that single barcode iteratively.
    :param possible_nuc_sequences:
    :param oligo_row:
    :param existing_barcodes:
    :return:
    """
    for possible_nuc_sequence in possible_nuc_sequences:
        oligo_row['nuc_sequence'] = possible_nuc_sequence
        get_all_barcodes(oligo_row)
        test_pos = [oligo_row[f'barcode_{i}'] in existing_barcodes[f'barcode_{i}'].values for i in
                    range(len(BARCODE_NUC_LENGTHS))].index(True)
        pos = [int(sum(BARCODE_NUC_LENGTHS[:test_pos]) / 3), math.ceil(sum(BARCODE_NUC_LENGTHS[:test_pos + 1]) / 3)]
        if not BARCODE_IN_5_PRIME_END:
            pos = [OLIGO_AA_LENGTH - pos[1], OLIGO_AA_LENGTH - pos[0]]
        res = iterative_barcode_construction(oligo_row['oligo_aa_sequence'][pos[0]:pos[1]],
                                             oligo_row['nuc_sequence'][:pos[0] * 3],
                                             oligo_row['nuc_sequence'][pos[1] * 3:],
                                             existing_barcodes)
        if res is not None:
            oligo_row['nuc_sequence'] = res
            get_all_barcodes(oligo_row)
            return oligo_row
    return None


def get_all_barcodes(oligo_row):
    start_location = 0
    for i, barcode_length in enumerate(BARCODE_NUC_LENGTHS):
        oligo_row[f'barcode_{i}'] = get_barcode_from_nuc_seq(oligo_row['nuc_sequence'], start_location,
                                                             barcode_length)
        start_location += barcode_length
    return oligo_row


def create_new_nuc_sequence(oligo_row, existing_barcodes, num_tries=50, try_random_probs=False, bad_aa_seqs=[]):
    global cnt_tries

    aa_range_to_recode = math.ceil(sum(BARCODE_NUC_LENGTHS) / 3)
    if BARCODE_IN_5_PRIME_END:
        aa_to_recode = oligo_row['oligo_aa_sequence'][:aa_range_to_recode]
        nuc_seq_to_maintain = oligo_row['nuc_sequence'][3 * aa_range_to_recode:]
    else:
        aa_to_recode = oligo_row['oligo_aa_sequence'][-aa_range_to_recode:]
        nuc_seq_to_maintain = oligo_row['nuc_sequence'][:-3 * aa_range_to_recode]
    if aa_to_recode not in bad_aa_seqs:
        nuc_sequences_with_least_failed = []
        least_failed_barcodes = 6
        if try_random_probs:
            codon_prob_opts = [True, False]
        else:
            codon_prob_opts = [True]
        for by_codon_probabilities in codon_prob_opts:
            for _ in range(num_tries):
                new_nuc_barcode_area = code_one_aa_sequence_to_nuc(aa_to_recode, num_tries=1,
                                                                   by_codon_frequencies=by_codon_probabilities, log=False)
                if new_nuc_barcode_area is None:
                    continue
                if BARCODE_IN_5_PRIME_END:
                    oligo_row['nuc_sequence'] = new_nuc_barcode_area + nuc_seq_to_maintain
                else:
                    oligo_row['nuc_sequence'] = nuc_seq_to_maintain + new_nuc_barcode_area
                if has_no_restricted_sequences(oligo_row['nuc_sequence']):
                    start_location = 0
                    for i, barcode_length in enumerate(BARCODE_NUC_LENGTHS):
                        oligo_row[f'barcode_{i}'] = get_barcode_from_nuc_seq(oligo_row['nuc_sequence'], start_location,
                                                                             barcode_length)
                        start_location += barcode_length
                    if all([oligo_row[f'barcode_{i}'] not in existing_barcodes[f'barcode_{i}'].values for i in
                            range(len(BARCODE_NUC_LENGTHS))]):
                        return oligo_row, int(not by_codon_probabilities), None
                    number_of_failed_barcodes = sum(
                        [oligo_row[f'barcode_{i}'] in existing_barcodes[f'barcode_{i}'].values for i in
                         range(len(BARCODE_NUC_LENGTHS))])
                    if number_of_failed_barcodes < least_failed_barcodes:
                        least_failed_barcodes = number_of_failed_barcodes
                        nuc_sequences_with_least_failed = [oligo_row['nuc_sequence']]
                    elif number_of_failed_barcodes == least_failed_barcodes:
                        nuc_sequences_with_least_failed.append(oligo_row['nuc_sequence'])

        # If num_tries random cases failed and there is only a single barcode that fails
        # try to re-create the barcode iteratively.
        if least_failed_barcodes == 1:
            res = iterative_correction_of_single_barcode(nuc_sequences_with_least_failed, oligo_row, existing_barcodes)
            if res is not None:
                return oligo_row, 2, None
        logging.debug(
            f"Failed to barcode sequence from amino acids to nucleotides without adding restricted sequences "
            f"{aa_to_recode}")
    else:
        logging.debug(f"Skipping failed {aa_to_recode}")
        aa_to_recode = None

    oligo_row['nuc_sequence'] = None
    for i, barcode_length in enumerate(BARCODE_NUC_LENGTHS):
        oligo_row[f'barcode_{i}'] = None
    return oligo_row, 3, aa_to_recode


def get_barcode_from_nuc_seq(nuc_seq, start_location, barcode_length):
    if BARCODE_IN_5_PRIME_END:
        return nuc_seq[start_location:start_location + barcode_length]
    else:
        return nuc_seq[::-1][start_location:start_location + barcode_length][::-1]


def barcode_sequences(oligo_sequences):
    start_location = 0
    for i, barcode_length in enumerate(BARCODE_NUC_LENGTHS):
        oligo_sequences[f'barcode_{i}'] = oligo_sequences['nuc_sequence'].apply(
            lambda nuc: get_barcode_from_nuc_seq(nuc, start_location, barcode_length))
        start_location += barcode_length
    existing_barcodes = read_barcoded_nucleotide_files()
    cols = existing_barcodes.columns

    inds = oligo_sequences.index
    for i, _ in enumerate(BARCODE_NUC_LENGTHS):
        inds = oligo_sequences.loc[inds][~np.isin(oligo_sequences.loc[inds]['barcode_%d' % i],
                                                  existing_barcodes['barcode_%d' % i])].index
    for i, _ in enumerate(BARCODE_NUC_LENGTHS):
        inds = oligo_sequences.loc[inds].drop_duplicates(['barcode_%d' % i]).index
    print("%d of %d added with no correction" % (len(inds), len(oligo_sequences)))

    existing_barcodes = pd.concat([existing_barcodes, oligo_sequences.loc[inds, cols]])
    oligo_sequences = oligo_sequences.drop(inds)

    uncoded_oligos = []
    bad_aa_seqs = []
    # Add barcoded oligos one at a time
    cnt = [0, 0, 0, 0, 0]
    for oligo_id, oligo_row in oligo_sequences.iterrows():
        if (cnt[0] % 100) == 0:
            print("At %d of %d (%d recode, %d random probs, %d correctd, %d failed)" %
                  (cnt[0], len(oligo_sequences), cnt[1], cnt[2], cnt[3], cnt[4]),
                  time.ctime())
        cnt[0] += 1
        if all([oligo_row[f'barcode_{i}'] not in existing_barcodes[f'barcode_{i}'].values for i in
                range(len(BARCODE_NUC_LENGTHS))]):
            existing_barcodes = existing_barcodes.append(oligo_row[cols])
        else:
            # Try to create a new barcode section
            new_oligo_row, ret, bad_aa = create_new_nuc_sequence(oligo_row, existing_barcodes, bad_aa_seqs=bad_aa_seqs)
            cnt[ret + 1] += 1
            if new_oligo_row['nuc_sequence'] is not None:
                existing_barcodes = existing_barcodes.append(new_oligo_row[cols])
            else:
                uncoded_oligos.append(new_oligo_row.copy())
            if bad_aa is not None:
                bad_aa_seqs.append(bad_aa)
    print("Finished (%d recode, %d random probs, %d correctd, %d failed)" % (cnt[1], cnt[2], cnt[3], cnt[4]),
          time.ctime())
    # Save barcodes
    existing_barcodes.to_csv(BARCODED_NUC_FILE)
    if len(uncoded_oligos) > 0:
        return existing_barcodes, pd.concat(uncoded_oligos, axis=1).T
    else:
        return existing_barcodes, pd.DataFrame()


def update_unconverted_oligos_file(unconverted_oligos):
    current_unconverted = read_unconverted_sequences()
    pd.concat([current_unconverted, unconverted_oligos[['oligo_aa_sequence']]]).to_csv(UNCONVERTED_SEQUENCES_FILE)


def aa_to_nuc(oligos_aa_sequences):
    pool = Pool(NUM_MAPPING_THREADS)
    oligos_aa_sequences['nuc_sequence'] = pool.starmap(code_one_aa_sequence_to_nuc,
                                                       oligos_aa_sequences[['oligo_aa_sequence']].values)
    pool.close()
    unconverted_oligos = oligos_aa_sequences[oligos_aa_sequences['nuc_sequence'].isnull()]
    if len(unconverted_oligos) > 0:
        update_unconverted_oligos_file(unconverted_oligos)
    oligos_aa_sequences = oligos_aa_sequences[oligos_aa_sequences['nuc_sequence'].notnull()].copy()
    ret, unconverted_oligos = barcode_sequences(oligos_aa_sequences)
    if len(unconverted_oligos) > 0:
        update_unconverted_oligos_file(unconverted_oligos)
    return ret
