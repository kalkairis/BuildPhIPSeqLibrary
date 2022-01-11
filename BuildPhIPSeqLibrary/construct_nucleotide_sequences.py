import logging
import math

import numpy as np
import pandas as pd

from BuildPhIPSeqLibrary.config import AMINO_INFO, BARCODE_IN_5_PRIME_END, BARCODE_NUC_LENGTHS, RESTRICTED_SEQUENCES, \
    UNCONVERTED_SEQUENCES_FILE, BARCODED_NUC_FILE
from BuildPhIPSeqLibrary.read_pipeline_files import read_barcoded_nucleotide_files, read_unconverted_sequences


def has_no_restricted_sequences(nuc_seq):
    return all(list(map(lambda restricted: restricted not in nuc_seq, RESTRICTED_SEQUENCES)))


def code_one_aa_sequence_to_nuc(aa_seq, num_tries=10):
    amino_info = AMINO_INFO.set_index('amino_acid')
    for _ in range(num_tries):
        ret = ''
        for aa in aa_seq:
            amino_acid_df = amino_info.loc[aa]
            if isinstance(amino_acid_df, pd.DataFrame):
                ret += np.random.choice(amino_acid_df['codon'].values,
                                        p=amino_acid_df['corrected_relative_frequency'].values)
            else:
                ret += amino_acid_df['codon']
        # Ensure no restricted sequences are in the nucleotide sequence
        if has_no_restricted_sequences(ret):
            return ret
    logging.warning(f"Failed to convert sequence from amino acids to nucleotides without adding restricted sequences"
                    f"{aa}")
    return None


def create_new_nuc_sequence(oligo_row, existing_barcodes, num_tries=20):
    aa_range_to_recode = math.ceil(sum(BARCODE_NUC_LENGTHS) / 3)
    succeeded_recoding = False
    for _ in range(num_tries):
        if BARCODE_IN_5_PRIME_END:
            aa_to_recode = oligo_row['oligo_aa_sequence'][:aa_range_to_recode]
            new_nuc_barcode_area = code_one_aa_sequence_to_nuc(aa_to_recode)
            oligo_row['nuc_sequence'] = new_nuc_barcode_area + oligo_row['nuc_sequence'][len(new_nuc_barcode_area):]

        else:
            aa_to_recode = oligo_row['oligo_aa_sequence'][-aa_range_to_recode:]
            new_nuc_barcode_area = code_one_aa_sequence_to_nuc(aa_to_recode)
            oligo_row['nuc_sequence'] = oligo_row['nuc_sequence'][:-len(new_nuc_barcode_area)] + new_nuc_barcode_area
        if has_no_restricted_sequences(oligo_row['nuc_sequence']):
            start_location = 0
            for i, barcode_length in enumerate(BARCODE_NUC_LENGTHS):
                oligo_row[f'barcode_{i}'] = get_barcode_from_nuc_seq(oligo_row['nuc_sequence'], start_location,
                                                                     barcode_length)
                start_location += barcode_length
            if all([oligo_row[f'barcode_{i}'] not in existing_barcodes[f'barcode_{i}'] for i in
                    range(len(BARCODE_NUC_LENGTHS))]):
                return oligo_row
    # TODO: run the full option case
    oligo_row['nuc_sequence'] = None
    return oligo_row


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
    uncoded_oligos = []
    # Add barcoded oligos one at a time
    for oligo_id, oligo_row in oligo_sequences.iterrows():
        if False: #all([oligo_row[f'barcode_{i}'] not in existing_barcodes[f'barcode_{i}'] for i in
                #range(len(BARCODE_NUC_LENGTHS))]):
            existing_barcodes = existing_barcodes.append(oligo_row[cols])
        else:
            # Try to create a new barcode section
            new_oligo_row = create_new_nuc_sequence(oligo_row, existing_barcodes)
            if new_oligo_row['nuc_sequence'] is not None:
                existing_barcodes = existing_barcodes.append(new_oligo_row[cols])
            else:
                uncoded_oligos.append(new_oligo_row.copy())
    # Save barcodes
    existing_barcodes.to_csv(BARCODED_NUC_FILE)
    return oligo_sequences, pd.concat(uncoded_oligos, axis=1).T


def update_unconverted_oligos_file(unconverted_oligos):
    current_unconverted = read_unconverted_sequences()
    pd.concat([current_unconverted, unconverted_oligos['oligo_aa_sequence']]).to_csv(UNCONVERTED_SEQUENCES_FILE)


def aa_to_nuc(oligos_aa_sequences):
    oligos_aa_sequences['nuc_sequence'] = oligos_aa_sequences.oligo_aa_sequence.apply(code_one_aa_sequence_to_nuc)
    unconverted_oligos = oligos_aa_sequences[oligos_aa_sequences['nuc_sequence'].isnull()]
    if len(unconverted_oligos) > 0:
        update_unconverted_oligos_file(unconverted_oligos)
    oligos_aa_sequences = oligos_aa_sequences[oligos_aa_sequences['nuc_sequence'].notnull()].copy()
    ret, unconverted_oligos = barcode_sequences(oligos_aa_sequences)
    if len(unconverted_oligos) > 0:
        update_unconverted_oligos_file(unconverted_oligos)
    return ret
