import time

import pandas as pd
from numpy.random import choice

from BuildPhIPSeqLibrary.config import OLIGO_AA_LENGTH, OLIGO_AA_OVERLAP, AMINO_ACIDS, OLIGO_SEQUENCES_FILE
from BuildPhIPSeqLibrary.read_pipeline_files import read_oligo_sequences_to_file


def split_single_sequence_to_oligos(sequence):
    """
    This function received a single AA sequence and splits it into a list of oligos.
    :param sequence: AA sequence as string
    :return: A pandas DataFrame containing position in oligo, AA sequence
    """
    ret = []
    if len(sequence) < OLIGO_AA_LENGTH:
        ret = [[0, sequence + '*' + ''.join(choice(AMINO_ACIDS, OLIGO_AA_LENGTH - len(sequence) - 1))]]
    elif len(sequence) == OLIGO_AA_LENGTH:
        ret = [[0, sequence]]
    else:
        for position in range(0, len(sequence), OLIGO_AA_LENGTH - OLIGO_AA_OVERLAP):
            if (position + OLIGO_AA_LENGTH) > len(sequence):
                if len(sequence) - OLIGO_AA_LENGTH <= ret[-1][0]:
                    continue
                ret.append([len(sequence) - OLIGO_AA_LENGTH, sequence[-OLIGO_AA_LENGTH:]])
            else:
                ret.append([position, sequence[position: position + OLIGO_AA_LENGTH]])
    ret = pd.DataFrame(data=ret, columns=['position', 'oligo_aa_sequence'])
    return ret


def split_sequences_to_oligos(id_to_sequences):
    """
    This function performs the cut of sequences into smaller oligos.
    It ensures overlap and oligo sizes as defined in config.py.
    :param id_to_sequences: A dictionary from sequence ID to AA sequences.
    :return: A dictionary of oligo sequence, list of origins
    """
    ret = []
    for seq_id, seq in id_to_sequences.items():
        single_seq_table = split_single_sequence_to_oligos(seq)
        single_seq_table['seq_id'] = seq_id
        ret.append(single_seq_table)
    ret = pd.concat(ret, axis=0, ignore_index=True)
    ret = ret.groupby('oligo_aa_sequence').apply(
        lambda col: list(zip(list(col['seq_id']), list(col['position'])))).to_frame().rename(columns={0: 'origins'})
    return ret


def merge_sequences(oligos_aa_sequences, id_to_sequences):
    print("Merge oligos sequences, this takes time", time.ctime())
    base_df = read_oligo_sequences_to_file()

    # Stage 1: add amino acid oligos already in origins into origins list and remove from new sequences
    existing_oligos = oligos_aa_sequences.index.intersection(base_df['oligo_aa_sequence'].values)
    for existing_oligo in existing_oligos:
        base_df.loc[base_df['oligo_aa_sequence'].eq(existing_oligo).idxmax(), 'origins'] += \
            oligos_aa_sequences.loc[existing_oligo]['origins']
    oligos_aa_sequences.drop(index=existing_oligos, inplace=True)

    # Stage 2: add running oligo index to newly added oligos
    oligos_aa_sequences.reset_index(inplace=True)
    running_index = 0 if len(base_df) == 0 else base_df.index.str.split('_').str[1].astype(int).max() + 1
    oligos_aa_sequences['oligo_id'] = list(
        map(lambda idx: f'oligo_{idx}', range(running_index, running_index + len(oligos_aa_sequences))))
    oligos_aa_sequences.set_index('oligo_id', inplace=True)

    # Stage 3: concatenate dfs
    ret = pd.concat([base_df, oligos_aa_sequences], axis=0, ignore_index=False)
    ret.to_csv(OLIGO_SEQUENCES_FILE)
    return ret, oligos_aa_sequences


def split_and_map_new_sequences(id_to_sequences):
    oligo_aa_sequences = split_sequences_to_oligos(id_to_sequences)
    all_sequences, new_sequences = merge_sequences(oligo_aa_sequences, id_to_sequences)
    return all_sequences, new_sequences
