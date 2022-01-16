import pandas as pd
import os

from BuildPhIPSeqLibrary.config import seq_AA_col, seq_ID_col, AMINO_INFO, SEQUENCES_IDS_FILE
from BuildPhIPSeqLibrary.read_pipeline_files import read_sequence_ids_file


def is_amino_acid_sequence(sequence):
    amino_acid_letters = set(AMINO_INFO['amino_acid'].unique())
    amino_acid_letters.remove('*')
    return all(list(map(lambda letter: letter in amino_acid_letters, list(sequence))))


def add_sequences_to_files_list(sequences, filename, output_path=None):
    """
    This file adds a list of sequences from a newly read file to the sequences table.
    It will return only newly-added sequences.
    :param output_path:
    :param sequences: dictionary with sequence_ID -> AA_sequence
    :param filename: Name of file where sequences originated from.
        This allows us to make sequence IDs non-unique between files.
    :return: Newly added sequences with their seq_ID (running identifier, not identical to sequence_ID).
    """
    if output_path is None:
        output_path = SEQUENCES_IDS_FILE
    # Read existing table of sequences
    sequences_df = read_sequence_ids_file(output_path)
    if len(sequences_df)>0:
        running_index_ID = sequences_df.index.str.split('_').str[1].astype(int).max() + 1
    else:
        running_index_ID = 0
    sequences_dict = sequences_df.reset_index()[['seq_ID', seq_AA_col]].drop_duplicates().set_index(seq_AA_col)[
        'seq_ID'].to_dict()

    # Go over each sequence, if exists add with same ID else add with new ID
    new_sequences = []
    return_sequences = {}
    for sequence, sequence_AA in sequences.items():
        if sequence_AA in sequences_dict.keys():
            sequence_ID = sequences_dict[sequence_AA]
        else:
            assert is_amino_acid_sequence(
                sequence_AA), f"Sequence {sequence}, in {filename} is not an amino acid sequence"
            sequence_ID = '_'.join(['seq', str(running_index_ID)])
            running_index_ID += 1
            return_sequences[sequence_ID] = sequence_AA
        new_sequences.append(
            {'seq_ID': sequence_ID, seq_AA_col: sequence_AA, seq_ID_col: sequence, 'input_file': filename})
    new_sequences = pd.DataFrame(new_sequences).set_index('seq_ID')
    sequences_df = pd.concat([sequences_df, new_sequences])
    sequences_df.to_csv(output_path)
    return return_sequences
