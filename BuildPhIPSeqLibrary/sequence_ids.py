import pandas as pd
import os

from BuildPhIPSeqLibrary.config import OUTPUT_DIR, seq_AA_col, seq_ID_col, AMINO_INFO


def is_amino_acid_sequence(sequence):
    amino_acid_letters = set(AMINO_INFO['amino_acid'].unique())
    amino_acid_letters.remove('*')
    return all(list(map(lambda letter: letter in amino_acid_letters, list(sequence))))


def add_sequences_to_files_list(sequences, filename):
    """
    This file adds a list of sequences from a newly read file to the sequences table.
    It will return only newly-added sequences.
    :param sequences: dictionary with sequence_ID -> AA_sequence
    :param filename: Name of file where sequences originated from.
        This allows us to make sequence IDs non-unique between files.
    :return: Newly added sequences with their seq_ID (running identifier, not identical to sequence_ID).
    """
    output_path = os.path.join(OUTPUT_DIR, 'sequences_ids.csv')
    # Read existing table of sequences
    if os.path.exists(output_path):
        sequences_df = pd.read_csv(output_path, index_col=0)
        running_index_ID = sequences_df.index.str.split('_').str[1].astype(int).max() + 1
    else:
        sequences_df = pd.DataFrame(columns=['seq_ID', seq_AA_col, seq_ID_col, 'input_file']).set_index('seq_ID')
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
    sequences_df = pd.concat([sequences_df, new_sequences])  # TODO: check axis
    sequences_df.to_csv(output_path)
    return return_sequences
