import os.path

import pandas as pd

from BuildPhIPSeqLibrary.config import OLIGO_SEQUENCES_FILE, SEQUENCES_IDS_FILE, seq_AA_col, seq_ID_col, \
    BARCODED_NUC_FILE, BARCODE_NUC_LENGTHS, UNCONVERTED_SEQUENCES_FILE


def read_oligo_sequences_to_file(file_path=None):
    """
    Reads the oligo sequences file. If does not exist returns an empty DataFrame
    :param file_path: path of oligo sequence file. By default, OLIGO_SEQUENCES_FILE
    :return: pandas DataFrame of oligos sequences with origins and mappings
    """
    if file_path is None:
        file_path = OLIGO_SEQUENCES_FILE
    if os.path.exists(file_path):
        ret = pd.read_csv(file_path, index_col=0)
        ret['origins'] = ret['origins'].apply(eval)
        ret['mapped'] = ret['mapped'].apply(eval)
        return ret
    else:
        return pd.DataFrame(columns=['origins', 'mapped', 'oligo_aa_sequence', 'oligo_id']).set_index('oligo_id')


def read_sequence_ids_file(file_path=None):
    """
    Reads sequences IDS file
    :param file_path:
    :return:
    """
    if file_path is None:
        file_path = SEQUENCES_IDS_FILE
    if os.path.exists(file_path):
        sequences_df = pd.read_csv(file_path, index_col=0)
    else:
        sequences_df = pd.DataFrame(columns=['seq_ID', seq_AA_col, seq_ID_col, 'input_file']).set_index('seq_ID')
    return sequences_df


def read_unconverted_sequences(file_path=None):
    """
    Reads table of barcoded nucleotides and oligo IDs
    :param file_path:
    :return:
    """
    if file_path is None:
        file_path = UNCONVERTED_SEQUENCES_FILE
    if os.path.exists(file_path):
        unconverted_sequences_df = pd.read_csv(file_path, index_col=0)
    else:
        unconverted_sequences_df = pd.DataFrame(columns=['oligo_id', 'oligo_aa_sequence']).set_index('oligo_id')
    return unconverted_sequences_df


def read_barcoded_nucleotide_files(file_path=None):
    """
    Reads table of barcoded nucleotides and oligo IDs
    :param file_path:
    :return:
    """
    if file_path is None:
        file_path = BARCODED_NUC_FILE
    if os.path.exists(file_path):
        barcoded_nuc_df = pd.read_csv(file_path, index_col=0)
    else:
        barcoded_nuc_df = pd.DataFrame(columns=['oligo_id', 'nuc_sequence'] + list(
            map(lambda i: f"barcode_{i}", range(len(BARCODE_NUC_LENGTHS))))).set_index('oligo_id')
    return barcoded_nuc_df
