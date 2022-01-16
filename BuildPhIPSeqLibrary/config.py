import os

import pandas as pd

if os.path.exists(__file__.replace('.py', '_local.py')):
    from config_local import *

# Directory paths
DATA_DIR = globals().get('DATA_DIR', os.path.join(os.path.dirname(os.path.dirname(__file__)), 'Data'))
INPUT_DIR = globals().get('INPUT_DIR', os.path.join(DATA_DIR, 'Input'))
OUTPUT_DIR = globals().get('OUTPUT_DIR', os.path.join(DATA_DIR, 'Output'))
MOCK_DATA_DIR = globals().get('MOCK_DATA_DIR', os.path.join(os.path.dirname(os.path.dirname(__file__)), 'Mock'))

# File paths
FILES_INPUT_HASH_FILE = globals().get('FILES_INPUT_HASH_FILE', os.path.join(OUTPUT_DIR, 'files_hash.csv'))
SEQUENCES_IDS_FILE = globals().get('SEQUENCES_IDS_FILE', os.path.join(OUTPUT_DIR, 'sequences_ids.csv'))
OLIGO_SEQUENCES_FILE = globals().get('OLIGO_SEQUENCES_FILE', os.path.join(OUTPUT_DIR, 'oligos_sequence.csv'))
BARCODED_NUC_FILE = globals().get('BARCODED_NUC_FILE', os.path.join(OUTPUT_DIR, 'barcoded_nuc_file.csv'))
UNCONVERTED_SEQUENCES_FILE = globals().get('UNCONVERTED_SEQUENCES_FILE',
                                           os.path.join(OUTPUT_DIR, 'unconverted_sequences.csv'))
ORDER_FILE = globals().get('ORDER_FILE', os.path.join(OUTPUT_DIR, 'order_file.txt'))

# Column names
seq_ID_col = 'sequence_ID'
seq_AA_col = 'AA_sequence'

# Oligo lengths
OLIGO_AA_LENGTH = globals().get('OLIGO_AA_LENGTH', 64)
OLIGO_AA_OVERLAP = globals().get('OLIGO_AA_OVERLAP', 20)
BARCODE_IN_5_PRIME_END = globals().get('BARCODE_IN_5_PRIME_END', True)  # TODO: check this with Thomas
BARCODE_NUC_LENGTHS = globals().get('BARCODE_NUC_LENGTHS', [15, 15, 15, 15, 15])
# TODO: check with Thomas this does not change if we choose to change the barcodes to the beginning of the oligo

# Amino acid coding and frequencies
AMINO_ACID_PATH = globals().get('AMINO_ACID_PATH', os.path.join(os.path.dirname(__file__), 'amino_acids_config.csv'))
# Codon usage frequencies from https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=83333&aa=1&style=N
AMINO_INFO = pd.read_csv(AMINO_ACID_PATH, index_col=0)
AMINO_ACIDS = AMINO_INFO[AMINO_INFO['amino_acid'].ne('*')]['amino_acid'].unique()
RESTRICTED_SEQUENCES = globals().get('RESTRICTED_SEQUENCES',
                                     ['GAATTC', 'AAGCTT', 'AGGAGG', 'CCTCCT', 'CTTTCA', 'TGAAAG',
                                      'AACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG', 'A' * 8, 'T' * 8, 'C' * 8, 'G' * 8,
                                      'AC' * 8, 'AG' * 8, 'AT' * 8, 'CA' * 8, 'CG' * 8, 'CT' * 8, 'GC' * 8, 'GA' * 8,
                                      'GT' * 8, 'TC' * 8, 'TG' * 8, 'TA' * 8])

# Promoters
PREFIX_PROMOTER = globals().get('PREFIX_PROMOTER', "GATGCGCCGTGGGAATTCT")
SUFFIX_PROMOTER = globals().get('SUFFIX_PROMOTER', "TGAAAGCTTGCCACCCGAC")