import os
from unittest import TestCase, mock

import pandas

from BuildPhIPSeqLibrary.config import MOCK_DATA_DIR
from BuildPhIPSeqLibrary.read_input_files import read_file
from BuildPhIPSeqLibrary.sequence_ids import add_sequences_to_files_list, is_amino_acid_sequence


class Test(TestCase):
    def setUp(self) -> None:
        self.seqs_table = os.path.join(MOCK_DATA_DIR, 'Output', 'sequences_ids.csv')
        os.makedirs(os.path.join(MOCK_DATA_DIR, 'Output'), exist_ok=True)
        if os.path.exists(self.seqs_table):
            os.remove(self.seqs_table)

    def doCleanups(self) -> None:
        if os.path.exists(self.seqs_table):
            os.remove(self.seqs_table)

    @mock.patch('BuildPhIPSeqLibrary.sequence_ids.OUTPUT_DIR', os.path.join(MOCK_DATA_DIR, 'Output'))
    def test_add_sequences_to_files_list(self):
        filename = os.path.join(MOCK_DATA_DIR, 'Input', 'sample_input.csv')
        ret = add_sequences_to_files_list(read_file(filename), filename)
        self.assertEqual(len(ret), 21)
        seqs_table = pandas.read_csv(self.seqs_table)
        self.assertEqual(len(seqs_table), 21)
        self.assertEqual(seqs_table.seq_ID.nunique(), 21)
        filename = os.path.join(MOCK_DATA_DIR, 'Input', 'sample_input.csv')
        ret = add_sequences_to_files_list(read_file(filename), filename)
        self.assertEqual(len(ret), 0)
        seqs_table = pandas.read_csv(self.seqs_table)
        self.assertEqual(len(seqs_table), 42)
        self.assertEqual(seqs_table.seq_ID.nunique(), 21)

    @mock.patch('BuildPhIPSeqLibrary.sequence_ids.OUTPUT_DIR', os.path.join(MOCK_DATA_DIR, 'Output'))
    def test_check_is_amino_acid_sequence(self):
        amino_acid = 'FLLMVKN'
        self.assertTrue(is_amino_acid_sequence(amino_acid))
        amino_acid = 'xFLfd'
        self.assertFalse(is_amino_acid_sequence(amino_acid))
        self.assertTrue(is_amino_acid_sequence(''))
