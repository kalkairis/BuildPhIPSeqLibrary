import os
from unittest import TestCase

from BuildPhIPSeqLibrary.config import MOCK_DATA_DIR
from BuildPhIPSeqLibrary.read_input_files import read_file
from BuildPhIPSeqLibrary.sequence_ids import add_sequences_to_files_list, is_amino_acid_sequence


class Test(TestCase):
    def test_add_sequences_to_files_list(self):
        filename = os.path.join(MOCK_DATA_DIR, 'Input', 'sample_input.csv')
        ret = add_sequences_to_files_list(read_file(filename), filename)
        self.assertEqual(len(ret), 21)

    def test_check_is_amino_acid_sequence(self):
        amino_acid = 'FLLMVKN'
        self.assertTrue(is_amino_acid_sequence(amino_acid))
        amino_acid = 'xFLfd'
        self.assertFalse(is_amino_acid_sequence(amino_acid))
        self.assertTrue(is_amino_acid_sequence(''))
