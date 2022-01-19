import os
import shutil
from unittest import TestCase, mock

from numpy.random import choice

from BuildPhIPSeqLibrary.config import OLIGO_AA_LENGTH, AMINO_ACIDS, OLIGO_AA_OVERLAP, MOCK_DATA_DIR
from BuildPhIPSeqLibrary.read_pipeline_files import read_sequence_ids_file
from BuildPhIPSeqLibrary.split_sequences_to_oligos import split_single_sequence_to_oligos, split_sequences_to_oligos, \
    merge_sequences


def create_random_aa(len_of_seq):
    return ''.join(choice(AMINO_ACIDS, len_of_seq))


class Test(TestCase):
    def setUp(self) -> None:
        shutil.rmtree(os.path.join(MOCK_DATA_DIR, 'Output'), ignore_errors=True)
        os.makedirs(os.path.join(MOCK_DATA_DIR, 'Output'))

    def doCleanups(self) -> None:
        shutil.rmtree(os.path.join(MOCK_DATA_DIR, 'Output'), ignore_errors=True)
        os.makedirs(os.path.join(MOCK_DATA_DIR, 'Output'))

    def test_split_single_sequence_to_oligos(self):
        # Test a short sequence:
        sequence = create_random_aa(OLIGO_AA_LENGTH - 4)
        ret = split_single_sequence_to_oligos(sequence)
        self.assertEqual(len(ret), 1)
        self.assertIn('*', ret.iloc[0]['oligo_aa_sequence'])
        self.assertEqual(0, ret.iloc[0]['position'])
        self.assertEqual(len(ret.iloc[0]['oligo_aa_sequence']), OLIGO_AA_LENGTH)

        # Test a sequence of length OLIGO_AA_LENGTH
        sequence = create_random_aa(OLIGO_AA_LENGTH)
        ret = split_single_sequence_to_oligos(sequence)
        self.assertEqual(len(ret), 1)
        self.assertNotIn('*', ret.iloc[0]['oligo_aa_sequence'])
        self.assertEqual(0, ret.iloc[0]['position'])
        self.assertEqual(len(ret.iloc[0]['oligo_aa_sequence']), OLIGO_AA_LENGTH)

        # Test a sequence of length 2*OLIGO_AA_LENGTH-OLIGO_AA_OVERLAP-4
        sequence = create_random_aa(2 * OLIGO_AA_LENGTH - OLIGO_AA_OVERLAP - 4)
        ret = split_single_sequence_to_oligos(sequence)
        self.assertEqual(len(ret), 2)
        self.assertFalse(ret['oligo_aa_sequence'].str.contains('*', regex=False).any())
        self.assertEqual(0, ret.iloc[0]['position'])
        self.assertGreater(OLIGO_AA_LENGTH - OLIGO_AA_OVERLAP, ret.iloc[1]['position'])
        self.assertTrue(ret['oligo_aa_sequence'].apply(len).eq(OLIGO_AA_LENGTH).all())

        # Test a sequence of length 2*OLIGO_AA_LENGTH-OLIGO_AA_OVERLAP
        sequence = create_random_aa(2 * OLIGO_AA_LENGTH - OLIGO_AA_OVERLAP)
        ret = split_single_sequence_to_oligos(sequence)
        self.assertEqual(len(ret), 2)
        self.assertFalse(ret['oligo_aa_sequence'].str.contains('*', regex=False).any())
        self.assertEqual(0, ret.iloc[0]['position'])
        self.assertEqual(OLIGO_AA_LENGTH - OLIGO_AA_OVERLAP, ret.iloc[1]['position'])
        self.assertTrue(ret['oligo_aa_sequence'].apply(len).eq(OLIGO_AA_LENGTH).all())

    def test_split_sequences_to_oligos(self):
        # Test one sequence
        sequence_1 = create_random_aa(OLIGO_AA_LENGTH + OLIGO_AA_OVERLAP)
        sequence_dict = {'seq_1': sequence_1}
        ret = split_sequences_to_oligos(sequence_dict)
        self.assertEqual(len(ret), 2)
        self.assertTrue(ret['origins'].apply(len).eq(1).all())

        # Test two sequences with some overlap
        sequence_2 = sequence_1[:OLIGO_AA_LENGTH] + create_random_aa(OLIGO_AA_OVERLAP)
        sequence_dict = {'seq_1': sequence_1, 'seq_2': sequence_2}
        ret = split_sequences_to_oligos(sequence_dict)
        self.assertEqual(len(ret), 3)
        self.assertEqual(ret['origins'].apply(len).max(), 2)

    @mock.patch('BuildPhIPSeqLibrary.split_sequences_to_oligos.OLIGO_SEQUENCES_FILE',
                os.path.join(MOCK_DATA_DIR, 'Output', 'oligos_sequence.csv'))
    @mock.patch('BuildPhIPSeqLibrary.read_pipeline_files.OLIGO_SEQUENCES_FILE',
                os.path.join(MOCK_DATA_DIR, 'Output', 'oligos_sequence.csv'))
    @mock.patch('BuildPhIPSeqLibrary.read_pipeline_files.SEQUENCES_IDS_FILE',
                os.path.join(MOCK_DATA_DIR, 'Output', 'sequences_ids.csv'))
    def test_merge_sequences(self):
        # Simple case
        sequence_1 = create_random_aa(OLIGO_AA_LENGTH + OLIGO_AA_OVERLAP)
        sequence_dict = {'seq_1': sequence_1}
        oligos_df = split_sequences_to_oligos(sequence_dict)
        ret, _ = merge_sequences(oligos_df, sequence_dict)
        self.assertEqual(len(ret), 2)
        self.assertTrue(ret['origins'].apply(len).eq(1).all())

        # Add the existing sequence
        sequence_dict = {'seq_2': sequence_1}
        oligos_df = split_sequences_to_oligos(sequence_dict)
        ret, _ = merge_sequences(oligos_df, sequence_dict)
        self.assertEqual(len(ret), 2)
        self.assertTrue(ret['origins'].apply(len).eq(2).all())

        # Add shifted oligo
        sequence_3 = create_random_aa(2) + sequence_1[:OLIGO_AA_LENGTH + 2]
        sequence_dict = {'seq_3': sequence_3}
        oligos_df = split_sequences_to_oligos(sequence_dict)
        ret, _ = merge_sequences(oligos_df, sequence_dict)
        self.assertEqual(len(ret), 4)
        self.assertEqual(ret['origins'].apply(len).max(), 2)

    @mock.patch('BuildPhIPSeqLibrary.split_sequences_to_oligos.OLIGO_SEQUENCES_FILE',
                os.path.join(MOCK_DATA_DIR, 'Output', 'oligos_sequence.csv'))
    @mock.patch('BuildPhIPSeqLibrary.read_pipeline_files.OLIGO_SEQUENCES_FILE',
                os.path.join(MOCK_DATA_DIR, 'Output', 'oligos_sequence.csv'))
    @mock.patch('BuildPhIPSeqLibrary.read_pipeline_files.SEQUENCES_IDS_FILE',
                os.path.join(MOCK_DATA_DIR, 'PipelineFiles', 'sequences_ids.csv'))
    def test_merge_and_map_sequences_with_existing_files(self):
        sequence_dict = read_sequence_ids_file()['AA_sequence'].to_dict()
        oligos_df = split_sequences_to_oligos(sequence_dict)
        ret, _ = merge_sequences(oligos_df, sequence_dict)
        self.assertTrue(ret['origins'].apply(len).eq(1).all())

        #
        sequence_1 = create_random_aa(OLIGO_AA_LENGTH + OLIGO_AA_OVERLAP)
        sequence_dict = {'seq_10': sequence_1, 'seq_11': sequence_dict['seq_5'][2:-2]}
        oligos_df = split_sequences_to_oligos(sequence_dict)
        with mock.patch('BuildPhIPSeqLibrary.split_sequences_to_oligos.OLIGO_SEQUENCES_FILE',
                        os.path.join(MOCK_DATA_DIR, 'PipelineFiles', 'oligos_sequence.csv')):
            ret, _ = merge_sequences(oligos_df, sequence_dict)
        self.assertTrue(ret['origins'].apply(len).eq(1).all())
