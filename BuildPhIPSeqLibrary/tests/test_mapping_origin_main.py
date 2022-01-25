import os
from collections import Counter
from unittest import TestCase, mock

from BuildPhIPSeqLibrary.config import MOCK_DATA_DIR
from BuildPhIPSeqLibrary.mapping_origin_main import map_single_oligo_to_sequences_list, run_mapping_of_all_files
from BuildPhIPSeqLibrary.read_pipeline_files import read_oligo_sequences_to_file


class Test(TestCase):

    def test_map_single_oligo_to_sequences_list(self):
        oligo_id = 'oligo_1'
        oligo = 'LLL'
        sequences_dict = {'seq_1': 'LLLAMR', 'seq_2': 'LLLLADWEE'}
        ret_id, ret_locations = map_single_oligo_to_sequences_list(oligo_id, oligo, sequences_dict)
        self.assertEqual(len(ret_locations), 3)
        counter = Counter(list(map(lambda ret: ret[0], ret_locations)))
        self.assertEqual(counter['seq_1'], 1)
        self.assertEqual(counter['seq_2'], 2)

    @mock.patch('BuildPhIPSeqLibrary.read_pipeline_files.SEQUENCES_IDS_FILE',
                os.path.join(MOCK_DATA_DIR, 'PipelineFiles', 'sequences_ids.csv'))
    @mock.patch('BuildPhIPSeqLibrary.read_pipeline_files.OLIGO_SEQUENCES_FILE',
                os.path.join(MOCK_DATA_DIR, 'PipelineFiles', 'oligos_sequence.csv'))
    @mock.patch('pandas.DataFrame.to_csv')
    def test_run_mapping_of_all_files(self, to_csv_mock):
        ret = run_mapping_of_all_files()
        self.assertTrue(
            ret.apply(lambda row: all(map(lambda origin: origin in row['mapped'], row['origins'])), axis=1).all())
        self.assertFalse(ret['origins'].eq(ret['mapped']).all())
        aa_oligos = read_oligo_sequences_to_file()
        self.assertEqual(len(ret), len(aa_oligos))
