import logging
import os
import shutil
from unittest import TestCase, mock

import pandas
import pandas as pd

from BuildPhIPSeqLibrary.config import MOCK_DATA_DIR, seq_ID_col, seq_AA_col
from BuildPhIPSeqLibrary.read_input_files import get_input_files, read_file


class Test(TestCase):
    def setUp(self) -> None:
        self.mock_output = os.path.join(MOCK_DATA_DIR, 'Output')
        shutil.rmtree(self.mock_output, ignore_errors=True)
        os.makedirs(self.mock_output)

    def doCleanups(self) -> None:
        shutil.rmtree(self.mock_output, ignore_errors=True)

    @mock.patch('BuildPhIPSeqLibrary.read_input_files.INPUT_DIR', os.path.join(MOCK_DATA_DIR, 'Input'))
    @mock.patch('BuildPhIPSeqLibrary.read_input_files.OUTPUT_DIR', os.path.join(MOCK_DATA_DIR, 'Output'))
    def test_get_input_files(self):
        ret = get_input_files()
        self.assertEqual(len(ret), 2)
        self.assertTrue(os.path.exists(os.path.join(self.mock_output, 'files_hash.csv')))
        ret = get_input_files()
        self.assertEqual(len(ret), 0)
        with mock.patch('pandas.read_csv',
                        return_value=pandas.DataFrame(index=['fake_file'], data={'0': ["fake_file_hash"]})):
            ret = get_input_files()
        files_hash = pandas.read_csv(os.path.join(self.mock_output, 'files_hash.csv'), index_col=0)
        files_hash.iloc[0, 0] = 'fake_file_hash'
        with mock.patch('pandas.read_csv', return_value=files_hash):
            self.assertRaises(AssertionError, get_input_files)

    def test_read_file(self):
        file_path = os.path.join(MOCK_DATA_DIR, 'Input', 'sample_input.csv')
        ret = read_file(file_path)
        self.assertEqual(len(ret), 21)
        with mock.patch('pandas.read_csv',
                        return_value=pd.DataFrame({seq_ID_col: ['a', 'a', 'b'], seq_AA_col: ['S1', 'S2', 'S3']})):
            self.assertRaises(AssertionError, read_file, file_path)
        with mock.patch('pandas.read_csv', return_value=pd.DataFrame(columns=[seq_ID_col, seq_AA_col])):
            ret = read_file(file_path)
            self.assertEqual(len(ret), 0)
        with mock.patch('pandas.read_csv', return_value=pd.DataFrame()):
            self.assertRaises(KeyError, read_file, file_path)
