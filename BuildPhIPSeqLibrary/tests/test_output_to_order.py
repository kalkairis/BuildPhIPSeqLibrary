import os
import shutil
from unittest import TestCase, mock

import pandas as pd

from BuildPhIPSeqLibrary.config import MOCK_DATA_DIR, PREFIX, SUFFIX
from BuildPhIPSeqLibrary.output_to_order import transfer_to_order


class Test(TestCase):
    def doCleanups(self) -> None:
        shutil.rmtree(os.path.join(MOCK_DATA_DIR, 'Output'), ignore_errors=True)
        os.makedirs(os.path.join(MOCK_DATA_DIR, 'Output'))

    @mock.patch('BuildPhIPSeqLibrary.output_to_order.ORDER_FILE',
                os.path.join(MOCK_DATA_DIR, 'Output', 'order_file.txt'))
    def test_transfer_to_order(self):
        output_oligos = pd.DataFrame(index=['oligo_1', 'oligo_2'],
                                     data={'nuc_sequence': ['AGTAGTAGT', 'TGATGATGA'], 'barcode_0': ['A', 'B']})
        transfer_to_order(output_oligos)
        num_lines = 0
        with open(os.path.join(MOCK_DATA_DIR, 'Output', 'order_file.txt'), 'r') as file:
            for line in file:
                self.assertIn(line.strip()[len(PREFIX):-len(SUFFIX)], output_oligos['nuc_sequence'].values)
                num_lines += 1
        self.assertEqual(num_lines, len(output_oligos))
        transfer_to_order(output_oligos.iloc[:1])
        num_lines = 0
        with open(os.path.join(MOCK_DATA_DIR, 'Output', 'order_file.txt'), 'r') as file:
            for line in file:
                self.assertIn(line.strip()[len(PREFIX):-len(SUFFIX)], output_oligos['nuc_sequence'].values)
                num_lines += 1
        self.assertEqual(num_lines, 1)
