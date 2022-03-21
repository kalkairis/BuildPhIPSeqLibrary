import os
from unittest import TestCase

import numpy

from BuildPhIPSeqLibrary.ID_from_barcode import find_and_output
from BuildPhIPSeqLibrary.config import BARCODE_IN_5_PRIME_END, BARCODE_NUC_LENGTHS, MOCK_DATA_DIR
from BuildPhIPSeqLibrary.read_pipeline_files import read_barcoded_nucleotide_files


class Test(TestCase):
    def test_id_barcoding(self):
        num_rounds = 10
        df_barcodes = read_barcoded_nucleotide_files(
            os.path.join(MOCK_DATA_DIR, 'PipelineFiles', 'barcoded_nuc_file.csv'))

        for i in range(num_rounds):
            oli = numpy.random.choice(df_barcodes.index)
            if BARCODE_IN_5_PRIME_END:
                input_df = df_barcodes.loc[oli].nuc_sequence[:sum(BARCODE_NUC_LENGTHS)]
            else:
                input_df = df_barcodes.loc[oli].nuc_sequence[-sum(BARCODE_NUC_LENGTHS):][::-1]
            ID, indels, errs, _, _ = find_and_output(input_df, False, df_barcodes)
            self.assertEqual(ID, oli)
            self.assertTrue((indels == 0) and (errs == 0))

            pos = numpy.random.randint(0, sum(BARCODE_NUC_LENGTHS))
            ID, indels, errs, _, _ = find_and_output(input_df[:pos] + "N" + input_df[pos + 1:], False, df_barcodes)
            self.assertEqual(ID, oli)
            self.assertTrue((indels == 0) and (errs == 1))

            pos2 = numpy.random.randint(0, sum(BARCODE_NUC_LENGTHS))
            tmp_input = input_df[:pos] + "N" + input_df[pos + 1:]
            tmp_input = tmp_input[:pos2] + "N" + tmp_input[pos2 + 1:]
            ID, indels, errs, _, _ = find_and_output(tmp_input, False, df_barcodes)
            self.assertEqual(ID, oli)
            self.assertTrue((indels == 0) and (errs <= 2))

            ID, indels, errs, _, _ = find_and_output(input_df[:pos] + input_df[pos + 1:] + "N", True, df_barcodes)
            self.assertEqual(ID, oli)
            self.assertTrue(((indels <= 1) and ((indels + errs) <= 1)) or ((indels == 0) and (errs <= 2)))

            tmp_input = input_df[:pos] + "N" + input_df[pos + 1:]
            tmp_input = tmp_input[:pos2] + tmp_input[pos2 + 1:] + "N"
            ID, indels, errs, _, _ = find_and_output(tmp_input, True, df_barcodes)
            self.assertEqual(ID, oli)
            self.assertTrue((indels <= 1) and ((indels + errs) <= 2))

            tmp_input = input_df[:pos] + input_df[pos + 1:] + "N"
            tmp_input = tmp_input[:pos2] + tmp_input[pos2 + 1:] + "N"
            ID, indels, errs, _, _ = find_and_output(tmp_input, True, df_barcodes)
            self.assertEqual(ID, oli)
            self.assertTrue((indels <= 2) and ((indels + errs) <= 2))
