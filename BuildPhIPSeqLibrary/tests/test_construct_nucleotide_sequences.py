import math
import os
import shutil
from multiprocessing.dummy import Pool
from unittest import TestCase, mock

import numpy as np
import pandas as pd

from BuildPhIPSeqLibrary.config import RESTRICTED_SEQUENCES, AMINO_ACIDS, AMINO_INFO, BARCODE_NUC_LENGTHS, \
    OLIGO_AA_LENGTH, MOCK_DATA_DIR
from BuildPhIPSeqLibrary.construct_nucleotide_sequences import has_no_restricted_sequences, code_one_aa_sequence_to_nuc, \
    get_barcode_from_nuc_seq, create_new_nuc_sequence, barcode_sequences, iterative_barcode_construction, \
    get_all_barcodes, iterative_correction_of_single_barcode, aa_to_nuc
from BuildPhIPSeqLibrary.read_pipeline_files import read_oligo_sequences_to_file


class Test(TestCase):
    def setUp(self) -> None:
        shutil.rmtree(os.path.join(MOCK_DATA_DIR, 'Output'), ignore_errors=True)
        os.makedirs(os.path.join(MOCK_DATA_DIR, 'Output'))

    def doCleanups(self) -> None:
        shutil.rmtree(os.path.join(MOCK_DATA_DIR, 'Output'), ignore_errors=True)
        os.makedirs(os.path.join(MOCK_DATA_DIR, 'Output'))

    def test_has_no_restricted_sequences(self):
        self.assertTrue(has_no_restricted_sequences('ASDSAFEWRSDSKLFJDSFIERJFLKDSJFSAKLDRJEFD'))
        for restricted in RESTRICTED_SEQUENCES:
            self.assertFalse(has_no_restricted_sequences(restricted))
            self.assertFalse(has_no_restricted_sequences(
                ''.join(np.random.choice(['A', 'G', 'C', 'T'], 10)) + restricted + ''.join(
                    np.random.choice(['A', 'G', 'C', 'T'], 4))))

    def test_code_one_aa_sequence_to_nuc(self):
        aa_to_nuc = AMINO_INFO.set_index('codon')['amino_acid'].to_dict()
        for _ in range(5):
            aa_seq = ''.join(np.random.choice(AMINO_ACIDS, 20))
            nuc_seq = code_one_aa_sequence_to_nuc(aa_seq)
            new_aa_seq = ''.join([aa_to_nuc[nuc_seq[i: i + 3]] for i in range(0, len(nuc_seq), 3)])
            self.assertEqual(aa_seq, new_aa_seq)

    def test_create_new_nuc_sequence(self):
        with mock.patch('BuildPhIPSeqLibrary.construct_nucleotide_sequences.BARCODE_IN_5_PRIME_END', True):
            oligo_rows = pd.DataFrame(index=['oligo_1'], data={
                'oligo_aa_sequence': ''.join(np.random.choice(AMINO_INFO['amino_acid'].unique(), 300))})
            oligo_rows['nuc_sequence'] = oligo_rows.oligo_aa_sequence.apply(code_one_aa_sequence_to_nuc)
            start_location = 0
            for i, barcode_length in enumerate(BARCODE_NUC_LENGTHS):
                oligo_rows[f'barcode_{i}'] = oligo_rows['nuc_sequence'].apply(
                    lambda nuc: get_barcode_from_nuc_seq(nuc, start_location, barcode_length))
                start_location += barcode_length
            ret, _ = create_new_nuc_sequence(oligo_rows.iloc[0].copy(), oligo_rows)
            barcoded_length = 3 * math.ceil(sum(BARCODE_NUC_LENGTHS) / 3)
            self.assertEqual(ret['nuc_sequence'][barcoded_length:],
                             oligo_rows.iloc[0]['nuc_sequence'][barcoded_length:])
            self.assertNotEqual(ret['nuc_sequence'],
                                oligo_rows.iloc[0]['nuc_sequence'])
            for i in range(len(BARCODE_NUC_LENGTHS)):
                self.assertNotEqual(ret[f'barcode_{i}'], oligo_rows.iloc[0][f'barcode_{i}'])

    def test_get_barcode_from_nuc_seq(self):
        nuc_seq = ''.join(np.random.choice(['A', 'G', 'C', 'T'], 50))
        with mock.patch('BuildPhIPSeqLibrary.construct_nucleotide_sequences.BARCODE_IN_5_PRIME_END', True):
            barcode = get_barcode_from_nuc_seq(nuc_seq, 0, 15)
            self.assertEqual(barcode, nuc_seq[:15])
            barcode = get_barcode_from_nuc_seq(nuc_seq, 20, 25)
            self.assertEqual(barcode, nuc_seq[20:20 + 25])
        with mock.patch('BuildPhIPSeqLibrary.construct_nucleotide_sequences.BARCODE_IN_5_PRIME_END', False):
            barcode = get_barcode_from_nuc_seq(nuc_seq, 0, 15)
            self.assertEqual(barcode, nuc_seq[-15:])
            barcode = get_barcode_from_nuc_seq(nuc_seq, 20, 25)
            self.assertEqual(barcode, nuc_seq[-20 - 25:-20])

    @mock.patch('BuildPhIPSeqLibrary.read_pipeline_files.BARCODED_NUC_FILE',
                os.path.join(MOCK_DATA_DIR, 'Output', 'barcoded_nuc_file.csv'))
    @mock.patch('BuildPhIPSeqLibrary.read_pipeline_files.UNCONVERTED_SEQUENCES_FILE',
                os.path.join(MOCK_DATA_DIR, 'Output', 'unconverted_sequences.csv'))
    @mock.patch('pandas.DataFrame.to_csv', print)
    def test_barcode_sequences(self):
        with mock.patch('BuildPhIPSeqLibrary.construct_nucleotide_sequences.BARCODE_IN_5_PRIME_END', True):
            barcode_aa_length = math.ceil(sum(BARCODE_NUC_LENGTHS) / 3)
            oligo_sequences = pd.Series(
                {'aa_sequence_1': AMINO_INFO.groupby('amino_acid')['codon'].count().sort_values().index[
                                      0] * barcode_aa_length + ''.join(
                    np.random.choice(AMINO_INFO['amino_acid'].unique(), OLIGO_AA_LENGTH - barcode_aa_length)),
                 'aa_sequence_2': AMINO_INFO.groupby('amino_acid')['codon'].count().sort_values().index[
                                      0] * barcode_aa_length + ''.join(
                     np.random.choice(AMINO_INFO['amino_acid'].unique(),
                                      OLIGO_AA_LENGTH - barcode_aa_length))}).to_frame().rename(
                columns={0: 'oligo_aa_sequence'})
            oligo_sequences['nuc_sequence'] = oligo_sequences.oligo_aa_sequence.apply(code_one_aa_sequence_to_nuc)
            ret, unconverted = barcode_sequences(oligo_sequences)
            self.assertEqual(len(ret), 1)
            self.assertEqual(len(unconverted), 1)
            barcode_nuc_lengths = [9]

            with mock.patch('BuildPhIPSeqLibrary.construct_nucleotide_sequences.BARCODE_NUC_LENGTHS',
                            barcode_nuc_lengths):
                with mock.patch('BuildPhIPSeqLibrary.read_pipeline_files.BARCODE_NUC_LENGTHS', barcode_nuc_lengths):
                    for _ in range(200):
                        num_repetitions = AMINO_INFO['amino_acid'].eq('L').sum()
                        aa_seq = 'L' * 3
                        num_aa_in_barcode = math.ceil(barcode_nuc_lengths[0] / 3)
                        oligos = pd.DataFrame(
                            index=[f'oligo_{i}' for i in range((num_repetitions ** num_aa_in_barcode) + 1)],
                            data={'oligo_aa_sequence': [aa_seq] * ((num_repetitions ** num_aa_in_barcode) + 1)})
                        oligos['nuc_sequence'] = oligos.oligo_aa_sequence.apply(code_one_aa_sequence_to_nuc)
                        ret, unconverted = barcode_sequences(oligos)
                        self.assertEqual(len(unconverted), 5)
                        self.assertEqual(len(ret) + len(unconverted), (num_repetitions ** num_aa_in_barcode) + 1)

    def test_iterative_barcode_construction(self):
        for barcode_nuc_lengths in [[3, 6], [4, 5], [3, 5]]:
            barcode_size_in_aa = math.ceil(sum(barcode_nuc_lengths) / 3)
            with mock.patch('BuildPhIPSeqLibrary.construct_nucleotide_sequences.BARCODE_NUC_LENGTHS',
                            barcode_nuc_lengths):
                with mock.patch('BuildPhIPSeqLibrary.construct_nucleotide_sequences.BARCODE_IN_5_PRIME_END', True):
                    aa_seq = AMINO_INFO.groupby('amino_acid')['codon'].count().sort_values().index[0] * 6
                    nuc_seq = code_one_aa_sequence_to_nuc(aa_seq)
                    existing_barcodes = pd.DataFrame(columns=['oligo_id', 'nuc_sequence'] + list(
                        map(lambda i: f"barcode_{i}", range(len(barcode_nuc_lengths)))),
                                                     data=[['oligo_1', nuc_seq, nuc_seq[:barcode_nuc_lengths[0]],
                                                            nuc_seq[barcode_nuc_lengths[0]:sum(
                                                                barcode_nuc_lengths)]]]).set_index(
                        'oligo_id')
                    ret = iterative_barcode_construction(aa_seq[:barcode_size_in_aa], '',
                                                         nuc_seq[3 * barcode_size_in_aa:],
                                                         existing_barcodes)
                    self.assertIsNone(ret)

                    aa_seq = 'F' + aa_seq[1:]
                    nuc_seq = code_one_aa_sequence_to_nuc(aa_seq)
                    ret = iterative_barcode_construction(aa_seq[:barcode_size_in_aa], '',
                                                         nuc_seq[3 * barcode_size_in_aa:],
                                                         existing_barcodes)
                    self.assertIsNone(ret)

                    aa_seq = AMINO_INFO.groupby('amino_acid')['codon'].count().sort_values().index[-1] * 6
                    nuc_seq = code_one_aa_sequence_to_nuc(aa_seq)
                    existing_barcodes = pd.DataFrame(columns=['oligo_id', 'nuc_sequence'] + list(
                        map(lambda i: f"barcode_{i}", range(len(barcode_nuc_lengths)))),
                                                     data=[['oligo_1', nuc_seq, nuc_seq[:barcode_nuc_lengths[0]],
                                                            nuc_seq[barcode_nuc_lengths[0]:sum(
                                                                barcode_nuc_lengths)]]]).set_index(
                        'oligo_id')
                    ret = iterative_barcode_construction(aa_seq[:barcode_size_in_aa], '',
                                                         nuc_seq[3 * barcode_size_in_aa:],
                                                         existing_barcodes)
                    self.assertIsNotNone(ret)
                    self.assertEqual(len(ret), len(nuc_seq))

                with mock.patch('BuildPhIPSeqLibrary.construct_nucleotide_sequences.BARCODE_IN_5_PRIME_END', False):
                    aa_seq = AMINO_INFO.groupby('amino_acid')['codon'].count().sort_values().index[0] * 6
                    nuc_seq = code_one_aa_sequence_to_nuc(aa_seq)
                    existing_barcodes = pd.DataFrame(columns=['oligo_id', 'nuc_sequence'] + list(
                        map(lambda i: f"barcode_{i}", range(len([3, 6])))),
                                                     data=[['oligo_1', nuc_seq, nuc_seq[-barcode_nuc_lengths[0]:],
                                                            nuc_seq[-(sum(barcode_nuc_lengths)):-barcode_nuc_lengths[
                                                                0]]]]).set_index(
                        'oligo_id')
                    ret = iterative_barcode_construction(aa_seq[-barcode_size_in_aa:],
                                                         nuc_seq[:-(3 * barcode_size_in_aa)], '',
                                                         existing_barcodes)
                    self.assertIsNone(ret)

                    aa_seq = aa_seq[:-1] + 'F'
                    nuc_seq = code_one_aa_sequence_to_nuc(aa_seq)
                    ret = iterative_barcode_construction(aa_seq[-barcode_size_in_aa:],
                                                         nuc_seq[:-(3 * barcode_size_in_aa)], '',
                                                         existing_barcodes)
                    self.assertIsNone(ret)

                    aa_seq = AMINO_INFO.groupby('amino_acid')['codon'].count().sort_values().index[-1] * 6
                    nuc_seq = code_one_aa_sequence_to_nuc(aa_seq)
                    existing_barcodes = pd.DataFrame(columns=['oligo_id', 'nuc_sequence'] + list(
                        map(lambda i: f"barcode_{i}", range(len([3, 6])))),
                                                     data=[['oligo_1', nuc_seq, nuc_seq[-barcode_nuc_lengths[0]:],
                                                            nuc_seq[-(sum(barcode_nuc_lengths)):-barcode_nuc_lengths[
                                                                0]]]]).set_index(
                        'oligo_id')
                    ret = iterative_barcode_construction(aa_seq[-barcode_size_in_aa:],
                                                         nuc_seq[:-(3 * barcode_size_in_aa)], '',
                                                         existing_barcodes)
                    self.assertIsNotNone(ret)
                    self.assertEqual(len(ret), len(nuc_seq))

    def test_iterative_correction_of_single_barcode(self):
        for barcode_nuc_lengths in [[3, 6, 3], [4, 5, 3]]:
            with mock.patch('BuildPhIPSeqLibrary.construct_nucleotide_sequences.BARCODE_NUC_LENGTHS',
                            barcode_nuc_lengths):
                for barcode_at_5_prime in [True, False]:
                    with mock.patch('BuildPhIPSeqLibrary.construct_nucleotide_sequences.BARCODE_IN_5_PRIME_END',
                                    barcode_at_5_prime):
                        aa_seq = AMINO_INFO.groupby('amino_acid')['codon'].count().sort_values().index[
                                     -1] * OLIGO_AA_LENGTH
                        nuc_seq = code_one_aa_sequence_to_nuc(aa_seq)
                        oligo_row = {'oligo_id': 'oligo_1', 'nuc_sequence': nuc_seq, 'oligo_aa_sequence': aa_seq}
                        oligo_row = get_all_barcodes(oligo_row)
                        existing_barcodes = pd.Series(oligo_row).to_frame().T.set_index('oligo_id')
                        nuc_seq = oligo_row['barcode_0']
                        if barcode_at_5_prime:
                            nuc_seq += code_one_aa_sequence_to_nuc(aa_seq)[barcode_nuc_lengths[0]:]
                        else:
                            nuc_seq = code_one_aa_sequence_to_nuc(aa_seq)[:-barcode_nuc_lengths[0]] + nuc_seq
                        oligo_row['oligo_id'] = 'oligo_2'
                        ret = iterative_correction_of_single_barcode([nuc_seq], oligo_row, existing_barcodes)
                        length_of_changed_barcode = 3 * math.ceil(barcode_nuc_lengths[0] / 3)
                        self.assertEqual(
                            get_subset_of_sequence(nuc_seq, length_of_changed_barcode, len(nuc_seq),
                                                   barcode_at_5_prime),
                            get_subset_of_sequence(ret['nuc_sequence'], length_of_changed_barcode, len(nuc_seq),
                                                   barcode_at_5_prime))
                        self.assertNotEqual(
                            get_subset_of_sequence(nuc_seq, 0, length_of_changed_barcode, barcode_at_5_prime),
                            get_subset_of_sequence(ret['nuc_sequence'], 0, length_of_changed_barcode,
                                                   barcode_at_5_prime))
                        self.assertNotIn(oligo_row['barcode_0'], existing_barcodes['barcode_0'].values)

    @mock.patch('BuildPhIPSeqLibrary.read_pipeline_files.OLIGO_SEQUENCES_FILE',
                os.path.join(MOCK_DATA_DIR, 'PipelineFiles', 'oligos_sequence.csv'))
    @mock.patch('BuildPhIPSeqLibrary.read_pipeline_files.BARCODED_NUC_FILE',
                os.path.join(MOCK_DATA_DIR, 'PipelineFiles', 'barcoded_nuc_file.csv'))
    @mock.patch('BuildPhIPSeqLibrary.construct_nucleotide_sequences.update_unconverted_oligos_file', print)
    @mock.patch('pandas.DataFrame.to_csv', print)
    @mock.patch('multiprocessing.Pool', Pool)
    def test_aa_to_nuc(self):
        aa_oligos = read_oligo_sequences_to_file()
        ret = aa_to_nuc(aa_oligos)
        self.assertEqual(len(aa_oligos), len(ret))


def get_subset_of_sequence(seq, start, end, not_reverse):
    if not_reverse:
        return seq[start: end]
    else:
        return seq[::-1][start: end][::-1]
