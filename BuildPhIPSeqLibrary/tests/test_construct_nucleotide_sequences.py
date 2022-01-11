from unittest import TestCase, mock

import numpy as np

from BuildPhIPSeqLibrary.config import RESTRICTED_SEQUENCES, AMINO_ACIDS, AMINO_INFO
from BuildPhIPSeqLibrary.construct_nucleotide_sequences import has_no_restricted_sequences, code_one_aa_sequence_to_nuc, \
    get_barcode_from_nuc_seq


class Test(TestCase):
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
        self.fail()

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

    def test_barcode_sequences(self):
        self.fail()

    def test_aa_to_nuc(self):
        self.fail()
