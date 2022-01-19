from BuildPhIPSeqLibrary.config import PREFIX, SUFFIX ORDER_FILE


def transfer_to_order(oligo_barcoded_sequences):
    nuc_sequences = '\n'.join((PREFIX_PROMOTER + oligo_barcoded_sequences['nuc_sequence'] + SUFFIX_PROMOTER).values)
    with open(ORDER_FILE, 'w') as file:
        file.write(nuc_sequences)
