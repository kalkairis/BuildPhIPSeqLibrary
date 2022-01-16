from BuildPhIPSeqLibrary.config import PREFIX_PROMOTER, SUFFIX_PROMOTER, ORDER_FILE


def transfer_to_order(oligo_barcoded_sequences):
    nuc_sequences = '\n'.join((PREFIX_PROMOTER + oligo_barcoded_sequences['nuc_sequence'] + SUFFIX_PROMOTER).values)
    with open(ORDER_FILE, 'w') as file:
        file.write(nuc_sequences)
