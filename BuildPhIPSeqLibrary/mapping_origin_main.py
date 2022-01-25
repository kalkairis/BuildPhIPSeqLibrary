"""
This code maps the origin as well as the possible mappings of each oligo.
Note that this is NOT part of the main library building code.
"""
import time
from multiprocessing import Pool

import pandas
import regex as re

from BuildPhIPSeqLibrary.config import NUM_MAPPING_THREADS, MAPPED_OLIGO_SEQUENCES_FILE
from BuildPhIPSeqLibrary.read_pipeline_files import read_oligo_sequences_to_file, read_sequence_ids_file


def map_single_oligo_to_sequences_list(oligo_id, oligo, sequences_dict):
    """
    Maps a single oligo to all sequences in input.
    :param oligo_id:
    :param oligo: An amino sequence of the oligo
    :param sequences_dict: A dictionary from seq_ID to amino acid sequences.
    :return:
    """
    ret = []
    if '*' in oligo:
        oligo = oligo[:oligo.index('*')]
    for sequence_id, sequence in sequences_dict.items():
        for position in re.finditer(oligo, sequence, overlapped=True):
            ret.append((sequence_id, position.span()[0]))
    return [oligo_id, ret]


def run_mapping_of_all_files():
    pool = Pool(NUM_MAPPING_THREADS)

    aa_seqs = read_sequence_ids_file().AA_sequence.to_dict()
    aa_olis = read_oligo_sequences_to_file()
    lst = []
    for oligo_id in aa_olis.index:
        lst.append((oligo_id, aa_olis.loc[oligo_id].oligo_aa_sequence, aa_seqs))
    print("Running maps on %d oligos, with %d threads" % (len(lst), NUM_MAPPING_THREADS), time.ctime())
    maps = pool.starmap(map_single_oligo_to_sequences_list, lst)
    aa_olis['mapped'] = pandas.DataFrame(maps, columns=['oligo_id', 'mapped']).set_index('oligo_id')
    aa_olis.to_csv(MAPPED_OLIGO_SEQUENCES_FILE)
    print("Mapping done on %d oligos" % len(lst), time.ctime())
    pool.close()
    return aa_olis


if __name__ == "__main__":
    run_mapping_of_all_files()
