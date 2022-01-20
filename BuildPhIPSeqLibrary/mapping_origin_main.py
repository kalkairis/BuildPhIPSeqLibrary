"""
This code maps the origin as well as the possible mappings of each oligo.
Note that this is NOT part of the main library building code.
"""
import pandas
import re
import time
from multiprocessing import Pool

from BuildPhIPSeqLibrary.read_pipeline_files import read_oligo_sequences_to_file, read_sequence_ids_file
from BuildPhIPSeqLibrary.config import NUM_MAPPING_THREADS, MAPPED_OLIGO_SEQUENCES_FILE


def map_single_oligo_to_sequences_list(inp):
    id, oligo, sequences_dict = inp
    ret = []
    if '*' in oligo:
        oligo = oligo[:oligo.index('*')]
    for sequence_id, sequence in sequences_dict.items():
        for position in re.finditer(oligo, sequence):
            ret.append((sequence_id, position.span()[0]))
    return [id, ret]


if __name__ == "__main__":
    p = Pool(NUM_MAPPING_THREADS)

    aa_seqs = read_sequence_ids_file().AA_sequence.to_dict()
    aa_olis = read_oligo_sequences_to_file()
    lst = []
    for id in aa_olis.index:
        lst.append((id, aa_olis.loc[id].oligo_aa_sequence, aa_seqs))
    print("Running maps on %d oligos, with %d threads" % (len(lst), NUM_MAPPING_THREADS), time.ctime())
    maps = p.map(map_single_oligo_to_sequences_list, lst)
    aa_olis['mapped'] = pandas.DataFrame(maps, columns=['oligo_id', 'mapped']).set_index('oligo_id')
    aa_olis.to_csv(MAPPED_OLIGO_SEQUENCES_FILE)
    print("Mapping done on %d oligos" % len(lst), time.ctime())
    p.close()
