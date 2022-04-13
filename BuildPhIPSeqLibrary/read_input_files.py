import glob
import hashlib
import logging
import os.path

import pandas
import pandas as pd
from Bio import SeqIO

from BuildPhIPSeqLibrary.config import INPUT_DIR, seq_ID_col, seq_AA_col, FILES_INPUT_HASH_FILE


def get_input_files(files_hash_path=None, **kwargs):
    """
    Lists files to process for library construction.
    Asserts no file has been changed.
    Will warn if file has been removed from INPUT_DIR.
    Saves MD5 hash of all files in OUTPUT_DIR/files_hash.csv
    :param files_hash_path:
    :param kwargs: unused
    :return: List of new files to process
    """
    if files_hash_path is None:
        files_hash_path = FILES_INPUT_HASH_FILE
    input_files = glob.glob(os.path.join(INPUT_DIR, "*.csv")) + glob.glob(os.path.join(INPUT_DIR, "*.fa"))
    new_added_files = set(input_files)
    input_hashes = {}
    for filename in input_files:
        hash_md5 = hashlib.md5()
        with open(filename, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        input_hashes[os.path.basename(filename)] = hash_md5.hexdigest()

    if os.path.exists(files_hash_path):
        files_hash = pd.read_csv(files_hash_path, index_col=0)['0'].to_dict()
        for filename, file_hash in files_hash.items():
            if os.path.join(INPUT_DIR, filename) not in input_files:
                logging.debug(
                    f"File {filename} exists in previous version of the library construction,"
                    f" but is absent from {INPUT_DIR}.")
            else:
                assert input_hashes[filename] == file_hash, f"File {filename} changed content. " \
                                                            f"Re-run library construction from scratch."
                new_added_files.remove(os.path.join(INPUT_DIR, filename))
        input_hashes.update(files_hash)
    pd.Series(input_hashes).to_csv(files_hash_path, header=True)
    if len(new_added_files) == 0:
        logging.warning("No new files are added.")
    return new_added_files


def read_file(file_path):
    """
    Reads a single CSV file to process. Ensuring columns 'sequence_ID', and 'AA_sequence' are in the file.
    Ensures 'sequence_ID' is a unique key.
    :param file_path:
    :return: Dict from 'sequence_ID' to 'AA_sequence'.
    """
    if os.path.splitext(file_path)[1] == '.csv':
        ret = pd.read_csv(file_path, usecols=[seq_ID_col, seq_AA_col])
    else:
        ret = []
        for rec in SeqIO.parse(file_path, 'fasta'):
            ret.append([rec.name, str(rec.seq)])
        ret = pandas.DataFrame(ret, columns=[seq_ID_col, seq_AA_col])
    value_counts = ret[seq_ID_col].value_counts()
    value_counts = value_counts[value_counts.gt(1)]
    error_values = '\n'.join(value_counts.index.values)
    assert len(value_counts) == 0, IOError(
        f"Repeating sequence_ID in {file_path}, for sequences: {error_values}")
    return ret.set_index(seq_ID_col)[seq_AA_col].to_dict()
