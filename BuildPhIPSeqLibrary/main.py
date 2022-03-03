import argparse
import logging
import os.path
import time

from BuildPhIPSeqLibrary.config import OUTPUT_DIR, BARCODED_NUC_FILE, ORDER_FILE
from BuildPhIPSeqLibrary.construct_nucleotide_sequences import aa_to_nuc
from BuildPhIPSeqLibrary.output_to_order import transfer_to_order
from BuildPhIPSeqLibrary.read_input_files import get_input_files, read_file
from BuildPhIPSeqLibrary.sequence_ids import add_sequences_to_files_list
from BuildPhIPSeqLibrary.split_sequences_to_oligos import split_and_map_new_sequences

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Building new library main.")
    parser.add_argument('--overwrite', dest='overwrite', default=False, action='store_true')
    args = parser.parse_args()
    if args.overwrite and os.path.exists(OUTPUT_DIR):
        assert len([filename for filename in os.listdir(OUTPUT_DIR) if
                    filename != 'README.md' and not filename.startswith(
                        '.')]) == 0, f"""In order to overwrite you must empty the output dir {OUTPUT_DIR}
Consider running:
for filename in os.listdir('{OUTPUT_DIR}'):
    if filename != 'README.md':
        os.remove(os.path.join('{OUTPUT_DIR}', filename))""""""
        """
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    files = get_input_files()
    logging.basicConfig(level=logging.INFO)
    for filename in files:
        logging.info(f"Working on file {filename}, {time.ctime()}")
        logging.info(f"{filename}: reading file, {time.ctime()}")
        seq_id_to_sequences = read_file(filename)
        logging.info(f"{filename}: Got {len(seq_id_to_sequences)} sequences, {time.ctime()}")
        logging.info(f"{filename}: Identifying new sequences, {time.ctime()}")
        seq_id_to_sequences = add_sequences_to_files_list(seq_id_to_sequences, filename)
        logging.info(f"{filename}: Got {len(seq_id_to_sequences)} new sequences, {time.ctime()}")
        logging.info(f"{filename}: Converting sequences to oligos, {time.ctime()}")
        all_oligos_aa_sequences, new_oligos_aa_sequences = split_and_map_new_sequences(seq_id_to_sequences)
        logging.info(f"{filename}: "
              f"Converted {len(new_oligos_aa_sequences)} new oligos "
              f"overall {len(all_oligos_aa_sequences)} oligos, {time.ctime()}")
        logging.info(f"{filename}: Barcoding oligos, {time.ctime()}")
        oligo_barcoded_sequences = aa_to_nuc(new_oligos_aa_sequences)
        logging.info(
            f"{filename}: Finished barcoding oligos. "
            f"Current number of oligos is {len(oligo_barcoded_sequences)}, {time.ctime()}")
    if len(files) > 0:
        logging.info(f"Finished creating sequences. Find them in {BARCODED_NUC_FILE}, {time.ctime()}")
        logging.info(f"Converting all sequences to order, {time.ctime()}")
        transfer_to_order(oligo_barcoded_sequences)
        logging.info(f"Find file of order in {ORDER_FILE}, {time.ctime()}")
    else:
        logging.info("No new files")
