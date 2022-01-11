import argparse
import os.path

from BuildPhIPSeqLibrary.config import OUTPUT_DIR
from BuildPhIPSeqLibrary.construct_nucleotide_sequences import aa_to_nuc
from BuildPhIPSeqLibrary.read_input_files import get_input_files, read_file
from BuildPhIPSeqLibrary.sequence_ids import add_sequences_to_files_list
from BuildPhIPSeqLibrary.split_sequences_to_oligos import split_sequences_to_oligos, split_and_map_new_sequences

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Building new library main.")
    parser.add_argument('--overwrite', dest='overwrite', default=False, action='store_true')
    args = parser.parse_args()
    if args.overwrite and os.path.exists(OUTPUT_DIR):
        # TODO: remove line
        for filename in os.listdir(OUTPUT_DIR):
            if filename != 'README.md':
                os.remove(os.path.join(OUTPUT_DIR, filename))
        assert len([filename for filename in os.listdir(OUTPUT_DIR) if
                    filename != 'README.md']) == 0, f"In order to overwrite you must empty the output dir {OUTPUT_DIR}"

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    files = get_input_files()
    for filename in files:
        seq_id_to_sequences = read_file(filename)
        seq_id_to_sequences = add_sequences_to_files_list(seq_id_to_sequences, filename)
        all_oligos_aa_sequences, new_oligos_aa_sequences = split_and_map_new_sequences(seq_id_to_sequences)
        oligo_non_barcoded_nuc_sequences = aa_to_nuc(new_oligos_aa_sequences)
        print("here")
