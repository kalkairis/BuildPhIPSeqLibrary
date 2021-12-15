import argparse
import os.path

from BuildPhIPSeqLibrary.config import OUTPUT_DIR
from BuildPhIPSeqLibrary.read_input_files import get_input_files, read_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Building new library main.")
    parser.add_argument('--overwrite', type=bool, default=False, action='store_true')
    args = parser.parse_args()
    if args['overwrite'] and os.path.exists(OUTPUT_DIR):
        assert len([filename for filename in os.listdir(OUTPUT_DIR) if
                    filename != 'README.md']), f"In order to overwrite you must empty the output dir {OUTPUT_DIR}"

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    files = get_input_files()
    for filename in files:
        seq_id_to_sequences = read_file(filename)
        seq_id_to_sequences = add_sequences_to_files_list(seq_id_to_sequences, filename)


