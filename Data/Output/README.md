# Output directory

In this directory all output files will be added during the run. This include:

- `files_hash.csv`: A hash table to ensure files have not been changed. Changes should be made by newly added files or
  by rerunning the code from scratch.
- `sequences_ids.csv`: A unique running index identifier based on origin file and the given sequence id.
- `oligos_sequence.csv`: A table were each nucleotide oligo has its origin sequences (this can be more than one and is
  therefore a list of sequence IDs and positions), and a list of mapped locations (again, can be more than one).
- `barcoded_nuc_file.csv`: A table of oligo_id, nucleotide sequence, and barcodes of the specific oligo. Oligos are
  uniquely barcoded for each one of the barcodes, allowing them to be identified by multiple sequences.
- `unconverted_sequences.csv`: This is the list of oligos which have failed to be barcoded. Reasons vary from inability
  to create the sequence without restriction sites, to having no viable unique barcode for the specific oligo.

  __Note__: if this file does not exist, it means that all oligos have successfully been added to the barcodes file.
  This is not likely in a large library but in smaller ones is a possibility. 
