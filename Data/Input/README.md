# Input directory

In this directory all input files should be added. It is possible to add new files and re-run the entire script.
However, it is not possible to change existing file without erasing the entire output directory first.

File format should be either: 

.csv files, where the following columns MUST exist
- `sequence_ID` (must be unique within file)
- `AA_sequence` the entire AA sequence to add into the library. This must include string sequences from the 20 AA
  alphabet. No stop codons or unknown AAs (lines with unknown symbols will be ignored, with a warning message).   

.fa files where
- record name will be treated as 'sequence_ID' 
- record sequence will be treated as 'AA_sequence'
