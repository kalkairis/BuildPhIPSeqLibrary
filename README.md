# BuildPhIPSeqLibrary

Building a new Phage IP Sequencing library.

Input files are described in:
BuildPhIPSeqLibrary/Data/Input/README.md

Output files are described in:
BuildPhIPSeqLibrary/Data/Output/README.md

Run process:
1. Run tests, it is always better to make sure things run correctly on a new platform.
2. Go over config.py and change any needed parameters.
3. Create your formatted input files.
4. Run BuildPhIPSeqLibrary/main.py to create most output files.
5. Run BuildPhIPSeqLibrary/mapping_origin_main.py (can be done offline, takes time to run) to create the 
   mapped_oligos_sequence output file.
