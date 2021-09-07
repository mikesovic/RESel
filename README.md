# RESel

requires...
fasta-splitter.pl (script from http://kirill-kryukov.com/study/tools/fasta-splitter/)
format_fasta.pl
digest.sh
digest_single.sh
analyze_dig.R
analyze_dig_single.R
parallel_digest.sh
parallel_digest_single.sh

paths for all of the above are hard-coded in RESel and parallel_RESel scripts; 

Also need GNU grep as default grep on system (BSD grep - likely default on Macs, won't work)

Usage...

There are two version of RESel - a standard version that can be run on a personal computer, and a parallel version that is set up for HPC.

Each require...
fasta genome file
Restriction enzyme file

Each allow for optional features file in bed format

The standard version runs all pairwise combinations of restriction enzymes defined in the restriction enzyme file in succession.
The parallel version sends a separate job for each enzyme pair.

After runs for all combinations are complete, can run...
RESel_cleanup.sh *fasta base name* (fasta file without extension)

This produces *summary* files in the '_out' dir.



