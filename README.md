# blast-common-ancestor

A python script for assigning taxonomy to sequences in fasta files using E-utilities and command-line blast.

Takes the top 10 blastn hits for each sequence, fetches their full lineage with E-utilities, and finds the smallest grouping that each hit shares. Outputs a corresponding .tax file for each fasta file with sequence name and taxonomic assignment.

To run use:
```
python blastcommonancestor.py [path to fasta or directory of fastas]
```
