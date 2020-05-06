# stabilityMotifs

The python (version 2.7) script analyzedMotifsEnrichments.py compares the occurrences of several
stability motifs within a set of tested sequences to that obtained by selecting an equally
populated set of random sequences from the human Hg38 database

The software depends on biopython and scipy which can be installed with the following commands

pip install biopython 
pip install scipy

The script can be used by typing the following example command:

python analyzedStabiltyFactorsRandom.py cds_test.fasta 3utr_test.fasta 5utr_test.fasta 1000

where cds_test.fasta, 3utr_test.fasta and 5utr_test.fasta are the coding sequences, 3'UTR 
and 5'UTR sequences respectively of the transcript dataset one wants to test the stability
motif on. The last number (e.g. 5 in the example command) is the number of randomizations
to perform. 

The script returns the following output file:

randomizationResults.txt

This file reports for each transcript portion (cds, 3UTR, 5UTR or cDNA) the number of times
the number of found motifs is higher or lower with the respect to the number of time the
same motif was found in the random dataset.



