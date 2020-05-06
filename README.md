# stabilityMotifs

The python (version 2.7) script analyzedStabiltyFactorsRandom.py compare the occurrences of several
stability motifs with a set of sequences to that obtained by randomly selecting an equally
populated set of sequences from the human Hg38 database

The software depends on biopython and scipy that can be installed by the following commands

pip install biopython 
pip install scipy

The script can be used by typing the following example command:

python analyzedStabiltyFactorsRandom.py cds_test.fasta 3utr_test.fasta 5utr_test.fasta 5

where cds_test.fasta, 3utr_test.fasta and 5utr_test.fasta are the coding sequences, 3'UTR 
and 5'UTR sequences respectively of the transcript dataset one wants to test the stability
motif on. The last number (e.g. 5 in the example command) is the number of randomizations
to perform. 

The script returns the following output files:

randomizationResults.txt
this reports, for each randomization step, the transcripts that have been randomly selected 
together with the counts of all the analyzed stability motifs

randomizationAverageResults.txt
this reports, for each randomization step, the average value of the counts of all the 
analyzed motifs

anovaResults.txt
This reports, for each motif, the number of time the test dataset showed, on average, a higher
number of occurrences, the number time this difference was significant, the number of
times the test dataset showed a lower number of occurrences and the number ot time this
difference was significant (significancy tested by t-test)

