from Bio import SeqIO
import sys
import math
import random as rd
import numpy as nm
from scipy import stats

def occurrences(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count

queryCDS = SeqIO.to_dict(SeqIO.parse(sys.argv[1],"fasta"))
query3UTR = SeqIO.to_dict(SeqIO.parse(sys.argv[2],"fasta"))
query5UTR = SeqIO.to_dict(SeqIO.parse(sys.argv[3],"fasta"))
randomArg = int(sys.argv[4])
querycDNA = {}

allCDS = SeqIO.to_dict(SeqIO.parse("cds_lv.fasta","fasta"))
all3UTR = SeqIO.to_dict(SeqIO.parse("3utr_lv.fasta","fasta"))
all5UTR = SeqIO.to_dict(SeqIO.parse("5utr_lv.fasta","fasta"))
allcDNA = {}


m_2d_motifs_base = 'ATTTATTTA'
m_2d_motifs = []
for a in ['A','T']:
    for b in ['A','T']:
        for c in ['A','T']:
            for d in ['A','T']:
                m_2d_motifs.append(a+b+m_2d_motifs_base+c+d)

m_2e_motifs_base = 'ATTTATTTA'
m_2e_motifs = []
for a in ['A','T']:
    for b in ['A','T']:
        for c in ['A','T']:
            for d in ['A','T']:
                for e in ['A','T']:
                    for f in ['A','T']:
                        for g in ['A','T']:
                            for h in ['A','T']:
                                m_2e_motifs.append(a+b+c+d+m_2d_motifs_base+e+f+g+h)

motifsList = [['AATTTAA','TATTTAT','AATTTAT','TATTTAA'],['ATTTATTTATTTATTTATTTA'],['ATTTATTTATTTATTTA'],['AATTTATTTATTTAA','TATTTATTTATTTAT','AATTTATTTATTTAT','TATTTATTTATTTAA'],m_2d_motifs,m_2e_motifs,["TTATTTATT"],["TATTTAT"],["TTTTTTT"],["TTTTTAAA"],["TTGTAAATA"],["TTTTAAAT"],["TTTTAATTT"],["AAATATTTT"],["AATATTTTT"],["CCGCCTC"],["CCAGCCTC"],["GGGCCTGG"],["CCCAGCCCC"]]
motifsNames = ["1","2a","2b","2c","2d","2e","meg","megshort","h1",'h2',"h3","b1","b2","b3","b4","h_1","h_2","b_1","b_2"]


#Calculate the initial values
initial_5utr = {}
initial_cds = {}
initial_3utr = {}
initial_cDNA = {}

for a in motifsNames:
    if not a in initial_3utr:
        initial_3utr[a] = 0
        initial_cds[a] = 0
        initial_5utr[a] = 0
        initial_cDNA[a] = 0
    

numCDS = 0
num5UTR = 0
num3UTR = 0

for transcript in queryCDS:
    numCDS += 1
    if not transcript in querycDNA:
        querycDNA[transcript] = str(queryCDS[transcript].seq)
    if transcript in query5UTR:
        num5UTR+=1
        querycDNA[transcript] = str(query5UTR[transcript].seq)+querycDNA[transcript]
    if transcript in query3UTR:
        num3UTR+=1
        querycDNA[transcript] = querycDNA[transcript]+str(query3UTR[transcript].seq)

#Calculate occurrences
for a in range(len(motifsNames)):
    for motif in motifsList[a]:
        for transcript in query3UTR:
            initial_3utr[motifsNames[a]] += float(occurrences(str(query3UTR[transcript].seq),motif))
        for transcript in query5UTR:
            initial_5utr[motifsNames[a]] += float(occurrences(str(query5UTR[transcript].seq),motif))
        for transcript in queryCDS:
            initial_cds[motifsNames[a]] += float(occurrences(str(queryCDS[transcript].seq),motif))
        for transcript in querycDNA:
            initial_cDNA[motifsNames[a]] += float(occurrences(str(querycDNA[transcript]),motif))



#Calculate averages
for a in range(len(motifsNames)):
    initial_3utr[motifsNames[a]] = initial_3utr[motifsNames[a]] / num3UTR
    initial_5utr[motifsNames[a]] = initial_5utr[motifsNames[a]] / num5UTR
    initial_cds[motifsNames[a]] = initial_cds[motifsNames[a]] / numCDS
    initial_cDNA[motifsNames[a]] = initial_cDNA[motifsNames[a]] / numCDS




numHigher_3UTR = {}
numLower_3UTR = {}
numHigher_5UTR = {}
numLower_5UTR = {}
numHigher_CDS = {}
numLower_CDS = {}
numHigher_cDNA = {}
numLower_cDNA = {}

for a in range(len(motifsNames)):
    if not motifsNames[a] in numHigher_3UTR:
        numHigher_3UTR[motifsNames[a]] = 0
        numLower_3UTR[motifsNames[a]] = 0
        numHigher_5UTR[motifsNames[a]] = 0
        numLower_5UTR[motifsNames[a]] = 0
        numHigher_CDS[motifsNames[a]] = 0
        numLower_CDS[motifsNames[a]] = 0
        numHigher_cDNA[motifsNames[a]] = 0
        numLower_cDNA[motifsNames[a]] = 0



for numRandomization in range(randomArg):
    #Calculate the random values
    random_5utr = {}
    random_cds = {}
    random_3utr = {}
    random_cDNA = {}
    


    for a in motifsNames:
        if not a in random_3utr:
            random_3utr[a] = 0
            random_cds[a] = 0
            random_5utr[a] = 0
            random_cDNA[a] = 0
    
    print "Randomization",numRandomization
    foundGenes = 0
    all_cDNA = {}

    while foundGenes <numCDS:
        randomTranscript = rd.choice(allCDS.keys())
        if randomTranscript in allCDS and randomTranscript in all3UTR and randomTranscript in all5UTR:
            foundGenes+=1
        
            if not randomTranscript in all_cDNA:
                all_cDNA[randomTranscript] = str(all5UTR[randomTranscript].seq)+str(allCDS[randomTranscript].seq)+str(all3UTR[randomTranscript].seq)


            #Calculate occurrences
            for a in range(len(motifsNames)):
                for motif in motifsList[a]:
                    random_3utr[motifsNames[a]] += float(occurrences(str(all3UTR[randomTranscript].seq),motif))
                    random_5utr[motifsNames[a]] += float(occurrences(str(all5UTR[randomTranscript].seq),motif))
                    random_cds[motifsNames[a]] += float(occurrences(str(allCDS[randomTranscript].seq),motif))
                    random_cDNA[motifsNames[a]] += float(occurrences(str(all_cDNA[randomTranscript]),motif))

     

    #Calculate averages
    for a in range(len(motifsNames)):
        random_3utr[motifsNames[a]] = random_3utr[motifsNames[a]] / numCDS
        random_5utr[motifsNames[a]] = random_5utr[motifsNames[a]] / numCDS
        random_cds[motifsNames[a]] = random_cds[motifsNames[a]] / numCDS
        random_cDNA[motifsNames[a]] = random_cDNA[motifsNames[a]] / numCDS



    #Calculate statistics
    for a in range(len(motifsNames)):
        if initial_3utr[motifsNames[a]] > random_3utr[motifsNames[a]]:
            numHigher_3UTR[motifsNames[a]] += 1
        if initial_3utr[motifsNames[a]] < random_3utr[motifsNames[a]] :
            numLower_3UTR[motifsNames[a]] += 1

        if initial_5utr[motifsNames[a]] > random_5utr[motifsNames[a]]:
            numHigher_5UTR[motifsNames[a]] += 1
        if initial_5utr[motifsNames[a]] < random_5utr[motifsNames[a]]:
            numLower_5UTR[motifsNames[a]] += 1

        if initial_cds[motifsNames[a]] > random_cds[motifsNames[a]]:
            numHigher_CDS[motifsNames[a]] += 1
        if initial_cds[motifsNames[a]] < random_cds[motifsNames[a]]:
            numLower_CDS[motifsNames[a]] += 1

        if initial_cDNA[motifsNames[a]] > random_cDNA[motifsNames[a]]:
            numHigher_cDNA[motifsNames[a]] += 1
        if initial_cDNA[motifsNames[a]] < random_cDNA[motifsNames[a]]:
            numLower_cDNA[motifsNames[a]] += 1

outfile = open("randomizationResults.txt","w")

outfile.write("Statics for the 3UTR region\n")
outfile.write("Motif\tNum_higherInQuery\tNum_lowerInQuery\n")
for a in range(len(motifsNames)):
    outfile.write(motifsNames[a]+"\t"+str(numHigher_3UTR[motifsNames[a]])+"\t"+str(numLower_3UTR[motifsNames[a]])+"\n")

outfile.write("\nStatics for the 5UTR region\n")
outfile.write("Motif\tNum_higherInQuery\tNum_lowerInQuery\n")
for a in range(len(motifsNames)):
    outfile.write(motifsNames[a]+"\t"+str(numHigher_5UTR[motifsNames[a]])+"\t"+str(numLower_5UTR[motifsNames[a]])+"\n")

outfile.write("\nStatics for the CDS region\n")
outfile.write("Motif\tNum_higherInQuery\tNum_lowerInQuery\n")
for a in range(len(motifsNames)):
    outfile.write(motifsNames[a]+"\t"+str(numHigher_CDS[motifsNames[a]])+"\t"+str(numLower_CDS[motifsNames[a]])+"\n")

            
outfile.write("\nStatics for the cDNA region\n")
outfile.write("Motif\tNum_higherInQuery\tNum_lowerInQuery\n")
for a in range(len(motifsNames)):
    outfile.write(motifsNames[a]+"\t"+str(numHigher_cDNA[motifsNames[a]])+"\t"+str(numLower_cDNA[motifsNames[a]])+"\n")


outfile.close()



    

