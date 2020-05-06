from Bio import SeqIO
import sys
import math
import random as rd
import numpy as nm
from scipy import stats




downGenes = set()
upGenes = set()

downGenesTranscripts = {}
upGenesTranscripts = {}

are2 = ["TTATTTAAA" ,"TTATTTATT", "TTATTTAAT","TTATTTATA" ]
puf = ["TGTAAATA","TGTATATA","TGTAGATA","TGTACATA"]

cdsTotest = sys.argv[1]
utr3ToTest = sys.argv[2]
utr5ToTest = sys.argv[3]
randomArg = int(sys.argv[4])

mobilityCDS = SeqIO.to_dict(SeqIO.parse(cdsTotest,"fasta"))
mobility3UTR = SeqIO.to_dict(SeqIO.parse(utr3ToTest,"fasta"))
mobility5UTR = SeqIO.to_dict(SeqIO.parse(utr5ToTest,"fasta"))

allCDS = SeqIO.to_dict(SeqIO.parse("cds_lv.fasta","fasta"))
all3UTR = SeqIO.to_dict(SeqIO.parse("3utr_lv.fasta","fasta"))
all5UTR = SeqIO.to_dict(SeqIO.parse("5utr_lv.fasta","fasta"))

numExons = {}
geneLen = {}

numExFile = open("Homo_sapiens.GRCh38.genes.all_length.txt")
while True:
    line = numExFile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    if not fields[0] in geneLen:
        geneLen[fields[0]] = fields[1]
numExFile.close()



#Collect number of exons
numExFile = open("numExons.txt")
while True:
    line = numExFile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    if not fields[0] in numExons:
        numExons[fields[0]] = int(fields[1])
numExFile.close()

def get_CpG(sequence):
    numGpC = 0
    for a in range(len(sequence)):
        if sequence[a:a+2] == "CG":
            numGpC += 1
    return numGpC

def get_ARE4(sequence):
    numARE4 = 0
    for a in range(len(sequence)-12):
        if sequence[a:a+12].count("G") + sequence[a:a+12].count("C") <=1:
            numARE4 += 1
            a +=12
    return numARE4

def get_ARE2(sequence):
    numARE2 = 0
    for a in range(len(sequence)-9):
        if sequence[a:a+9] in are2:
            numARE2 += 1
            a += 9
    return numARE2

def get_puf(sequence):
    numPUF = 0
    for a in range(len(sequence)-8):
        if sequence[a:a+8] in puf:
            numPUF += 1
            a +=8
    return numPUF

def get_numExons(seqID):
    return numExons[(seqID.split("_"))[1]]

def get_geneLength(seqID):
    return geneLen[(seqID.split("_"))[0]]

def get_GCcont(sequence):
    gc = float(sequence.count("G") + sequence.count("C") + sequence.count("g") + sequence.count("c"))/float(len(sequence))
    return gc

def get_CPEB(sequence):
    numCpeb = 0
    for a in range(len(sequence)-7):
        if sequence[a:a+7] == "TTTTTAT":
            numCpeb += 1
            a +=7
    return numCpeb


def get_StabilityValue(numExons,geneLength,mRNALength,ARE4_3utr,ARE4_mRNA,ARE2_3utr,CpG_5utr,Puf_3utr,Puf_cds):
    stability = -0.592 -0.179*(math.log( float(numExons)/float(geneLength) *1000)) -0.02*math.log(mRNALength)   
    stability += 0.054*math.log(ARE4_3utr+1)
    stability += 0.022*math.log(ARE4_mRNA+1)
    stability += 0.040*math.log(ARE2_3utr+1)
    stability += 0.051*math.log(CpG_5utr+1)
    stability += 0.071*math.log(Puf_3utr+1) 
    stability += 0.044*math.log(Puf_cds+1)
    
    
    
    return stability

CpG_5utr_higher = 0
ARE4_5utr_higher = 0
ARE4_3utr_higher = 0
ARE2_3utr_higher = 0
Puf_3utr_higher = 0
ARE4_cds_higher = 0
Puf_cds_higher = 0
NumExons_higher = 0
GC_mRNA_higher = 0
GC_cds_higher = 0
GC_3utr_higher = 0
GC_5utr_higher = 0
Stability_higher = 0
Cpeb_3utr_higher = 0
Cpeb_cds_higher = 0
Cpeb_5utr_higher = 0


CpG_5utr_lower = 0
ARE4_5utr_lower = 0
ARE4_3utr_lower = 0
ARE2_3utr_lower = 0
Puf_3utr_lower = 0
ARE4_cds_lower = 0
Puf_cds_lower = 0
NumExons_lower = 0
GC_mRNA_lower = 0
GC_cds_lower = 0
GC_3utr_lower = 0
GC_5utr_lower = 0
Stability_lower = 0
Cpeb_3utr_lower = 0
Cpeb_cds_lower = 0
Cpeb_5utr_lower = 0


CpG_5utr_higher_sign = 0
ARE4_5utr_higher_sign = 0
ARE4_3utr_higher_sign = 0
ARE2_3utr_higher_sign = 0
Puf_3utr_higher_sign = 0
ARE4_cds_higher_sign = 0
Puf_cds_higher_sign = 0
NumExons_higher_sign = 0
GC_mRNA_higher_sign = 0
GC_cds_higher_sign = 0
GC_3utr_higher_sign = 0
GC_5utr_higher_sign = 0
Stability_higher_sign = 0
Cpeb_3utr_higher_sign = 0
Cpeb_cds_higher_sign = 0
Cpeb_5utr_higher_sign = 0



CpG_5utr_lower_sign = 0
ARE4_5utr_lower_sign = 0
ARE4_3utr_lower_sign = 0
ARE2_3utr_lower_sign = 0
Puf_3utr_lower_sign = 0
ARE4_cds_lower_sign = 0
Puf_cds_lower_sign = 0
NumExons_lower_sign = 0
GC_mRNA_lower_sign = 0
GC_cds_lower_sign = 0
GC_3utr_lower_sign = 0
GC_5utr_lower_sign = 0
Stability_lower_sign = 0
Cpeb_3utr_lower_sign = 0
Cpeb_cds_lower_sign = 0
Cpeb_5utr_lower_sign = 0








outfile = open("randomizationResults.txt","w")
outfile2 = open("randomizationAverageResults.txt","w")
outfile3 = open("anovaResults.txt","w")
outfile.write("Randomization\tTranscript\tCpG_5utr\tARE4_5utr\tARE4_3utr\tARE2_3utr\tPuf_3utr\tARE4_cds\tPuf_cds\tNumExons\t3UTR_GC\t5UTR_GC\tCDS_GC\tmRNA_GC\tStability\tCpeb_3utr\tCpeb_5utr\tCpeb_cds\n")
outfile2.write("Randomization\tCpG_5utr\tARE4_5utr\tARE4_3utr\tARE2_3utr\tPuf_3utr\tARE4_cds\tPuf_cds\tNumExons\t3UTR_GC\t5UTR_GC\tCDS_GC\tmRNA_GC\tStability\tCpeb_3utr\tCpeb_5utr\tCpeb_cds\n")
outfile3.write("Parameter\tHigher\tHigherSign\tLower\tLowerSign\n")


#Calculate initial values for the 108 (now 104 genes)
CpG_5utr_mobility_list = []
ARE4_5utr_mobility_list = []
ARE4_3utr_mobility_list = []
ARE2_3utr_mobility_list = []
Puf_3utr_mobility_list = []
ARE4_cds_mobility_list = []
Puf_cds_mobility_list = []
NumExons_mobility_list = []
GC_mRNA_mobility_list = []
GC_cds_mobility_list = []
GC_3utr_mobility_list = []
GC_5utr_mobility_list = []
stability_mobility_list = []
Cpeb_3utr_mobility_list = [] 
Cpeb_5utr_mobility_list = [] 
Cpeb_cds_mobility_list = [] 


for transcript in mobilityCDS:

    #5UTR statistics
    if transcript in all5UTR:
        CpG_5utr_mobility = get_CpG(all5UTR[transcript].seq)
        CpG_5utr_mobility_list.append(CpG_5utr_mobility)

        ARE4_5utr_mobility = get_ARE4(all5UTR[transcript].seq)
        ARE4_5utr_mobility_list.append(ARE4_5utr_mobility)

        GC_5utr_mobility = get_GCcont(all5UTR[transcript].seq)
        GC_5utr_mobility_list.append(GC_5utr_mobility)

        Cpeb_5utr_mobility = get_CPEB(all5UTR[transcript].seq)
        Cpeb_5utr_mobility_list.append(Cpeb_5utr_mobility)

    #3UTR statistics
    if transcript in all3UTR:
        ARE4_3utr_mobility = get_ARE4(all3UTR[transcript].seq)
        ARE4_3utr_mobility_list.append(ARE4_3utr_mobility)

        ARE2_3utr_mobility = get_ARE2(all3UTR[transcript].seq)
        ARE2_3utr_mobility_list.append(ARE2_3utr_mobility)

        Puf_3utr_mobility = get_puf(all3UTR[transcript].seq)
        Puf_3utr_mobility_list.append(Puf_3utr_mobility)

        GC_3utr_mobility = get_GCcont(all3UTR[transcript].seq)
        GC_3utr_mobility_list.append(GC_3utr_mobility)

        Cpeb_3utr_mobility = get_CPEB(all3UTR[transcript].seq)
        Cpeb_3utr_mobility_list.append(Cpeb_3utr_mobility)

    #CDS statistics
    if transcript in allCDS:
        ARE4_cds_mobility = get_ARE4(allCDS[transcript].seq)
        ARE4_cds_mobility_list.append(ARE4_cds_mobility)

        Puf_cds_mobility = get_puf(allCDS[transcript].seq)
        Puf_cds_mobility_list.append(Puf_cds_mobility)

        mRNASeq = allCDS[transcript].seq
        if transcript in all3UTR:
            mRNASeq += all3UTR[transcript].seq
        if transcript in all5UTR:
            mRNASeq += all5UTR[transcript].seq

        GC_mRNA_mobility = get_GCcont(mRNASeq)
        GC_mRNA_mobility_list.append(GC_mRNA_mobility)

        GC_cds_mobility = get_GCcont(allCDS[transcript].seq)
        GC_cds_mobility_list.append(GC_cds_mobility)

        Cpeb_cds_mobility = get_CPEB(allCDS[transcript].seq)
        Cpeb_cds_mobility_list.append(Cpeb_cds_mobility)

    NumExons = get_numExons(transcript)
    NumExons_mobility_list.append(NumExons)

    geneLength = get_geneLength(transcript)

    stability_mobility = get_StabilityValue(NumExons,geneLength,len(allCDS[transcript].seq),ARE4_3utr_mobility,ARE4_cds_mobility,ARE2_3utr_mobility,CpG_5utr_mobility,Puf_3utr_mobility,Puf_cds_mobility)
    stability_mobility_list.append(stability_mobility)

    outfile.write("MobilityDataset"+"\t"+transcript+"\t" + str(CpG_5utr_mobility)+"\t"+str(ARE4_5utr_mobility)+"\t"+str(ARE4_3utr_mobility)+"\t"+str(ARE2_3utr_mobility)+"\t"+str(Puf_3utr_mobility)+"\t"+str(ARE4_cds_mobility)+"\t"+str(Puf_cds_mobility)+"\t"+str(NumExons)+"\t"+str(GC_3utr_mobility)+"\t"+str(GC_5utr_mobility)+"\t"+str(GC_cds_mobility)+"\t"+str(GC_mRNA_mobility)+"\t"+str(stability_mobility)+"\t"+str(Cpeb_3utr_mobility)+"\t"+str(Cpeb_5utr_mobility)+"\t"+str(Cpeb_cds_mobility)+"\n")


avMobilityCpG_5utr= nm.mean(CpG_5utr_mobility_list)
avMobilityARE4_5utr = nm.mean(ARE4_5utr_mobility_list)
avMobilityARE4_3utr = nm.mean(ARE4_3utr_mobility_list)
avMobilityARE2_3utr = nm.mean(ARE2_3utr_mobility_list)
avMobilityPuf_3utr = nm.mean(Puf_3utr_mobility_list)
avMobilityARE4_cds = nm.mean(ARE4_cds_mobility_list)
avMobilityPuf_cds = nm.mean(Puf_cds_mobility_list)
avMobilityNumExons = nm.mean(NumExons_mobility_list)
avMobility_mRNA_GC = nm.mean(GC_mRNA_mobility_list)
avMobility_3utr_GC = nm.mean(GC_3utr_mobility_list)
avMobility_5utr_GC = nm.mean(GC_5utr_mobility_list)
avMobility_cds_GC = nm.mean(GC_cds_mobility_list)
avMobility_stability = nm.mean(stability_mobility_list)
avMobility_Cpeb_3utr = nm.mean(Cpeb_3utr_mobility_list)
avMobility_Cpeb_5utr = nm.mean(Cpeb_5utr_mobility_list)
avMobility_Cpeb_cds = nm.mean(Cpeb_cds_mobility_list)



outfile2.write("MobilityDataset"+"\t"+str(avMobilityCpG_5utr)+"\t"+str(avMobilityARE4_5utr)+"\t"+str(avMobilityARE4_3utr)+"\t"+str(avMobilityARE2_3utr)+"\t"+str(avMobilityPuf_3utr)+"\t"+str(avMobilityARE4_cds)+"\t"+str(avMobilityPuf_cds)+"\t"+str(avMobilityNumExons)+"\t"+str(avMobility_3utr_GC)+"\t"+str(avMobility_5utr_GC)+"\t"+str(avMobility_cds_GC)+"\t"+str(avMobility_mRNA_GC)+"\t"+str(avMobility_stability)+"\n")


for numRandomization in range(randomArg):
    foundGenes_list = []
    CpG_5utr_list = []
    ARE4_5utr_list = []
    ARE4_3utr_list = []
    ARE2_3utr_list = []
    Puf_3utr_list = []
    ARE4_cds_list = []
    Puf_cds_list = []
    CpG_5utr_list = []
    CpG_3utr_list = []
    CpG_cds_list = []
    CpG_transcript_list = []
    NumExons_list = []
    GC_3utr_list = []
    GC_5utr_list = []
    GC_cds_list = []
    GC_mRNA_list = []
    stability_list = [] 
    Cpeb_3utr_list = [] 
    Cpeb_5utr_list = [] 
    Cpeb_cds_list = [] 





    print "Randomization",numRandomization
    foundGenes = 0
    while foundGenes <104:
        randomTranscript = rd.choice(allCDS.keys())
        if randomTranscript in allCDS and randomTranscript in all3UTR and randomTranscript in all5UTR:
            foundGenes+=1
    
            CpG_5utr = get_CpG(all5UTR[randomTranscript].seq)
            CpG_5utr_list.append(CpG_5utr)
            ARE4_5utr = get_ARE4(all5UTR[randomTranscript].seq)
            ARE4_5utr_list.append(ARE4_5utr)
            ARE4_3utr = get_ARE4(all3UTR[randomTranscript].seq)
            ARE4_3utr_list.append(ARE4_3utr)
            ARE2_3utr = get_ARE2(all3UTR[randomTranscript].seq)
            ARE2_3utr_list.append(ARE2_3utr)
            Puf_3utr = get_puf(all3UTR[randomTranscript].seq)
            Puf_3utr_list.append(Puf_3utr)
            ARE4_cds = get_ARE4(allCDS[randomTranscript].seq)
            ARE4_cds_list.append(ARE4_cds)
            Puf_cds = get_puf(allCDS[randomTranscript].seq)
            Puf_cds_list.append(Puf_cds)
            NumExons = get_numExons(randomTranscript)
            NumExons_list.append(NumExons)
            GC_3utr = get_GCcont(all3UTR[randomTranscript].seq)
            GC_3utr_list.append(GC_3utr)
            GC_5utr = get_GCcont(all5UTR[randomTranscript].seq)
            GC_5utr_list.append(GC_5utr)
            GC_cds = get_GCcont(allCDS[randomTranscript].seq)
            GC_cds_list.append(GC_cds)
            mRNAseq = allCDS[randomTranscript].seq+all3UTR[randomTranscript].seq+all5UTR[randomTranscript].seq
            GC_mRNA = get_GCcont(mRNAseq)
            GC_mRNA_list.append(GC_mRNA)

            geneLength = get_geneLength(randomTranscript)

            stability = get_StabilityValue(NumExons,geneLength,len(allCDS[randomTranscript].seq),ARE4_3utr_mobility,ARE4_cds_mobility,ARE2_3utr_mobility,CpG_5utr_mobility,Puf_3utr_mobility,Puf_cds_mobility)
            stability_list.append(stability)

            Cpeb_3utr = get_CPEB(all3UTR[randomTranscript].seq)
            Cpeb_3utr_list.append(Cpeb_3utr)

            Cpeb_5utr = get_CPEB(all5UTR[randomTranscript].seq)
            Cpeb_5utr_list.append(Cpeb_5utr)

            Cpeb_cds = get_CPEB(allCDS[randomTranscript].seq)
            Cpeb_cds_list.append(Cpeb_cds)

            outfile.write("randomization"+str(numRandomization)+"\t"+randomTranscript+"\t" + str(CpG_5utr)+"\t"+str(ARE4_5utr)+"\t"+str(ARE4_3utr)+"\t"+str(ARE2_3utr)+"\t"+str(Puf_3utr)+"\t"+str(ARE4_cds)+"\t"+str(Puf_cds)+"\t"+str(NumExons)+"\t"+str(GC_3utr)+"\t"+str(GC_5utr)+"\t"+str(GC_cds)+"\t"+str(GC_mRNA)+"\t"+str(stability)+"\t"+str(Cpeb_3utr)+"\t"+str(Cpeb_5utr)+"\t"+str(Cpeb_cds)+"\n")


    avCpG_5utr= nm.mean(CpG_5utr_list)
    avARE4_5utr = nm.mean(ARE4_5utr_list)
    avARE4_3utr = nm.mean(ARE4_3utr_list)
    avARE2_3utr = nm.mean(ARE2_3utr_list)
    avPuf_3utr = nm.mean(Puf_3utr_list)
    avARE4_cds = nm.mean(ARE4_cds_list)
    avPuf_cds = nm.mean(Puf_cds_list)
    avNumExons = nm.mean(NumExons_list)
    av_mRNA_GC = nm.mean(GC_mRNA_list)
    av_3utr_GC = nm.mean(GC_3utr_list)
    av_5utr_GC = nm.mean(GC_5utr_list)
    av_cds_GC = nm.mean(GC_cds_list)
    avStability = nm.mean(stability_list)
    av_Cpeb_3utr = nm.mean(Cpeb_3utr_list)
    av_Cpeb_5utr = nm.mean(Cpeb_5utr_list)
    av_Cpeb_cds = nm.mean(Cpeb_cds_list)

    outfile2.write("randomization"+str(numRandomization)+"\t"+str(avCpG_5utr)+"\t"+str(avARE4_5utr)+"\t"+str(avARE4_3utr)+"\t"+str(avARE2_3utr)+"\t"+str(avPuf_3utr)+"\t"+str(avARE4_cds)+"\t"+str(avPuf_cds)+"\t"+str(avNumExons)+"\t"+str(av_3utr_GC)+"\t"+str(av_5utr_GC)+"\t"+str(av_cds_GC)+"\t"+str(av_mRNA_GC)+"\t"+str(avStability)+"\t"+str(av_Cpeb_3utr)+"\t"+str(av_Cpeb_5utr)+"\t"+str(av_Cpeb_cds)+"\n")
    
    #Start Statistical analysis
    #GC content
    if avMobility_3utr_GC > av_3utr_GC:
        GC_3utr_higher += 1
        (ttest,pvalue) = stats.ttest_ind(GC_3utr_mobility_list ,GC_3utr_list)
        if pvalue <= 0.05:
            GC_3utr_higher_sign += 1

    if avMobility_5utr_GC > av_5utr_GC:
        GC_5utr_higher += 1
        (ttest,pvalue) = stats.ttest_ind(GC_5utr_mobility_list ,GC_5utr_list)
        if pvalue <= 0.05:
            GC_5utr_higher_sign += 1

    if avMobility_cds_GC > av_cds_GC:
        GC_cds_higher += 1
        (ttest,pvalue) = stats.ttest_ind(GC_cds_mobility_list ,GC_cds_list)
        if pvalue <= 0.05:
            GC_cds_higher_sign += 1


    if avMobility_3utr_GC < av_3utr_GC:
        GC_3utr_lower += 1
        (ttest,pvalue) = stats.ttest_ind(GC_3utr_mobility_list ,GC_3utr_list)
        if pvalue <= 0.05:
            GC_3utr_lower_sign += 1

    if avMobility_5utr_GC < av_5utr_GC:
        GC_5utr_lower += 1
        (ttest,pvalue) = stats.ttest_ind(GC_5utr_mobility_list ,GC_5utr_list)
        if pvalue <= 0.05:
            GC_5utr_lower_sign += 1

    if avMobility_cds_GC < av_cds_GC:
        GC_cds_lower += 1
        (ttest,pvalue) = stats.ttest_ind(GC_cds_mobility_list ,GC_cds_list)
        if pvalue <= 0.05:
            GC_cds_lower_sign += 1


    if avMobility_mRNA_GC > av_mRNA_GC:
        GC_mRNA_higher += 1
        (ttest,pvalue) = stats.ttest_ind(GC_mRNA_mobility_list ,GC_mRNA_list)
        if pvalue <= 0.05:
            GC_mRNA_higher_sign += 1
    else:
        GC_mRNA_lower += 1
        (ttest,pvalue) = stats.ttest_ind(GC_mRNA_mobility_list ,GC_mRNA_list)
        if pvalue <= 0.05:
            GC_mRNA_lower_sign += 1



    #CpG_5utr
    if avMobilityCpG_5utr > avCpG_5utr:
        CpG_5utr_higher += 1
        (ttest,pvalue) = stats.ttest_ind(CpG_5utr_mobility_list ,CpG_5utr_list )
        if pvalue <= 0.05:
            CpG_5utr_higher_sign += 1
    else:
        CpG_5utr_lower += 1
        (ttest,pvalue) = stats.ttest_ind(CpG_5utr_mobility_list ,CpG_5utr_list )
        if pvalue <= 0.05:
            CpG_5utr_lower_sign += 1

    #ARE4_5utr
    if avMobilityARE4_5utr > avARE4_5utr:
        ARE4_5utr_higher += 1
        (ttest,pvalue) = stats.ttest_ind(ARE4_5utr_mobility_list ,ARE4_5utr_list )
        if pvalue <= 0.05:
            ARE4_5utr_higher_sign += 1
    else:
        ARE4_5utr_lower += 1
        (ttest,pvalue) = stats.ttest_ind(ARE4_5utr_mobility_list ,ARE4_5utr_list )
        if pvalue <= 0.05:
            ARE4_5utr_lower_sign += 1

    #ARE4_3utr
    if avMobilityARE4_3utr > avARE4_3utr:
        ARE4_3utr_higher += 1
        (ttest,pvalue) = stats.ttest_ind(ARE4_3utr_mobility_list ,ARE4_3utr_list )
        if pvalue <= 0.05:
            ARE4_3utr_higher_sign += 1
    else:
        ARE4_3utr_lower += 1
        (ttest,pvalue) = stats.ttest_ind(ARE4_3utr_mobility_list ,ARE4_3utr_list )
        if pvalue <= 0.05:
            ARE4_3utr_lower_sign += 1

    #ARE4_cds
    if avMobilityARE4_cds > avARE4_cds:
        ARE4_cds_higher += 1
        (ttest,pvalue) = stats.ttest_ind(ARE4_cds_mobility_list ,ARE4_cds_list )
        if pvalue <= 0.05:
            ARE4_cds_higher_sign += 1
    else:
        ARE4_cds_lower += 1
        (ttest,pvalue) = stats.ttest_ind(ARE4_cds_mobility_list ,ARE4_cds_list )
        if pvalue <= 0.05:
            ARE4_cds_lower_sign += 1

    #ARE2_3utr
    if avMobilityARE2_3utr > avARE2_3utr:
        ARE2_3utr_higher += 1
        (ttest,pvalue) = stats.ttest_ind(ARE2_3utr_mobility_list ,ARE2_3utr_list )
        if pvalue <= 0.05:
            ARE2_3utr_higher_sign += 1
    else:
        ARE2_3utr_lower += 1
        (ttest,pvalue) = stats.ttest_ind(ARE2_3utr_mobility_list ,ARE2_3utr_list )
        if pvalue <= 0.05:
            ARE2_3utr_lower_sign += 1

    #Puf_3utr
    if avMobilityPuf_3utr > avPuf_3utr:
        Puf_3utr_higher += 1
        (ttest,pvalue) = stats.ttest_ind(Puf_3utr_mobility_list ,Puf_3utr_list )
        if pvalue <= 0.05:
            Puf_3utr_higher_sign += 1
    else:
        Puf_3utr_lower += 1
        (ttest,pvalue) = stats.ttest_ind(Puf_3utr_mobility_list ,Puf_3utr_list )
        if pvalue <= 0.05:
            Puf_3utr_lower_sign += 1

    #Puf_cds
    if avMobilityPuf_cds > avPuf_cds:
        Puf_cds_higher += 1
        (ttest,pvalue) = stats.ttest_ind(Puf_cds_mobility_list ,Puf_cds_list )
        if pvalue <= 0.05:
            Puf_cds_higher_sign += 1
    else:
        Puf_cds_lower += 1
        (ttest,pvalue) = stats.ttest_ind(Puf_cds_mobility_list ,Puf_cds_list )
        if pvalue <= 0.05:
            Puf_cds_lower_sign += 1

    #numExons
    if avMobilityNumExons > avNumExons:
        NumExons_higher += 1
        (ttest,pvalue) = stats.ttest_ind(NumExons_mobility_list ,NumExons_list )
        if pvalue <= 0.05:
            NumExons_higher_sign += 1
    else:
        NumExons_lower += 1
        (ttest,pvalue) = stats.ttest_ind(NumExons_mobility_list ,NumExons_list )
        if pvalue <= 0.05:
            NumExons_lower_sign += 1


    #stability
    if avMobility_stability > avStability:
        Stability_higher += 1
        (ttest,pvalue) = stats.ttest_ind(stability_mobility_list ,stability_list )
        if pvalue <= 0.05:
            Stability_higher_sign += 1
    else:
        Stability_lower += 1
        (ttest,pvalue) = stats.ttest_ind(stability_mobility_list ,stability_list )
        if pvalue <= 0.05:
            Stability_lower_sign += 1

    #Cpeb_3utr
    if avMobility_Cpeb_3utr > av_Cpeb_3utr:
        Cpeb_3utr_higher += 1
        (ttest,pvalue) = stats.ttest_ind(Cpeb_3utr_mobility_list ,Cpeb_3utr_list )
        if pvalue <= 0.05:
            Cpeb_3utr_higher_sign += 1
    else:
        Cpeb_3utr_lower += 1
        (ttest,pvalue) = stats.ttest_ind(Cpeb_3utr_mobility_list ,Cpeb_3utr_list )
        if pvalue <= 0.05:
            Cpeb_3utr_lower_sign += 1

    #Cpeb_5utr
    if avMobility_Cpeb_5utr > av_Cpeb_5utr:
        Cpeb_5utr_higher += 1
        (ttest,pvalue) = stats.ttest_ind(Cpeb_5utr_mobility_list ,Cpeb_5utr_list )
        if pvalue <= 0.05:
            Cpeb_5utr_higher_sign += 1
    else:
        Cpeb_5utr_lower += 1
        (ttest,pvalue) = stats.ttest_ind(Cpeb_5utr_mobility_list ,Cpeb_5utr_list )
        if pvalue <= 0.05:
            Cpeb_5utr_lower_sign += 1

    
    #Cpeb_cds
    if avMobility_Cpeb_cds > av_Cpeb_cds:
        Cpeb_cds_higher += 1
        (ttest,pvalue) = stats.ttest_ind(Cpeb_cds_mobility_list ,Cpeb_cds_list )
        if pvalue <= 0.05:
            Cpeb_cds_higher_sign += 1
    else:
        Cpeb_cds_lower += 1
        (ttest,pvalue) = stats.ttest_ind(Cpeb_cds_mobility_list ,Cpeb_cds_list )
        if pvalue <= 0.05:
            Cpeb_cds_lower_sign += 1

   




outfile3.write("GC_3utr\t"+str(GC_3utr_higher)+"\t"+str(GC_3utr_higher_sign)+"\t"+str(GC_3utr_lower)+"\t"+str(GC_3utr_lower_sign)+"\n")
outfile3.write("GC_5utr\t"+str(GC_5utr_higher)+"\t"+str(GC_5utr_higher_sign)+"\t"+str(GC_5utr_lower)+"\t"+str(GC_5utr_lower_sign)+"\n")
outfile3.write("CpG_5utr\t"+str(CpG_5utr_higher)+"\t"+str(CpG_5utr_higher_sign)+"\t"+str(CpG_5utr_lower)+"\t"+str(CpG_5utr_lower_sign)+"\n")
outfile3.write("GC_cds\t"+str(GC_cds_higher)+"\t"+str(GC_cds_higher_sign)+"\t"+str(GC_cds_lower)+"\t"+str(GC_cds_lower_sign)+"\n")
outfile3.write("GC_mRNA\t"+str(GC_mRNA_higher)+"\t"+str(GC_mRNA_higher_sign)+"\t"+str(GC_mRNA_lower)+"\t"+str(GC_mRNA_lower_sign)+"\n")
outfile3.write("ARE4_3utr\t"+str(ARE4_3utr_higher)+"\t"+str(ARE4_3utr_higher_sign)+"\t"+str(ARE4_3utr_lower)+"\t"+str(ARE4_3utr_lower_sign)+"\n")
outfile3.write("ARE4_5utr\t"+str(ARE4_5utr_higher)+"\t"+str(ARE4_5utr_higher_sign)+"\t"+str(ARE4_5utr_lower)+"\t"+str(ARE4_5utr_lower_sign)+"\n")
outfile3.write("ARE4_cds\t"+str(ARE4_cds_higher)+"\t"+str(ARE4_cds_higher_sign)+"\t"+str(ARE4_cds_lower)+"\t"+str(ARE4_cds_lower_sign)+"\n")
outfile3.write("ARE2_3utr\t"+str(ARE2_3utr_higher)+"\t"+str(ARE2_3utr_higher_sign)+"\t"+str(ARE2_3utr_lower)+"\t"+str(ARE2_3utr_lower_sign)+"\n")
outfile3.write("PUF_3utr\t"+str(Puf_3utr_higher)+"\t"+str(Puf_3utr_higher_sign)+"\t"+str(Puf_3utr_lower)+"\t"+str(Puf_3utr_lower_sign)+"\n")
outfile3.write("PUF_cds\t"+str(Puf_cds_higher)+"\t"+str(Puf_cds_higher_sign)+"\t"+str(Puf_cds_lower)+"\t"+str(Puf_cds_lower_sign)+"\n")
outfile3.write("Stability\t"+str(Stability_higher)+"\t"+str(Stability_higher_sign)+"\t"+str(Stability_lower)+"\t"+str(Stability_lower_sign)+"\n")
outfile3.write("Cpeb_3utr\t"+str(Cpeb_3utr_higher)+"\t"+str(Cpeb_3utr_higher_sign)+"\t"+str(Cpeb_3utr_lower)+"\t"+str(Cpeb_3utr_lower_sign)+"\n")
outfile3.write("Cpeb_5utr\t"+str(Cpeb_5utr_higher)+"\t"+str(Cpeb_5utr_higher_sign)+"\t"+str(Cpeb_5utr_lower)+"\t"+str(Cpeb_5utr_lower_sign)+"\n")
outfile3.write("Cpeb_cds\t"+str(Cpeb_cds_higher)+"\t"+str(Cpeb_cds_higher_sign)+"\t"+str(Cpeb_cds_lower)+"\t"+str(Cpeb_cds_lower_sign)+"\n")
outfile3.write("NumExons\t"+str(NumExons_higher)+"\t"+str(NumExons_higher_sign)+"\t"+str(NumExons_lower)+"\t"+str(NumExons_lower_sign)+"\n")





outfile.close()
outfile2.close()
outfile3.close()




        

        








