from Bio import SeqIO
import csv

# lambda that generates tels output file
tels_file_path = (lambda self, sample_name, extension : 
             self.source_prefix + sample_name + self.source_suffix + extension)


def fileDict(sample_name):
    sample = ""
    seqPlatform = ""

    #Determine organism
    if sample_name[0] == 'B':
        sample = "Bovine fecal"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'H':
        sample = "Human fecal"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'M':
        sample = "Mock"
        sample_name = sample_name[2:]
    else: #sample_name[0] = 'S'
        sample = "Soil"
        sample_name = sample_name[1:]

    #Determine whether V2 or XT
    if sample_name[0] == 'V':
        sample = sample + " (V2)"
        sample_name = sample_name[2:]
    else: #sample_name[0] = 'X'
        sample = sample + " (XT)"
        sample_name = sample_name[2:]

    #Determine probes
    if sample_name[0:2] == "AM":
        sample = sample + " + ARG-MGE probe"
        seqPlatform = "TELSeq"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'A':
        sample = sample + " + ARG probe"
        seqPlatform = "TELSeq"
        sample_name = sample_name[1:]
    elif sample_name[0] == 'M':
        sample = sample + " + MGE probe"
        seqPlatform = "TELSeq"
        sample_name = sample_name[1:]
    else: #sample_name[0] == 'N'
        seqPlatform = "PacBio"
        return (sample, seqPlatform)

    #Determine baits
    sample = sample + sample_name
    return (sample, seqPlatform)

def megares_analyzer(megaresFile):
    def sortList(mechanismsList):
        if len(mechanismsList) == 1:
            return mechanismsList
        else:
            listA = sortList(mechanismsList[0:int(len(mechanismsList)/2)])
            listB = sortList(mechanismsList[int(len(mechanismsList)/2):])
            toReturn = []
            for i in range(0,len(mechanismsList)):
                if len(listA) == 0:
                    toReturn.append(listB.pop(0))
                elif len(listB) == 0:
                    toReturn.append(listA.pop(0))
                elif listA[0][0] < listB[0][0]:
                    toReturn.append(listA.pop(0))
                else:
                    toReturn.append(listB.pop(0))
            return toReturn
    drugList = []
    otherList = []
    megares = open(megaresFile, "r")
    megares.readline()
    for line in megares:
        splitLine = line.split('|')
        tempTuple = (splitLine[2], splitLine[3])
        if splitLine[2] == "betalactams": tempTuple = ("Betalactams", splitLine[3])
        if splitLine[1] == "Drugs":
            if drugList.count(tempTuple) == 0:
                drugList.append(tempTuple)
        else:
            if otherList.count(tempTuple) == 0:
                otherList.append(tempTuple)
    megares.close()
    drugList = sortList(drugList)
    otherList = sortList(otherList)
    return (drugList, otherList)

def get_genes_length(fasta_file):
    gene_length_dict = dict()
    with open(fasta_file, "r") as fasta:
        fasta_reader = SeqIO.parse(fasta, 'fasta')
        for reference in fasta_reader:
            gene_length_dict.update({reference.id: len(reference.seq)})
    return gene_length_dict

# Returns sample name definition in tuple form:
# (Organism, Platform, Chemistry, Probe)
def get_sample_name_definition(sample_name):
    # Determine organism
    if sample_name[0] == 'B':
        organism = "Bovine"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'H':
        organism = "Human"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'M':
        organism = "Mock"
        sample_name = sample_name[2:]
    else: #sample_name[0] = 'S'
        organism = "Soil"
        sample_name = sample_name[1:]

    # Determine chemistry
    chemistry = sample_name[:2]
    sample_name = sample_name[2:]
    
    # Determine probe specificity and
    # sequencing platform
    if sample_name[0:2] == "AM":
        platform = "TELSeq"
        probe = "Combo"
    elif sample_name[0] == 'A':
        platform = "TELSeq"
        probe = "ARG"
    elif sample_name[0] == 'M':
        platform = "TELSeq"
        probe = "MGE"
    else: #sample_name[0] == 'N'
        platform = "PacBio"
        probe = None

    return (organism, platform, chemistry, probe)

def getSampleAndIndex(sample_name):
    sub_table = ""
    index = 0

    #Determine organism
    if sample_name[0] == 'B':
        sub_table = "Bovine"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'H':
        sub_table = "Human"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'M':
        sub_table = "Mock"
        sample_name = sample_name[2:]
    else: #sample_name[0] = 'S'
        sub_table = "Soil"
        sample_name = sample_name[1:]

    #Determine whether V2 or XT
    if sample_name[0] == 'V':
        sub_table = sub_table + "+V2"
        sample_name = sample_name[2:]
    else: #sample_name[0] = 'X'
        sub_table = sub_table + "+XT"
        sample_name = sample_name[2:]

    #Determine probes
    if sample_name[0:2] == "AM":
        index = 2
    elif sample_name[0] == 'A':
        index = 1
    elif sample_name[0] == 'M':
        index = 3
    else: #sample_name[0] == 'N'
        index = 4
    return (sub_table,index)

def get_mge_annot_dict(filepath):
    mge_annot = dict()
    with open(filepath, "r") as csv_file:
        csv_reader = csv.reader(csv_file)
        for line_num, line in enumerate(csv_reader):
            if line_num == 0: continue
            mge_annot[line[0]] = line[1]
    return mge_annot