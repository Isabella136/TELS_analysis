from math import log10

class indiv_stacked_abundance:
    def __init__(this):
        this.absolute_abundance = {}
        this.relative_abundance = {}
        this.classDict = {}

    def addToAbsolute(this, filepath, legend):
        argFile = open(filepath, "r")
        for line in argFile:
            blank_space = 0
            arg_header_list = line.split(',')[1:]
            for arg in arg_header_list:
                if arg == "":
                    blank_space += 1
                    if (blank_space % 2) == 0:
                        break
                    continue
                arg_header = arg.split("|")
                if arg_header[0] not in this.absolute_abundance:
                    this.absolute_abundance[arg_header[0]] = 0
                this.absolute_abundance[arg_header[0]] += 1
                className = arg_header[2]
                if arg_header[1] != "Drug":
                    className = arg_header[1] + " resistance"
                elif className == "betalactams":
                    className = "Betalactams"
                if className not in this.classDict:
                    this.classDict[className] = [arg_header[0]]
                elif arg_header[0] not in this.classDict[className]:
                    this.classDict[className].append(arg_header[0])
    
    def makeAbundanceRelative(this, file_size, gene_length):
        for arg, count in this.absolute_abundance.items():
            this.relative_abundance[arg] = log10((100.0 * float(count))/(float(gene_length[arg]*file_size)))

    def getAbundance(this):
        return this.relative_abundance