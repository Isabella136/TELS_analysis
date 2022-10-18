from math import log10

class indiv_abundance:
    def __init__(this, sample_type, x_axis):
        this.x_axis = x_axis
        this.sample_type = sample_type
        this.absolute_abundance = {}
        this.relative_abundance = {}

    def addToAbsolute(this, filepath):
        argFile = open(filepath, "r")
        lineNum = 0
        for line in argFile:
            lineNum += 1
            if lineNum < 20:
                continue
            line_list = line.split(',')
            arg_header = line_list[0].split("|")
            if arg_header[0] not in this.absolute_abundance:
                this.absolute_abundance[arg_header[0]] = 0
            this.absolute_abundance[arg_header[0]] += int(line_list[1])
    def makeAbundanceRelative(this, file_size, gene_length):
        for arg, count in this.absolute_abundance.items():
            this.relative_abundance[arg] = log10((100.0 * float(count))/(float(gene_length[arg]*file_size)))
    def getAbundance(this):
        return this.relative_abundance