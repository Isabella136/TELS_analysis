from math import log10

class indiv_stacked_abundance:
    def __init__(this):
        this.absolute_abundance = {}
        this.relative_abundance = {}
        this.classDict = {}
        this.classList = []
        this.statsFilepaths = {}

    def addToAbsolute(this, filepath, legend, stats):
        if legend not in this.statsFilepaths:
            this.statsFilepaths[legend] = []
        this.statsFilepaths[legend].append(stats)
        argFile = open(filepath, "r")
        if legend not in this.absolute_abundance:
            this.absolute_abundance[legend] = {}
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
                if arg_header[4] not in list(this.absolute_abundance[legend].keys()):
                    this.absolute_abundance[legend][arg_header[4]] = 0
                this.absolute_abundance[legend][arg_header[4]] += 1
                className = arg_header[2]
                if arg_header[1] != "Drugs":
                    className = arg_header[1] + " resistance"
                elif className == "betalactams":
                    className = "Betalactams"
                if arg_header[4] not in this.classDict:
                    this.classDict[arg_header[4]] = className
                if className not in this.classList:
                    this.classList.append(className)
    
    def makeAbundanceRelative(this):
        avgReads = {}
        for legend, list in this.statsFilepaths.items():
            for filepath in list:
                statFile = open(filepath, "r")
                statFile.readline()
                readCount = int(statFile.readline().split(',')[1][:-1])
                if legend not in avgReads:
                    avgReads[legend] = 0
                avgReads[legend]+=readCount
                statFile.close()
            avgReads[legend] = avgReads[legend]/3
        for legend, dict in this.absolute_abundance.items():
            for group, count in dict.items():
                if legend not in this.relative_abundance:
                    this.relative_abundance[legend] = {}
                this.relative_abundance[legend].update({group:log10(count/avgReads[legend]*1000000)})

    def getAbundance(this):
        toReturn = {}
        for arg in this.classDict:
            for legend, dict in this.relative_abundance.items():
                if legend not in toReturn:
                    toReturn[legend] = {}
                if arg not in dict:
                    toReturn[legend].update({arg:0})
                else:
                    toReturn[legend].update({arg:dict[arg]})
        return toReturn

    def getClassDict(this):
        this.classList.sort()
        return (this.classDict, this.classList)