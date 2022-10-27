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
        lineNum = 0
        for line in argFile:
            lineNum += 1
            if lineNum < 20:
                continue
            line_list = line.split(',')
            arg_header = line_list[0].split("|")
            if arg_header[4] not in list(this.absolute_abundance[legend].keys()):
                this.absolute_abundance[legend][arg_header[4]] = 0
            this.absolute_abundance[legend][arg_header[4]] += int(line_list[1])
            className = arg_header[2]
            if arg_header[1] != "Drugs":
                className = arg_header[1] + " resistance"
            elif className == "betalactams":
                className = "Betalactams"
            elif className == "Mycobacterium_tuberculosis-specific_Drug":
                className = "M_tuberculosis-specific_Drug"
            if arg_header[4] not in this.classDict:
                this.classDict[arg_header[4]] = className
            if className not in this.classList:
                this.classList.append(className)
    
    def makeAbundanceRelative(this):
        avgReads = {}
        for legend, list in this.statsFilepaths.items():
            for filepath in list:
                statFile = open(filepath, "r")
                readCount = int(statFile.readline().split(',')[1])
                if legend not in avgReads:
                    avgReads[legend] = 0
                avgReads[legend]+=readCount
                statFile.close()
            avgReads[legend] = avgReads[legend]/3
        total = {}
        tempRelative = {}
        for legend, dict in this.absolute_abundance.items():
            for group, count in dict.items():
                if legend not in tempRelative:
                    tempRelative[legend] = {}
                tempRelative[legend].update({group:log10(count/avgReads[legend]*1000000)})
                if group not in total:
                    total[group] = 0
                total[group] += tempRelative[legend][group]
        sortedTotal = {}
        while len(total) > 0:
            max = ('key', -1)
            for key,val in total.items():
                if val > max[1]:
                    max = (key,val)
            sortedTotal.update({max[0]:max[1]})
            total.pop(max[0])
        classDictTemp = this.classDict
        this.classDict = {}
        for key in sortedTotal:
            this.classDict[key] = classDictTemp[key]
            for legend in tempRelative:
                if legend not in this.relative_abundance:
                    this.relative_abundance[legend] = {}
                if key not in tempRelative[legend]:
                    this.relative_abundance[legend].update({key:0})
                else:
                    this.relative_abundance[legend].update({key:tempRelative[legend][key]})
        classDictTemp.clear()
        sortedTotal.clear()
        

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