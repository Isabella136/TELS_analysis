import numpy

class indiv_heatmap:
    def __init__(this, antimicrobial, sample_type, mechanismDicts):
        this.antimicrobial = antimicrobial
        this.sample_type = sample_type
        this.mechanism_dict = mechanismDicts[1]
        this.mech_list = list(this.mechanism_dict.keys())
        this.columns = [[0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict),
                        [0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict)]
        this.x_axis_list = []
        
    def addToMap(this, x_axis, filepath, dictBool):
        if this.x_axis_list.count(x_axis) == 0:
            this.x_axis_list.append(x_axis)
        argFile = open(filepath, "r")
        lineIndex = 0
        for line in argFile:
            lineIndex += 1
            if lineIndex < 20:
                continue
            splitLine = line.split('|')
            if this.antimicrobial == "Drugs":
                if splitLine[1] != "Drugs":
                    continue
            else:
                if splitLine[1] == "Drugs":
                    continue
            if len(splitLine) == 6:
                continue
            index = this.mech_list.index(splitLine[3])
            this.columns[len(this.x_axis_list)-1][index] = 1
            cl = splitLine[2]
            if cl == "betalactams": cl = "Betalactams"
            dictBool[(cl,splitLine[3])] = True
        argFile.close()
        return dictBool
    
    def makeMap(this, dictBool, classList):
        for tuple in dictBool:
            if not(dictBool[tuple]):
                index = this.mech_list.index(tuple[1])
                for column in this.columns:
                    column.pop(index)
                this.mech_list.pop(index)
        for column in this.columns:
            for i in range(0, len(column)):
                if column[i] == 1:
                    heatmapVal = classList.index(this.mechanism_dict[this.mech_list[i]]) + 1
                    column[i] = heatmapVal
        vals = numpy.array(this.columns).transpose()
        return vals, this.x_axis_list