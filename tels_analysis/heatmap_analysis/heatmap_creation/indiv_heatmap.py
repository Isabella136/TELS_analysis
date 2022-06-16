import numpy

class indiv_heatmap:
    def __init__(this, antimicrobial, sample_type, mechanismDicts):
        this.antimicrobial = antimicrobial
        this.sample_type = sample_type
        this.mechanism_dict = mechanismDicts[1]
        this.mech_list = list(this.mechanism_dict.keys())
        this.columns = [[0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict),
                        [0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict)]
        this.x_axis_list = []
        
    def addToMap(this, x_axis, filepath, dictBool):
        if this.x_axis_list.count(x_axis) == 0:
            this.x_axis_list.append(x_axis)
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
                arg_header = arg.split('|')
                if this.antimicrobial == "Drugs":
                    if arg_header[1] != "Drugs":
                        continue
                else:
                    if arg_header[1] == "Drugs":
                        continue
                index = this.mech_list.index(arg_header[3])
                this.columns[len(this.x_axis_list)-1][index] = 1
                cl = arg_header[2]
                if cl == "betalactams": cl = "Betalactams"
                dictBool[(cl,arg_header[3])] = True
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