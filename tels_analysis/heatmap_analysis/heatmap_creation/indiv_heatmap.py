from re import split
from matplotlib import pyplot
import seaborn, numpy

class indiv_heatmap:
    def __init__(this, antimicrobial, sample_type, mechanismDicts):
        this.antimicrobial = antimicrobial
        this.sample_type = sample_type
        this.mechanism_dict = mechanismDicts[1]
        this.class_dict = mechanismDicts[0]
        this.mech_list = list(this.mechanism_dict.keys())
        this.class_list = list(this.class_dict.keys())
        this.columns = [[0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict),
                        [0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict),[0]*len(this.mechanism_dict)]
        this.x_axis_list = []
        
    def addToMap(this, x_axis, filepath):
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
            index = this.mech_list.index(splitLine[3])
            heatmapVal = this.class_list.index(splitLine[2])
            this.columns[len(this.x_axis_list)-1][index] = heatmapVal
        argFile.close()
    
    def makeMap(this, filepath):
        vals = numpy.array(this.columns)
        seaborn.heatmap(vals, yticklabels=this.x_axis_list, xticklabels=False, cbar=False, cmap="viridis")
        pyplot.yticks(rotation=45)
        pyplot.savefig(filepath, dpi=100)
        pyplot.close()