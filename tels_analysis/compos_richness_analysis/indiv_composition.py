from tels_analysis import fileDict
from matplotlib import pyplot
from PIL import Image
import numpy


class indiv_composition:
    def __init__(this, fileName):
        this.name = fileName
        this.sample, this.seqPlatform = fileDict(fileName)
        this.arg_composition_data = [0,0,0,0] #in order: drug, metal, multi-compound, and biocide
        this.arg_class_richness = 0
        this.arg_mechanism_richness = 0
        this.arg_group_richness = 0
        this.mge_composition_data = [0,0,0,0,0,0] #in order: plasmid, phage, TE, IS, ICE, virus
        this.mge_richness = 0

    def getData(this):
        toReturn = this.sample
        toReturn = toReturn + ',' + this.seqPlatform
        toReturn = toReturn + ',,' + str(this.arg_class_richness)
        toReturn = toReturn + ',' + str(this.arg_mechanism_richness)
        toReturn = toReturn + ',' + str(this.arg_group_richness)
        toReturn = toReturn + ',,' + str(this.mge_richness)
        return toReturn

    def findAllData(this, filepath_ARG_composition, filepath_MGE, filepath_output):
        this.findARGComposition(filepath_ARG_composition)
        this.findMGEComposition(filepath_MGE)
        if this.arg_class_richness == 0:
            img = Image.new(mode = "RGB", size = (120, 48), color = (255, 255, 255))
            img.save(filepath_output)
        else:
            this.makePieChart(filepath_output)

    def findARGComposition(this, filepath):
        arg_composition_file = open(filepath, "r")
        lineNum = 0
        classList = []
        mechanismList = []
        groupList= []
        for line in arg_composition_file:
            lineNum += 1
            if lineNum < 20:
                continue
            arg_header = line.split(',')[0].split('|')
            if len(arg_header) == 6: continue
            if arg_header[1] == "Drugs":
                this.arg_composition_data[0] += 1
            elif arg_header[1] == "Metals":
                this.arg_composition_data[1] += 1
            elif arg_header[1] == "Multi-compound":
                this.arg_composition_data[2] += 1
            else: #arg_type == "Biocides"
                this.arg_composition_data[3] += 1
            if classList.count(arg_header[2]) == 0:
                classList.append(arg_header[2])
            if mechanismList.count(arg_header[3]) == 0:
                mechanismList.append(arg_header[3])
            if groupList.count(arg_header[4]) == 0:
                groupList.append(arg_header[4])
        arg_composition_file.close()
        this.arg_class_richness = len(classList)
        this.arg_group_richness = len(groupList)
        this.arg_mechanism_richness = len(mechanismList)

    def findMGEComposition(this, filepath):
        mge_composition_file = open(filepath)
        lineNum = 0
        for line in mge_composition_file:
            lineNum += 1
            if lineNum < 19:
                continue
            else:
                this.mge_richness = int(line.split(',')[1][:-1])
                break
        mge_composition_file.close()

    def makePieChart(this, filepath):
        vals = numpy.array(this.arg_composition_data)
        chartColors = ["#5891AD", "#004561", "#FF6F31", "#1C7685"]
        pyplot.figure(figsize=(120,48))
        pyplot.pie(vals, colors=chartColors, startangle=90)
        pyplot.savefig(filepath, dpi=10)
        pyplot.close()