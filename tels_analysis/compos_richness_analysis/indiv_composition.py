from tels_analysis import fileDict
from matplotlib import pyplot
from PIL import Image
import numpy


class indiv_composition:
    def __init__(this, fileName, mge_dict):
        this.name = fileName
        this.sample, this.seqPlatform = fileDict(fileName)
        this.arg_composition_data = [0,0,0,0] #in order: drug, metal, multi-compound, and biocide
        this.arg_class_richness = 0
        this.arg_mechanism_richness = 0
        this.arg_group_richness = 0
        this.mge_composition_data = [0,0,0,0,0,0,0] #in order: PLASMID, PHAGE, TE, IS, ICE, VIRUS, UNCLASSIFIED
        this.mge_richness = 0
        this.mge_dict = mge_dict

    def getData(this):
        toReturn = this.sample
        toReturn = toReturn + ',' + this.seqPlatform
        toReturn = toReturn + ',,' + str(this.arg_class_richness)
        toReturn = toReturn + ',' + str(this.arg_mechanism_richness)
        toReturn = toReturn + ',' + str(this.arg_group_richness)
        toReturn = toReturn + ',,' + str(this.mge_richness)
        return toReturn

    def findAllData(this, filepath_ARG_composition, filepath_MGE, filepath_arg_output, filepath_mge_output):
        this.findARGComposition(filepath_ARG_composition)
        this.findMGEComposition(filepath_MGE)
        this.makePieChart(filepath_arg_output, filepath_mge_output)

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
        mgeList = []
        for line in mge_composition_file:
            lineNum += 1
            if lineNum <= 20:
                continue
            mge = line.split(',')[0]
            if mge not in this.mge_dict:
                print(mge)
                continue
            if mgeList.count(mge) == 0:
                mgeList.append(mge)
            if this.mge_dict[mge] == "PLASMID":
                this.mge_composition_data[0] += 1
            elif this.mge_dict[mge] == "PROPHAGE":
                this.mge_composition_data[1] += 1
            elif this.mge_dict[mge] == "TE":
                this.mge_composition_data[2] += 1
            elif this.mge_dict[mge] == "IS":
                this.mge_composition_data[3] += 1
            elif this.mge_dict[mge] == "ICE":
                this.mge_composition_data[4] += 1
            elif this.mge_dict[mge] == "VIRUS":
                this.mge_composition_data[5] += 1
            elif this.mge_dict[mge] == "UNCLASSIFIED":
                this.mge_composition_data[6] += 1
        this.mge_richness = len(mgeList)
        mge_composition_file.close()

    def makePieChart(this, filepath_arg_output, filepath_mge_output):
        if this.arg_class_richness == 0:
            img = Image.new(mode = "RGB", size = (120, 48), color = (255, 255, 255))
            img.save(filepath_arg_output)
        else:
            vals = numpy.array(this.arg_composition_data)
            chartColors = ["#5891AD", "#004561", "#FF6F31", "#1C7685"]
            pyplot.figure(figsize=(120,48))
            pyplot.pie(vals, colors=chartColors, startangle=90)
            pyplot.savefig(filepath_arg_output, dpi=10)
            pyplot.close()

        if this.mge_richness == 0:
            img = Image.new(mode = "RGB", size = (120, 48), color = (255, 255, 255))
            img.save(filepath_mge_output)
        else:
            vals = numpy.array(this.mge_composition_data)
            chartColors = ["#34A853", "#FBBC04", "#FF6D01", "#EA4335", "#4285F4", "#46BDC6", "#FFC0CB"]
            pyplot.figure(figsize=(120,48))
            pyplot.pie(vals, colors=chartColors, startangle=90)
            pyplot.savefig(filepath_mge_output, dpi=10)
            pyplot.close()