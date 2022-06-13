from tels_analysis.heatmap_analysis.double_heatmap_creator import double_heatmap_creator
from tels_analysis.heatmap_analysis import megares_analyzer
from tels_analysis import heatmap_x_axis
from matplotlib import pyplot
import seaborn, numpy

class heatmap_analyzer:
    filePath = lambda this, fileName, extension : this.source_prefix + fileName + this.source_suffix + extension

    def __init__(this, SOURCE_PREFIX, SOURCE_SUFFIX, SHORT_AMR_DIV, SHORT_MGE, MEGARES):

        def fromListToDict(mechanismList):
            classDict = {}
            mechanismDict = {}
            for tuple in mechanismList:
                mechanismDict.update({tuple[1]:tuple[0]})
                if classDict.get(tuple[0], 0) == 0:
                    classDict.update({tuple[0]:1})
                else:
                    classDict[tuple[0]] += 1
            return (classDict, mechanismDict)
        def makeValsBools(mechanismList):
            toReturn = {}
            for tuple in mechanismList:
                toReturn.update({tuple:False})
            return toReturn
        drug_list, other_list = megares_analyzer(MEGARES)
        this.source_prefix = SOURCE_PREFIX
        this.source_suffix = SOURCE_SUFFIX
        this.amr_reads = SHORT_AMR_DIV
        this.mge_reads = SHORT_MGE
        this.drugClassDict, drugMechDict = fromListToDict(drug_list)
        this.otherClassDict, otherMechDict = fromListToDict(other_list)
        this.drugBool = makeValsBools(drug_list)
        this.otherBool = makeValsBools(other_list)
        this.heatmap_list = [double_heatmap_creator("Bovine", (this.drugClassDict, drugMechDict), (this.otherClassDict, otherMechDict)), 
                             double_heatmap_creator("Human", (this.drugClassDict, drugMechDict), (this.otherClassDict, otherMechDict)),
                             double_heatmap_creator("Soil", (this.drugClassDict, drugMechDict), (this.otherClassDict, otherMechDict)),
                             double_heatmap_creator("Mock", (this.drugClassDict, drugMechDict), (this.otherClassDict, otherMechDict))]
        
    def addToMaps(this, fileName):
        tempTuple = heatmap_x_axis(fileName)
        index = 0
        if tempTuple[0] == "Human":
            index = 1
        elif tempTuple[0] == "Soil":
            index = 2
        elif tempTuple[0] == "Mock":
            index = 3
        this.drugBool, this.otherBool = this.heatmap_list[index].addToMaps(tempTuple[1], this.filePath(fileName, this.amr_reads), this.drugBool, this.otherBool)

    def makeMaps(this, OUTPUT_PREFIX, HEATMAP):
        for tuple in this.drugBool:
            if not(this.drugBool[tuple]):
                this.drugClassDict[tuple[0]] -= 1
                if this.drugClassDict[tuple[0]] == 0:
                    this.drugClassDict.pop(tuple[0])
        for tuple in this.otherBool:
            if not(this.otherBool[tuple]):
                this.otherClassDict[tuple[0]] -= 1
                if this.otherClassDict[tuple[0]] == 0:
                    this.otherClassDict.pop(tuple[0])
        drugsMatrixList = []
        otherMatrixList = []
        drugsXaxisList = []
        otherXaxisList = []
        for i in this.heatmap_list:
            drugMatrix, drugXaxis, otherMatrix, otherXaxis = i.makeMaps(this.drugBool, this.otherBool, list(this.drugClassDict.keys()), list(this.otherClassDict.keys()))
            drugsMatrixList.append(drugMatrix)
            otherMatrixList.append(otherMatrix)
            drugsXaxisList.append(drugXaxis)
            otherXaxisList.append(otherXaxis)


        def heatmapMaker(matrixList, xaxisList, classDict, type):
            labelMatrix = []
            labelPos = []
            index = 1
            mechTotal = 0
            for cl in classDict:
                for i in range(0,classDict[cl]):
                    labelMatrix.append(index)
                index += 1
                down = mechTotal
                up = down + classDict[cl]
                labelPos.append((up + down - 1)/2)
                mechTotal += classDict[cl]

            fig, axs = pyplot.subplots(1,5, gridspec_kw={'width_ratios': [1,7,7,7,7]}, figsize=(40, 25))
            fig.suptitle('ARG - ' + type, fontsize=50)

            seaborn.heatmap(numpy.array(labelMatrix).reshape(len(labelMatrix),1), ax = axs[0], xticklabels=False, cbar=False, cmap="viridis", vmin=0, vmax=numpy.max(labelMatrix)+1)
            pyplot.sca(axs[0])
            pyplot.yticks(labelPos, list(classDict.keys()), fontsize=20, rotation = 0)

            seaborn.heatmap(matrixList[0], ax=axs[1], xticklabels=xaxisList[0], yticklabels=False, cbar=False, cmap="viridis", vmin=0, vmax=numpy.max(labelMatrix)+1)
            pyplot.sca(axs[1])
            pyplot.xticks(fontsize=20)
            axs[1].set_title("Bovine", fontsize=40)

            seaborn.heatmap(matrixList[1], ax=axs[2], xticklabels=xaxisList[1], yticklabels=False, cbar=False, cmap="viridis", vmin=0, vmax=numpy.max(labelMatrix)+1)
            pyplot.sca(axs[2])
            pyplot.xticks(fontsize=20)
            axs[2].set_title("Human", fontsize=40)

            seaborn.heatmap(matrixList[2], ax=axs[3], xticklabels=xaxisList[2], yticklabels=False, cbar=False, cmap="viridis", vmin=0, vmax=numpy.max(labelMatrix)+1)
            pyplot.sca(axs[3])
            pyplot.xticks(fontsize=20)
            axs[3].set_title("Soil", fontsize=40)

            seaborn.heatmap(matrixList[3], ax=axs[4], xticklabels=xaxisList[3], yticklabels=False, cbar=False, cmap="viridis", vmin=0, vmax=numpy.max(labelMatrix)+1)
            pyplot.sca(axs[4])
            pyplot.xticks(fontsize=20)
            axs[4].set_title("Mock", fontsize=40)

            pyplot.gcf().subplots_adjust(bottom=0.20, left=0.20)
            pyplot.savefig(OUTPUT_PREFIX + "/" + type + HEATMAP)
            pyplot.close()

        heatmapMaker(drugsMatrixList, drugsXaxisList, this.drugClassDict, "Drugs")
        heatmapMaker(otherMatrixList, otherXaxisList, this.otherClassDict, "Metals and Biocides")