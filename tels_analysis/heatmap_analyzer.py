from tels_analysis.heatmap_analysis.double_heatmap_creator import double_heatmap_creator
from tels_analysis.heatmap_analysis import megares_analyzer
from tels_analysis import heatmap_x_axis

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
        drug_list, other_list = megares_analyzer(MEGARES)
        this.source_prefix = SOURCE_PREFIX
        this.source_suffix = SOURCE_SUFFIX
        this.amr_reads = SHORT_AMR_DIV
        this.mge_reads = SHORT_MGE
        this.drugDicts = fromListToDict(drug_list)
        this.otherDicts = fromListToDict(other_list)
        this.heatmap_list = [double_heatmap_creator("Bovine", this.drugDicts, this.otherDicts), 
                             double_heatmap_creator("Human", this.drugDicts, this.otherDicts),
                             double_heatmap_creator("Soil", this.drugDicts, this.otherDicts),
                             double_heatmap_creator("Mock", this.drugDicts, this.otherDicts)]
        
    def addToMaps(this, fileName):
        tempTuple = heatmap_x_axis(fileName)
        index = 0
        if tempTuple[0] == "Human":
            index = 1
        elif tempTuple[0] == "Soil":
            index = 2
        elif tempTuple[0] == "Mock":
            index = 3
        this.heatmap_list[index].addToMaps(tempTuple[1], this.filePath(fileName, this.amr_reads))

    def makeMaps(this, OUTPUT_PREFIX, HEATMAP):
        for i in this.heatmap_list:
            i.makeMaps(OUTPUT_PREFIX, HEATMAP)
