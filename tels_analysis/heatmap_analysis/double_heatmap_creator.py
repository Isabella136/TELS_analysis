from tels_analysis.heatmap_analysis.heatmap_creation.indiv_heatmap import indiv_heatmap

class double_heatmap_creator:
    def __init__(this, sample_type, drugDicts, otherDicts):
        this.sample_type = sample_type
        this.drugDicts = drugDicts
        this.otherDicts = otherDicts
        this.heatmapList = [indiv_heatmap("Drugs", this.sample_type, this.drugDicts),
                            indiv_heatmap("Metals/Biocides", this.sample_type, this.otherDicts)]
    def addToMaps(this, x_axis, filepath, drugBool, otherBool):
        drugBool = this.heatmapList[0].addToMap(x_axis, filepath, drugBool)
        otherBool = this.heatmapList[1].addToMap(x_axis, filepath, otherBool)
        return (drugBool, otherBool)

    def makeMaps(this, drugBool, otherBool, drugClassList, otherClassList):
        drugMatrix, drugXaxis = this.heatmapList[0].makeMap(drugBool, drugClassList)
        otherMatrix, otherXaxis = this.heatmapList[1].makeMap(otherBool, otherClassList)
        return (drugMatrix, drugXaxis, otherMatrix, otherXaxis)
