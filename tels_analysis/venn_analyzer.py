from tels_analysis.venn_analysis.indiv_venn import indiv_venn
from tels_analysis import getSampleAndIndex

class venn_analyzer:

    filePath = lambda this, fileName, extension : this.source_prefix + fileName + this.source_suffix + extension

    def __init__(this, ARG_SAM_ANALYSIS_SOURCE_PREFIX, SOURCE_SUFFIX, ARG_SAM_ANALYSIS, MGE_SAM_ANALYSIS):
        this.source_prefix = ARG_SAM_ANALYSIS_SOURCE_PREFIX
        this.source_suffix = SOURCE_SUFFIX
        this.amr_reads = ARG_SAM_ANALYSIS
        this.mge_reads = MGE_SAM_ANALYSIS
        this.venn_dict = {}

    def addToCount(this, fileName):
        sample, index = getSampleAndIndex(fileName)
        if sample not in this.venn_dict:
            this.venn_dict[sample] = indiv_venn(sample)
        this.venn_dict[sample].addToCount(this.filePath(fileName, this.amr_reads), index)

    def makeVenn(this, outputFolder, VENN):
        for sample in this.venn_dict:
            this.venn_dict[sample].findFinalCount()
            this.venn_dict[sample].makeFigure(outputFolder, VENN)