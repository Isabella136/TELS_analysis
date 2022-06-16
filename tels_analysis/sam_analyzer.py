from tels_analysis.sam_analysis.kegg_analyzer import kegg_analyzer
from tels_analysis.sam_analysis.megares_analyzer import megares_analyzer
from tels_analysis.sam_analysis.mge_analyzer import mge_analyzer

class sam_analyzer:
    def __init__(this, fileName, SOURCE_PREFIX, SOURCE_SUFFIX, A_TO_MEGARES, A_TO_MGES, A_TO_KEGG):
        this.megaresAnalyzer = megares_analyzer(fileName, SOURCE_PREFIX, SOURCE_SUFFIX, A_TO_MEGARES)
        this.mgeAnalyzer = mge_analyzer(fileName, SOURCE_PREFIX, SOURCE_SUFFIX, A_TO_MGES)
        this.keggAnalyzer = kegg_analyzer(fileName, SOURCE_PREFIX, SOURCE_SUFFIX, A_TO_KEGG)

    def all_genes_list(this, outputFolder, OUTPUT_PREFIX, OUTPUT_SUFFIX, SAM_ANALYSIS):
        this.megaresAnalyzer.megares_genes_list(outputFolder + "/megares/" + OUTPUT_PREFIX, OUTPUT_SUFFIX, SAM_ANALYSIS)
        this.mgeAnalyzer.mge_genes_list(outputFolder + "/mge/" + OUTPUT_PREFIX, OUTPUT_SUFFIX, SAM_ANALYSIS)
        this.keggAnalyzer.kegg_genes_list(outputFolder + "/kegg/" + OUTPUT_PREFIX, OUTPUT_SUFFIX, SAM_ANALYSIS)
