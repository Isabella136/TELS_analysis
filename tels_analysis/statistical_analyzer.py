from tels_analysis.statistical_analysis.indiv_stats import indiv_stats

class statistical_analyzer:
    filePath = lambda this, fileName, extension : this.source_prefix + fileName + this.source_suffix + extension

    def __init__(this, SOURCE_PREFIX, SOURCE_SUFFIX, ARG_SAM_ANALYSIS, MGE_SAM_ANALYSIS, STATS):
        this.source_prefix = SOURCE_PREFIX
        this.source_suffix = SOURCE_SUFFIX
        this.amr_reads = ARG_SAM_ANALYSIS
        this.mge_reads = MGE_SAM_ANALYSIS
        this.stats = STATS
        this.statList = []

    def analyzeFile(this, fileName):
        fileStats = indiv_stats(fileName)
        fileStats.findAllStats(this.filePath(fileName, this.stats), this.filePath(fileName, this.amr_reads), this.filePath(fileName, this.mge_reads))
        this.statList.append(fileStats.getStats())

    def printAnalysis(this, outputFolder, STATISTICAL_ANALYSIS):
        analysis = open(outputFolder + "/" + STATISTICAL_ANALYSIS, "w")
        analysis.write("Sample,Sequencing platform,Raw reads ,De-duplicated reads,Duplication (%),ARG On-target (%),MGE On-target (%)\n")
        for stat in this.statList:
            analysis.write(stat + "\n")
        analysis.close()