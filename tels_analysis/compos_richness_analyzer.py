from tels_analysis.compos_richness_analysis.indiv_composition import indiv_composition

class compos_richness_analyzer:
    filePath = lambda this, fileName, extension : this.source_prefix + fileName + this.source_suffix + extension

    def __init__(this, SOURCE_PREFIX, SOURCE_SUFFIX, SHORT_AMR_DIV, SHORT_MGE):
        this.source_prefix = SOURCE_PREFIX
        this.source_suffix = SOURCE_SUFFIX
        this.amr_reads = SHORT_AMR_DIV
        this.mge_reads = SHORT_MGE
        this.compos_richness_list = []
        this.MGE_list = []

    def analyzeFile(this, fileName, outputFolder, INDIV_COMPOS_CHART):
        fileComp = indiv_composition(fileName)
        fileComp.findAllData(this.filePath(fileName, this.amr_reads), this.filePath(fileName, this.mge_reads), outputFolder + "/" + fileName + INDIV_COMPOS_CHART)
        this.compos_richness_list.append(fileComp.getData())

        mge_file = open(this.filePath(fileName, this.mge_reads), "r")
        lineIndex = 0
        for line in mge_file:
            lineIndex += 1
            if lineIndex < 21:
                continue
            mge = line.split(',')[0]
            if this.MGE_list.count(mge) == 0:
                this.MGE_list.append(mge)

    def printAnalysis(this, outputFolder, COMPOS_RICHNESS_ANALYSIS):
        analysis = open(outputFolder + "/" + COMPOS_RICHNESS_ANALYSIS, "w")
        analysis.write("Sample,Sequencing platform,ARG composition,ARG Class richness,ARG Mechanism richness,ARG Group richness,MGE composition,MGE Accession richness\n")
        for c_r in this.compos_richness_list:
            analysis.write(c_r + "\n")
        analysis.close()

        analysis = open(outputFolder + "/MGE_list.csv", "w")
        for MGE in this.MGE_list:
            analysis.write(MGE + "\n")
        analysis.close()
