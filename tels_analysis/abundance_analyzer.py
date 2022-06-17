from tels_analysis.abundance_analysis.indiv_abundance import indiv_abundance
from tels_analysis import x_axis
from tels_analysis import getGenesLength
from matplotlib import pyplot
import seaborn, pandas

class abundance_analyzer:
    filePath = lambda this, fileName, extension, prefix : prefix + fileName + this.source_suffix + extension

    def __init__(this, ARG_SAM_ANALYSIS_SOURCE_PREFIX, SOURCE_SUFFIX, fileOfSizesPath, ARG_SAM_ANALYSIS, MGE_SAM_ANALYSIS, MEGARES):
        this.allFileSizes = {}
        fileOfSizes = open (fileOfSizesPath, "r")
        for line in fileOfSizes:
            this.allFileSizes.update({line.split(',')[0]:line.split(',')[1]})
        this.source_prefix = ARG_SAM_ANALYSIS_SOURCE_PREFIX
        this.source_suffix = SOURCE_SUFFIX
        this.amr_reads = ARG_SAM_ANALYSIS
        this.mge_reads = MGE_SAM_ANALYSIS
        this.genes_length = getGenesLength(MEGARES)
        this.abundance_dict = {"Bovine":{}, "Human":{}, "Soil":{}, "Mock":{}}
        this.initial_source_size = {"Bovine":{}, "Human":{}, "Soil":{}, "Mock":{}}

    def findAbsoluteAbundance(this, fileName):
        sample, x_axis_name = x_axis(fileName)
        if x_axis_name not in list(this.abundance_dict[sample].keys()):
            this.abundance_dict[sample].update({x_axis_name:indiv_abundance(sample,x_axis_name)})
            this.initial_source_size[sample].update({x_axis_name:0})
        this.abundance_dict[sample][x_axis_name].addToAbsolute(this.filePath(fileName, this.amr_reads, this.source_prefix))
        this.initial_source_size[sample][x_axis_name] += float(this.allFileSizes[fileName]) / (10.0**9)

    def makeViolinPlot(this, outputFolder, VIOLIN):
        def makeAbundanceRelative():
            for sample in this.abundance_dict:
                for x_axis_name in this.abundance_dict[sample]:
                    this.abundance_dict[sample][x_axis_name].makeAbundanceRelative(this.initial_source_size[sample][x_axis_name], this.genes_length)
        makeAbundanceRelative()
        fig, axs = pyplot.subplots(figsize=(50, 25))
        fig.suptitle("ARG Relative Abundance", fontsize=50)
        paletteList = ["#e9e5ef", "#d7c6dc", "#c7a6c8", "#ca86b6", "#d05d9d", "#c7397a", "#aa2457", "#801241",
                       "#f7f2c9", "#f2e1a8", "#ebc778", "#e6a853", "#db873a", "#c76928", "#a84d1a", "#813916",
                       "#e9f4f7", "#d6ebe8", "#b1d8cf", "#87c2b1", "#62ac8d", "#449567", "#287743", "#0d5b2c",
                       "#f2f2f2", "#dfdfdf", "#c6c6c6", "#a8a8a8", "#868686", "#686868", "#474747", "#212121"]
        x_axis_list = []
        abundance_list = {}
        for sample in this.abundance_dict:
            for x_axis_name in this.abundance_dict[sample]:
                abundance_list[x_axis_name + "_" + sample] = this.abundance_dict[sample][x_axis_name].getAbundance()
                x_axis_list.append(x_axis_name)
        df = pandas.DataFrame.from_dict(abundance_list)
        seaborn.set_style("whitegrid")
        seaborn.set_context("paper")
        ax = seaborn.violinplot(data=df, inner="box", palette=paletteList)
        ax.set_xlabel(xlabel='Samples',size=40)
        ax.set_ylabel(ylabel='Log Relative Abundance',size=40)
        ax.set_xticklabels(labels=x_axis_list,fontsize=20,rotation=90)
        pyplot.yticks(fontsize=20)
        pyplot.gcf().subplots_adjust(bottom=0.20)
        pyplot.savefig(outputFolder + "/arg" + VIOLIN)
        pyplot.close()

        