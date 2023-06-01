from tels_analysis.abundance_analysis.indiv_abundance import indiv_abundance
from tels_analysis import x_axis
from tels_analysis import getGenesLength
from matplotlib import pyplot
import seaborn, pandas, os

class abundance_analyzer:
    filePath = lambda this, fileName, extension, prefix : prefix + fileName + this.source_suffix + extension

    def __init__(this, SOURCE_PREFIX, SOURCE_SUFFIX, fileOfSizesPath, SHORT_AMR_DIV, SHORT_MGE, MEGARES):
        this.allFileSizes = {}
        fileOfSizes = open (fileOfSizesPath, "r")
        for line in fileOfSizes:
            this.allFileSizes.update({line.split(',')[0]:line.split(',')[1]})
        this.source_prefix = SOURCE_PREFIX
        this.source_suffix = SOURCE_SUFFIX
        this.amr_reads = SHORT_AMR_DIV
        this.mge_reads = SHORT_MGE
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
        fig, axs = pyplot.subplots(2,2,figsize=(70, 70),sharey="row")
        fig.suptitle("ARG Relative Abundance", fontsize=70)
        fig.add_subplot(111, frame_on=False)
        
        pyplot.tick_params(labelcolor="none", bottom=False, left=False)
        for ax in axs.flat:
            ax.set_ylabel ("Log Relative Abundance", fontsize=55)
            ax.label_outer()
            ax.set_ylim(-4,6)
        paletteList = ["PuRd", "YlOrBr", "BuGn", "Greys"]
        i = 0
        for sample in this.abundance_dict:
            pyplot.sca(axs[i // 2, i % 2])
            x_axis_list = []
            abundance_list = {}
            for x_axis_name in this.abundance_dict[sample]:
                abundance_list[x_axis_name] = this.abundance_dict[sample][x_axis_name].getAbundance()
                x_axis_list.append(x_axis_name)
            df = pandas.DataFrame.from_dict(abundance_list)
            seaborn.set_style("whitegrid")
            seaborn.set_context("paper")
            axs[i // 2, i % 2] = seaborn.violinplot(data=df, inner="box", palette=paletteList[i])
            axs[i // 2, i % 2].set_title(sample, fontsize=60)
            axs[i // 2, i % 2].set_xticklabels(labels=x_axis_list,fontsize=40,rotation=90)
            if i % 2 == 0:
                pyplot.yticks(fontsize=40)
            i+=1
        if not(os.path.exists(outputFolder)):
            os.makedirs(outputFolder)
        pyplot.gcf().subplots_adjust(bottom=0.20)
        pyplot.savefig(outputFolder + "/arg" + VIOLIN)
        pyplot.close()

        