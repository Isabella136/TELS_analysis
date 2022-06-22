from tels_analysis.stacked_abundance_analysis.indiv_stacked_abundance import indiv_stacked_abundance
from tels_analysis import getGenesLength
from tels_analysis import getSampleGroup
from matplotlib import pyplot
import numpy, seaborn

class stacked_abundance_analyzer:
    filePath = lambda this, fileName, extension : this.source_prefix + fileName + this.source_suffix + extension

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
        this.abundance_dict = {}
        this.initial_source_size = {}

    def findAbsoluteAbundance(this, fileName):
        sample, legend = getSampleGroup(fileName)
        if sample not in this.abundance_dict:
            this.abundance_dict[sample] = indiv_stacked_abundance()
            this.initial_source_size[sample] = {legend:0}
        this.abundance_dict[sample].addToAbsolute(this.filePath(fileName, this.amr_reads), legend)
        if legend not in list(this.initial_source_size[sample].keys()):
            this.initial_source_size[sample].update({legend:0})
        this.initial_source_size[sample][legend] += float(this.allFileSizes[fileName]) / (10.0**9)
    
    def makeStack(this, outputFolder, STACKED):
        def makeAbundanceRelative():
            for sample in this.abundance_dict:
                this.abundance_dict[sample].makeAbundanceRelative(this.initial_source_size[sample], this.genes_length)
        makeAbundanceRelative()
        fig, axs = pyplot.subplots(3,8, gridspec_kw={'height_ratios': [2,.1,1.5]}, figsize=(60, 20))
        fig.suptitle('Relative Abundance & ARG Richness', fontsize=50)

        i = 0
        for sample in this.abundance_dict:
            pyplot.sca(axs[0][i])
            j = 0
            colorList = ['lightcoral', 'seagreen', 'deepskyblue', 'gold']
            abundance = this.abundance_dict[sample].getAbundance()
            for legend, dict in abundance.items():
                if j == 0:
                    axs[0][i].bar(list(dict.keys()), numpy.array(list(dict.values())), width=1.0, label=legend, color=colorList[j])
                else:
                    axs[0][i].bar(list(dict.keys()), numpy.array(list(dict.values())), width=1.0, label=legend, color=colorList[j], bottom=list(list(abundance.values())[j-1].values()))
                j += 1
            if i == 0:
                axs[0][i].set_ylabel('Log Relative Abundance', size = 25)
                pyplot.yticks(fontsize=20)
            else:
                axs[0][i].sharey(axs[0][0])
            pyplot.xticks([])
            axs[0][i].set_title(sample, size = 30)
            axs[0][i].set_anchor('NE')
            if i == 7:
                axs[0][i].legend()

            pyplot.sca(axs[1][i])
            classDict, classList = this.abundance_dict[sample].getClassDict()
            xMatrix = []
            for arg in classDict:
                xMatrix.append(classList.index(classDict[arg]))
            xMatrix = numpy.array([xMatrix])
            seaborn.heatmap(xMatrix, ax = axs[1][i], xticklabels=False, yticklabels=False, cbar=False, cmap='Purples')
            axs[1][i].set_anchor('NE')

            pyplot.sca(axs[2][i])
            labelMatrix = []
            for j in range(0,len(classList)):
                labelMatrix.append(j)
            seaborn.heatmap(numpy.array(labelMatrix).reshape(len(labelMatrix),1), ax = axs[2][i], square=True, xticklabels=False, cbar=False, cmap='Purples', linewidths=1)
            pyplot.yticks(ticks=pyplot.yticks()[0], labels=classList, fontsize=20, rotation = 0)
            axs[2][i].set_anchor('NE')
            i+=1

        pyplot.gcf()
        pyplot.savefig(outputFolder + "/arg" + STACKED)
        pyplot.close()