from tels_analysis.stacked_abundance_analysis.indiv_stacked_abundance import indiv_stacked_abundance
from tels_analysis import getSampleGroup
from matplotlib import pyplot
import numpy, seaborn, os

class stacked_abundance_analyzer:
    filePath = lambda this, prefix, fileName, extension : prefix + fileName + this.source_suffix + extension

    def __init__(this, SOURCE_PREFIX, SOURCE_SUFFIX, SHORT_AMR_DIV, SHORT_MGE, STATS):
        this.source_prefix = SOURCE_PREFIX
        this.source_suffix = SOURCE_SUFFIX
        this.amr_reads = SHORT_AMR_DIV
        this.mge_reads = SHORT_MGE
        this.stats = STATS
        this.abundance_dict = {}

    def findAbsoluteAbundance(this, fileName):
        sample, legend = getSampleGroup(fileName)
        if sample not in this.abundance_dict:
            this.abundance_dict[sample] = indiv_stacked_abundance()
        this.abundance_dict[sample].addToAbsolute(this.filePath(this.source_prefix, fileName, this.amr_reads), legend, this.filePath(this.source_prefix, fileName, this.stats))
    
    def makeStack(this, outputFolder, STACKED):
        def makeAbundanceRelative():
            for sample in this.abundance_dict:
                this.abundance_dict[sample].makeAbundanceRelative()
        makeAbundanceRelative()
        fig, axs = pyplot.subplots(3,8, gridspec_kw={'height_ratios': [2,.1,1.5]}, figsize=(60, 20))
        fig.suptitle('Relative Abundance & ARG Richness', fontsize=50)
        i = 0
        for sample in this.abundance_dict:
            pyplot.sca(axs[0][i])
            j = 0
            colorList = ['red', 'forestgreen', 'royalblue', 'gold']
            abundance = this.abundance_dict[sample].getAbundance()
            bottom_array = None
            for legend, dict in abundance.items():
                current_array = numpy.array(list(dict.values()))
                if j == 0:
                    axs[0][i].bar(list(dict.keys()), current_array, width=1.0, label=legend, color=colorList[j], alpha = 0.5)
                    bottom_array = current_array
                    axs.sort
                else:
                    axs[0][i].bar(list(dict.keys()), current_array, width=1.0, label=legend, color=colorList[j], alpha = 0.5, bottom=bottom_array)
                    bottom_array = bottom_array + current_array
                j += 1
            if i == 0:
                axs[0][i].set_ylabel('Log Relative Abundance', size = 25)
                pyplot.yticks(fontsize=20)
            else:
                axs[0][i].sharey(axs[0][0])
            pyplot.xticks([])
            pyplot.margins(x=0)
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
            seaborn.heatmap(xMatrix, ax = axs[1][i], xticklabels=False, yticklabels=False, cbar=False, cmap='gist_rainbow')
            axs[1][i].set_anchor('NE')

            pyplot.sca(axs[2][i])
            labelMatrix = []
            for j in range(0,len(classList)):
                labelMatrix.append(j)
            seaborn.heatmap(numpy.array(labelMatrix).reshape(len(labelMatrix),1), ax = axs[2][i], square=True, xticklabels=False, cbar=False, cmap='gist_rainbow', linewidths=1)
            pyplot.yticks(ticks=pyplot.yticks()[0], labels=classList, fontsize=20, rotation = 0)
            axs[2][i].set_anchor('NE')
            i+=1

        if not(os.path.exists(outputFolder)):
            os.makedirs(outputFolder)
        pyplot.gcf()
        pyplot.savefig(outputFolder + "/arg" + STACKED)
        pyplot.close()