from tels_analysis.abundance_analysis.indiv_abundance import IndivAbundance
from tels_analysis import get_genes_length
from tels_analysis import x_axis
from matplotlib import pyplot
import seaborn
import pandas
import csv
import os

class AbundanceAnalyzer:

    # lambda that generates tels output file
    file_path = (lambda self, sample_name, extension : self.source_prefix + sample_name + self.source_suffix + extension)

    def __init__(self, SOURCE_PREFIX, SOURCE_SUFFIX, FILE_SIZE_OUTPUT,
                 AMR_DIV_EXT, MGE_EXT, ACLAME_DB, ICEBERG_DB, PLASMID_DB, MEGARES_DB):
        
        # Retrieve size of dedup.fastq files
        self.all_file_sizes = dict()
        with open(FILE_SIZE_OUTPUT, "r") as file_of_sizes:
            csvreader = csv.reader(file_of_sizes)
            for row in csvreader:
                self.all_file_sizes.update({row[0]:row[1]})

        # Get tels output file names information
        self.source_prefix = SOURCE_PREFIX
        self.source_suffix = SOURCE_SUFFIX
        self.amr_reads_ext = AMR_DIV_EXT
        self.mge_reads_ext = MGE_EXT

        # Get gene lengths for mges and megares
        self.mge_genes_length = get_genes_length(ACLAME_DB)
        self.mge_genes_length.update(get_genes_length(ICEBERG_DB))
        self.mge_genes_length.update(get_genes_length(PLASMID_DB))

        self.megares_genes_length = get_genes_length(MEGARES_DB)

        # Estabilsh dictionaries for absolute abundance        
        self.abundance_dict_amr = {"Bovine":{}, "Human":{}, "Soil":{}, "Mock":{}}
        self.abundance_dict_mge = {"Bovine":{}, "Human":{}, "Soil":{}, "Mock":{}}
        self.initial_source_size = {"Bovine":{}, "Human":{}, "Soil":{}, "Mock":{}}

    def find_absolute_abundance(self, sample_name, amr_analysis, mge_analysis):
        # Get sample source (Bovine, Human, Soil, or Mock)
        # + probe type (PacBio, XT, or V2 and ARG, ARG-MGE, or MGE)
        sample, probe = x_axis(sample_name)

        # If we are analyzing AMRs
        if amr_analysis:

            # This is the first time we see this probe for the sample for amr
            if probe not in list(self.abundance_dict_amr[sample].keys()):
                self.abundance_dict_amr[sample].update({probe:IndivAbundance(True)})

            self.abundance_dict_amr[sample][probe].add_to_absolute(
                self.file_path(sample_name, self.amr_reads_ext))

        # If we are analyzing MGEs
        if mge_analysis:

            # This is the first time we see this probe for the sample for mge
            if probe not in list(self.abundance_dict_mge[sample].keys()):
                self.abundance_dict_mge[sample].update({probe:IndivAbundance(False)})

            self.abundance_dict_mge[sample][probe].add_to_absolute(
                self.file_path(sample_name, self.mge_reads_ext))
            
        # Overall, tThis is the first time we see this probe for the sample
        if probe not in list(self.initial_source_size[sample].keys()):
            self.initial_source_size[sample].update({probe:0}) 

        self.initial_source_size[sample][probe] += (
            float(self.all_file_sizes[sample_name]) / (10.0**9))

    def make_violin_plot(self, output_folder, violin_ext, amr_analysis, mge_analysis):
        def violin_plot_per_analysis(abundance_dict, genes_length, element_name):

            # Go through absolute abundace information to make it relative
            for sample in abundance_dict:
                for probe in abundance_dict[sample]:
                    abundance_dict[sample][probe].make_abundance_relative(
                        self.initial_source_size[sample][probe], genes_length)

            # Set up main figure        
            fig, axs = pyplot.subplots(2,2,figsize=(70, 70),sharey="row")
            fig.suptitle(element_name + " Relative Abundance", fontsize=70)
            fig.add_subplot(111, frame_on=False)
            pyplot.grid(False)
            pyplot.tick_params(labelcolor="none", bottom=False, left=False)
            for ax in axs.flat:
                ax.set_ylabel ("Log Relative Abundance", fontsize=55)
                ax.label_outer()
                ax.set_ylim(-4,6)

            palette_list = ["PuRd", "YlOrBr", "BuGn", "Greys"]

            # Set up individual violin plot
            for index, sample in enumerate(abundance_dict):
                pyplot.sca(axs[index// 2, index% 2])
                x_axis_list = list()
                sample_abundance_dict = dict()
                for probe in abundance_dict[sample]:
                    sample_abundance_dict[probe] = (
                        abundance_dict[sample][probe].get_abundance())
                    x_axis_list.append(probe)
                df = pandas.DataFrame.from_dict(sample_abundance_dict)
                seaborn.set_style("whitegrid")
                seaborn.set_context("paper")

                axs[index// 2, index% 2] = seaborn.violinplot(
                        data=df, inner="box", palette=palette_list[index])
                axs[index// 2, index% 2].grid(False)
                axs[index// 2, index% 2].set_title(sample, fontsize=60)

                axs[index// 2, index% 2].set_xticklabels(
                    labels=x_axis_list,
                    fontsize=40,
                    rotation=90)
                
                if index% 2 == 0:
                    pyplot.yticks(fontsize=40)

            if not(os.path.exists(output_folder)):
                os.makedirs(output_folder)
            pyplot.gcf().subplots_adjust(bottom=0.20)
            pyplot.savefig(output_folder + "/" + element_name + violin_ext)
            pyplot.close()

        if amr_analysis:
            violin_plot_per_analysis(self.abundance_dict_amr, 
                                     self.megares_genes_length, 
                                     'ARG')           
        if mge_analysis:
            violin_plot_per_analysis(self.abundance_dict_mge, 
                                     self.mge_genes_length, 
                                     'MGE')