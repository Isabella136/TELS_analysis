from tels_analysis.abundance_analysis.indiv_abundance import IndivAbundance
from tels_analysis import get_sample_name_definition
from tels_analysis import get_mge_annot_dict
from tels_analysis import get_genes_length
from tels_analysis import tels_file_path
from matplotlib import pyplot
import seaborn
import pandas
import csv
import os

class AbundanceAnalyzer:
    
    def __init__(
            self, SOURCE_PREFIX, SOURCE_SUFFIX, FILE_SIZE_OUTPUT, AMR_DIV_EXT, 
            MGE_EXT, ACLAME_DB, ICEBERG_DB, PLASMID_DB, MEGARES_DB, MGES_ANNOTATION):
        
        # Retrieve size of dedup.fastq files
        self.all_file_sizes = dict()
        with open(FILE_SIZE_OUTPUT, "r") as file_of_sizes:
            csvreader = csv.reader(file_of_sizes)
            for row in csvreader:
                self.all_file_sizes.update({row[0]:row[1]})

        # Get MGE types
        self.mges_annot = get_mge_annot_dict(MGES_ANNOTATION)

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
        # Retrieve definition in tuple form Organism, Platform, Chemistry, Probe)
        # and save the organism (Bovine, Human, Soil, or Mock)
        # + probe type (PacBio, XT, or V2 and ARG, ARG-MGE, or MGE)
        sample_name_definition = get_sample_name_definition(sample_name)
        organism = sample_name_definition[0]
        if sample_name_definition[1] == 'PacBio':
            probe_type = sample_name_definition[1]
        else:
            probe_type = (sample_name_definition[1] + ' ' 
                          + sample_name_definition[2] + ' ' 
                          + sample_name_definition[3])

        # If we are analyzing AMRs
        if amr_analysis:

            # This is the first time we see this probe for the organism for amr
            if probe_type not in list(self.abundance_dict_amr[organism].keys()):
                self.abundance_dict_amr[organism].update({probe_type:IndivAbundance(True)})

            self.abundance_dict_amr[organism][probe_type].add_to_absolute(
                tels_file_path(self, sample_name, self.amr_reads_ext))

        # If we are analyzing MGEs
        if mge_analysis:

            # This is the first time we see this probe for the organism for mge
            if probe_type not in list(self.abundance_dict_mge[organism].keys()):
                self.abundance_dict_mge[organism].update({probe_type:IndivAbundance(
                    False, self.mges_annot)})

            self.abundance_dict_mge[organism][probe_type].add_to_absolute(
                tels_file_path(self, sample_name, self.mge_reads_ext))
            
        # Overall, tThis is the first time we see this probe for the organism
        if probe_type not in list(self.initial_source_size[organism].keys()):
            self.initial_source_size[organism].update({probe_type:0}) 

        self.initial_source_size[organism][probe_type] += (
            float(self.all_file_sizes[sample_name]) / (10.0**9))

    def make_violin_plot(
            self, output_folder, violin_ext, amr_analysis, mge_analysis):

        def violin_plot_per_analysis(abundance_dict, genes_length, element_name):

            # Go through absolute abundace information to make it relative
            for organism in abundance_dict:
                for probe_type in abundance_dict[organism]:
                    abundance_dict[organism][probe_type].make_abundance_relative(
                        self.initial_source_size[organism][probe_type], genes_length)

            # Set up main figure        
            fig, axs = pyplot.subplots(
                nrows=1,
                ncols=4,
                sharey="row",
                figsize=(80, 20),
                layout="constrained")
            fig.suptitle(element_name + " Relative Abundance\n", fontsize=80)
            fig.add_subplot(111, frame_on=False)
            pyplot.grid(False)
            pyplot.tick_params(labelcolor="none", bottom=False, left=False)
            for ax in axs.flat:
                ax.set_ylabel ("Log Relative Abundance", fontsize=60)
                ax.label_outer()
                ax.set_ylim(-4,6)

            # Set up individual violin plot
            for index, organism in enumerate(abundance_dict):
                pyplot.sca(axs[index])
                x_axis_list = list()
                sample_abundance_dict = dict()
                for probe_type in abundance_dict[organism]:
                    sample_abundance_dict[probe_type] = (
                        abundance_dict[organism][probe_type].get_abundance())
                    x_axis_list.append(probe_type)
                df = pandas.DataFrame.from_dict(sample_abundance_dict)
                seaborn.set_style("whitegrid")
                seaborn.set_context("paper")

                axs[index] = seaborn.violinplot(
                        data=df, inner="box", palette="Greys")
                axs[index].grid(False)
                axs[index].set_title(organism, fontsize=70)

                if (index == 3) and (element_name == 'ARG'):
                    artists = axs[index].get_children()
                    legend_artists = [artists[0], artists[4], artists[8],
                               artists[12], artists[16], artists[20],
                               artists[24]]
                    pyplot.legend(
                        handles=legend_artists,
                        labels=x_axis_list,
                        loc='upper right',
                        fontsize=40,
                        ncols=2)

                axs[index].set_xticklabels(labels=[])
                
                if index == 0:
                    pyplot.yticks(fontsize=50)

            if not(os.path.exists(output_folder)):
                os.makedirs(output_folder)
            pyplot.gcf()
            pyplot.savefig(output_folder + element_name + violin_ext)
            pyplot.close()

        if amr_analysis:
            violin_plot_per_analysis(
                self.abundance_dict_amr, self.megares_genes_length, 'ARG')
        if mge_analysis:
            violin_plot_per_analysis(
                self.abundance_dict_mge, self.mge_genes_length, 'MGE')