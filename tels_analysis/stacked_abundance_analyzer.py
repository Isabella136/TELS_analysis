from tels_analysis.stacked_abundance_analysis.indiv_stacked_abundance import IndivStackedAbundance
from tels_analysis import get_sample_name_definition
from tels_analysis import get_mge_annot_dict
from tels_analysis import tels_file_path
from matplotlib import pyplot
import seaborn
import numpy
import os

class StackedAbundanceAnalyzer:
    
    def __init__(
            self, SOURCE_PREFIX, SOURCE_SUFFIX, AMR_DIV_EXT, 
            MGE_EXT, STATS_EXT, MGES_ANNOTATION):
        
        # Get tels output file names information
        self.source_prefix = SOURCE_PREFIX
        self.source_suffix = SOURCE_SUFFIX
        self.amr_reads_ext = AMR_DIV_EXT
        self.mge_reads_ext = MGE_EXT
        self.stats_ext = STATS_EXT

        # Get MGE types
        self.mges_annot = get_mge_annot_dict(MGES_ANNOTATION)

        # Estabilsh dictionaries for absolute abundance 
        self.abundance_dict_amr = dict()
        self.abundance_dict_mge = dict()

    def find_absolute_abundance(self, sample_name, amr_analysis, mge_analysis):

        amr_filepath = tels_file_path(self, sample_name, self.amr_reads_ext)
        mge_filepath = tels_file_path(self, sample_name, self.mge_reads_ext)
        stats_filepath = tels_file_path(self, sample_name, self.stats_ext)

        # Retrieve definition in tuple form Organism, Platform, Chemistry, Probe)
        # which must be used to declare the subtables (Organism + Chemistry) 
        # and the different groups in the legend (Platform + Probe)
        sample_name_definition = get_sample_name_definition(sample_name)
        subtable = sample_name_definition[0] + ' ' + sample_name_definition[2]
        legend = (sample_name_definition[1] if sample_name_definition[1] == 'PacBio'
                  else sample_name_definition[1] + ' ' + sample_name_definition[3])

        # If we are analyzing AMR
        if amr_analysis:

            # This is the first time we see this 
            # organism + chemistry combination for amr
            if subtable not in self.abundance_dict_amr:
                self.abundance_dict_amr[subtable] = IndivStackedAbundance(True)

            self.abundance_dict_amr[subtable].add_to_absolute(
                legend, amr_filepath, stats_filepath)

        # If we are analyzing MGEs
        if mge_analysis: 

        # This is the first time we see this 
        # organism + chemistry combination for mge
            if subtable not in self.abundance_dict_mge:
                self.abundance_dict_mge[subtable] = IndivStackedAbundance(
                    False, self.mges_annot)
                
            self.abundance_dict_mge[subtable].add_to_absolute(
                legend, mge_filepath, stats_filepath)
                
    def make_stacked_barplot(
            self, output_folder, stacked_ext, amr_analysis, mge_analysis):
        
        def superplot_per_analysis(abundance_dict, element_name):

            # Go through absolute abundace information to make it relative
            for subtable in abundance_dict:
                abundance_dict[subtable].make_abundance_relative()

            # Set up main figure
            fig, axs = pyplot.subplots(
                nrows=3,
                ncols=8,
                figsize=(60, 20),
                gridspec_kw={'height_ratios': [2,.1,1.5]})
            fig.suptitle(
                'Relative Abundance & ' + element_name + ' Richness', fontsize=50)

            # Go through each subtable
            for index, subtable in enumerate(abundance_dict):

                # Stacked bar plot
                pyplot.sca(axs[0][index])
                color_list = seaborn.color_palette("colorblind", n_colors=4)
                sub_abundance = abundance_dict[subtable].get_abundance()
                bottom_array = None
                for l_index, legend in enumerate(sub_abundance):
                    current_array = numpy.array(list(sub_abundance[legend].values()))
                    x_coords = list(sub_abundance[legend].keys())
                    if bottom_array is None:
                        axs[0][index].bar(x_coords, current_array, width=1.0,label=legend,
                                          color=color_list[l_index], alpha=0.5)
                        bottom_array = current_array
                    else:
                        axs[0][index].bar(x_coords, current_array, width=1.0, label=legend,
                                          color=color_list[l_index], alpha=0.5, bottom=bottom_array)
                        bottom_array = bottom_array + current_array
                if index == 0:
                    axs[0][index].set_ylabel('Log Relative Abundance', size=25)
                    pyplot.yticks(fontsize=20)
                else:
                    axs[0][index].sharey(axs[0][0])
                pyplot.xticks([])
                pyplot.margins(x=0)
                axs[0][index].set_title(subtable, size=30)
                axs[0][index].set_anchor('NE')
                if index == 7: axs[0][index].legend()   # only show legend color at the last subplot

                # Color-coded x-axis: 
                #   assign number to each gene based on category alphabetical sorting
                #   and use this number to determine color shown on x_axis
                pyplot.sca(axs[1][index])
                gene_to_category, category_list = abundance_dict[subtable].get_categories()
                x_matrix = list()
                for arg in gene_to_category:
                    x_matrix.append(category_list.index(gene_to_category[arg]))
                numpy_array = numpy.array([x_matrix])
                seaborn.heatmap(numpy_array, ax = axs[1][index], xticklabels=False, 
                                yticklabels=False, cbar=False, cmap='viridis')
                axs[1][index].set_anchor('NE')

                # Legend for color-coded x-axis
                pyplot.sca(axs[2][index])
                label_matrix = [*range(len(category_list))]
                numpy_array = numpy.array(label_matrix).reshape(len(label_matrix),1)
                seaborn.heatmap(numpy_array, ax = axs[2][index], square=True, 
                                xticklabels=False, cbar=False, cmap='viridis', 
                                linewidths=1)
                pyplot.yticks(
                    ticks=pyplot.yticks()[0], labels=category_list, fontsize=20, rotation=0)
                axs[2][index].set_anchor('NE')

            if not(os.path.exists(output_folder)):
                os.makedirs(output_folder)
            pyplot.gcf()
            pyplot.savefig(output_folder + element_name + stacked_ext)
            pyplot.close()

        if amr_analysis:
            superplot_per_analysis(self.abundance_dict_amr, 'ARG')
        if mge_analysis:
            superplot_per_analysis(self.abundance_dict_mge, 'MGE')
