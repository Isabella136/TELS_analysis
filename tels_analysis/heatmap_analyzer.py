from tels_analysis.heatmap_analysis.heatmap_creation.indiv_heatmap import IndivHeatmap
from tels_analysis.heatmap_analysis.double_heatmap_creator import DoubleHeatmapCreator
from tels_analysis import get_sample_name_definition
from tels_analysis import megares_analyzer
from tels_analysis import tels_file_path
from tels_analysis import mge_analyzer
from matplotlib import pyplot
import seaborn
import numpy
import os

class HeatmapAnalyzer:
    def __init__(
            self, SOURCE_PREFIX, SOURCE_SUFFIX, AMR_DIV_EXT,
            MGE_EXT, MEGARES_ANNOTATION, MGES_ANNOTATION,
            AMR_ANALYSIS, MGE_ANALYSIS):

        # Get tels output file names information
        self.source_prefix = SOURCE_PREFIX
        self.source_suffix = SOURCE_SUFFIX
        self.amr_reads = AMR_DIV_EXT
        self.mge_reads = MGE_EXT

        # Functions required to initialize object variables
        def paired_list_to_dict(paired_list):
            left_dict = dict()  #keeps track of class and type count
            right_dict = dict() #mech is in class ____; accession is ___ mge type
            for tuple in paired_list:
                right_dict.update({tuple[1]:tuple[0]})
                if left_dict.get(tuple[0], 0) == 0:
                    left_dict.update({tuple[0]:1})
                else:
                    left_dict[tuple[0]] += 1
            return (left_dict, right_dict)

        def paired_list_to_bool(paired_list):
            bool_dict = dict()
            for tuple in paired_list:
                bool_dict.update({tuple:False})
            return bool_dict

        # Set up variables for AMR analysis
        self.amr_analysis = AMR_ANALYSIS
        if self.amr_analysis:
            drug_list, other_list = megares_analyzer(MEGARES_ANNOTATION)
            self.drug_class_dict, drug_mech_dict = paired_list_to_dict(drug_list)
            self.other_class_dict, other_mech_dict = paired_list_to_dict(other_list)
            self.drug_bool_dict = paired_list_to_bool(drug_list)
            self.other_bool_dict = paired_list_to_bool(other_list)
            self.megares_heatmap_list = [
                DoubleHeatmapCreator("BF", drug_mech_dict, other_mech_dict),
                DoubleHeatmapCreator("FMT", drug_mech_dict, other_mech_dict),
                DoubleHeatmapCreator("PPS", drug_mech_dict, other_mech_dict),
                DoubleHeatmapCreator("MOCK", drug_mech_dict, other_mech_dict)]

        # Set up variables for MGE analysis
        self.mge_analysis = MGE_ANALYSIS
        if self.mge_analysis:
            mge_list = mge_analyzer(MGES_ANNOTATION)
            self.mge_annotations = {t[1]:t[0] for t in mge_list}
            self.mge_type_dict, mge_access_dict = paired_list_to_dict(mge_list)
            self.mge_bool_dict = {
                'BF': paired_list_to_bool(mge_list),
                'FMT': paired_list_to_bool(mge_list),
                'PPS': paired_list_to_bool(mge_list)}
            self.mge_heatmap_list = [
                IndivHeatmap("BF", mge_access_dict),
                IndivHeatmap("FMT", mge_access_dict),
                IndivHeatmap("PPS", mge_access_dict)]

    def add_to_maps(self, sample_name):
        sample_name_definition = get_sample_name_definition(sample_name, True)
        organism = sample_name_definition[0]
        index = ['BF', 'FMT', 'PPS', 'MOCK'].index(organism)
        if sample_name_definition[1] == 'PacBio':
            probe_type = sample_name_definition[1]
            if sample_name_definition[2] == 'V2':
                if sample_name[-1] == 'A':
                    duplicate = 'A'
                elif sample_name[-2:] == 'AM':
                    duplicate = 'B'
                else:
                    duplicate = 'C'
            else:
                if sample_name[-1] == 'A':
                    duplicate = 'D'
                elif sample_name[-2:] == 'AM':
                    duplicate = 'E'
                else:
                    duplicate = 'F'
        else:
            probe_type = (sample_name_definition[2] + ' ' + sample_name_definition[3])
            duplicate = sample_name[-1]

        if self.amr_analysis:
            self.megares_heatmap_list[index].add_to_maps(
                probe_type, tels_file_path(self, sample_name, self.amr_reads),
                duplicate, self.drug_bool_dict, self.other_bool_dict)
        if self.mge_analysis:
            if organism == 'MOCK' or sample_name == 'BFXTAMA':
                return None
            self.mge_heatmap_list[index].add_to_map(
                probe_type, tels_file_path(self, sample_name, self.mge_reads),
                duplicate, self.mge_bool_dict[organism], self.mge_annotations)

    def make_maps(self, output_folder, heatmap_ext):

        def heatmap_maker(
                matrix_list, x_axis_list, category_dict, file_prefix, element_name):
            
            # Calculate height of heatmap and position of labels
            label_matrix = list()
            label_pos = list()
            indiv_point_total = 0
            for category_index, category in enumerate(category_dict):
                for i in range(category_dict[category]):
                    label_matrix.append(category_index+1)
                down = indiv_point_total
                up = down + category_dict[category]
                if '\n' in category:
                    label_pos.append((up + down-.4)/2)
                else:
                    label_pos.append((up + down)/2)
                indiv_point_total += category_dict[category]

            # Create figure
            fig = pyplot.figure(figsize=(
                (65, 85) if 'ARG' in element_name
                else (25,85)))
            supersubfigs = fig.subfigures(
                nrows=2,
                ncols=1,
                height_ratios=[1,20]
            )
            supersubfigs[0].suptitle(element_name, fontsize=110, y=.75)

            if 'MGE' not in element_name:
                subfigs = supersubfigs[1].subfigures(
                    nrows=2,
                    ncols=7,
                    wspace=0,
                    hspace=0,
                    width_ratios=[12,2,20,20,20,20,1],
                    height_ratios=[1,80]
                )
            else:
                subfigs = supersubfigs[1].subfigures(
                    nrows=2,
                    #ncols=6,
                    ncols=4,
                    wspace=0,
                    hspace=0,
                    #width_ratios=[12,2,20,20,20,1],
                    width_ratios=[12,2,20,1],
                    height_ratios=[1,80]
                )

            # Create heatmap legend
            data = numpy.array(label_matrix).reshape(len(label_matrix),1)
            axs1 = subfigs[1][1].subplots(1,1)
            seaborn.heatmap(
                data,
                ax=axs1,
                cbar=False,
                cmap="viridis",
                xticklabels=False,
                vmin=(numpy.max(label_matrix) - 7 if 'MGE' in element_name else 0), 
                vmax=numpy.max(label_matrix))
            pyplot.sca(axs1)
            pyplot.yticks(
                label_pos, list(category_dict.keys()), fontsize=65, rotation=0)

            # Going through organsims
            if 'MGE' not in element_name:
                organism_list = ['BF', 'FMT', 'PPS', 'MOCK']
                for index, organism in enumerate(organism_list):
                    axs = subfigs[1][2+index].subfigures(
                        nrows=1,
                        ncols=2,
                        wspace=-0.1,
                        width_ratios=[6,1])
                    subfigs[0][index+2].suptitle(organism, fontsize=90, y=.95)
                    data = numpy.array(matrix_list[index])

                    sub_axs = axs[0].subplots(1,6)
                    axs[0].suptitle('TELSeq', fontsize=70)
                    pacbio_index_start = 18
                    for column in range(6):
                        seaborn.heatmap(
                            data[:,column*3:column*3+3],
                            ax=sub_axs[column],
                            cbar=False,
                            cmap="viridis",
                            yticklabels=False,
                            xticklabels=False,
                            vmin=0, 
                            vmax=numpy.max(label_matrix))
                        pyplot.sca(sub_axs[column])
                        pyplot.xlabel(x_axis_list[index][column], fontsize=65, rotation=90)

                    axs[0].subplots_adjust(wspace=0.05)

                    sub_axs1 = axs[1].subplots(1,1)
                    axs[1].suptitle('PacBio', fontsize=70, x=0.30)
                    seaborn.heatmap(
                        data[:,pacbio_index_start:],
                        ax=sub_axs1,
                        cbar=False,
                        cmap="viridis",
                        yticklabels=False,
                        xticklabels=False,
                        vmin=0, 
                        vmax=numpy.max(label_matrix))
            else:
                axs = subfigs[1][2].subfigures(
                    nrows=1,
                    ncols=2,
                    wspace=-0.1,
                    width_ratios=[6,1])
                data = numpy.array(matrix_list)

                sub_axs = axs[0].subplots(1,6)
                axs[0].suptitle('TELSeq', fontsize=70)
                pacbio_index_start = 18
                for column in range(6):
                    if ((column == 5)
                        and ('BF' in element_name)):    # Must remove outlier
                        pacbio_index_start = 17
                        seaborn.heatmap(
                            data[:,column*3:column*3+2],
                            ax=sub_axs[column],
                            cbar=False,
                            cmap="viridis",
                            yticklabels=False,
                            xticklabels=False,
                            vmin=numpy.max(label_matrix) - 7, 
                            vmax=numpy.max(label_matrix))
                    else:
                        seaborn.heatmap(
                            data[:,column*3:column*3+3],
                            ax=sub_axs[column],
                            cbar=False,
                            cmap="viridis",
                            yticklabels=False,
                            xticklabels=False,
                            vmin=numpy.max(label_matrix) - 7, 
                            vmax=numpy.max(label_matrix))
                    pyplot.sca(sub_axs[column])
                    pyplot.xlabel(x_axis_list[column], fontsize=65, rotation=90)

                axs[0].subplots_adjust(wspace=0.05)

                sub_axs1 = axs[1].subplots(1,1)
                axs[1].suptitle('PacBio', fontsize=70, x=0.30)
                seaborn.heatmap(
                    data[:,pacbio_index_start:],
                    ax=sub_axs1,
                    cbar=False,
                    cmap="viridis",
                    yticklabels=False,
                    xticklabels=False,
                    vmin=numpy.max(label_matrix) - 7, 
                    vmax=numpy.max(label_matrix))
            

            if not(os.path.exists(output_folder)):
                os.makedirs(output_folder)

            pyplot.gcf().subplots_adjust(top=0.95, bottom=0.1)
            pyplot.savefig(output_folder + file_prefix + heatmap_ext)
            pyplot.close()

        if self.amr_analysis:
            for tuple in self.drug_bool_dict:
                if not(self.drug_bool_dict[tuple]):
                    self.drug_class_dict[tuple[0]] -= 1
                    if self.drug_class_dict[tuple[0]] == 0:
                        self.drug_class_dict.pop(tuple[0])
            for tuple in self.other_bool_dict:
                if not(self.other_bool_dict[tuple]):
                    self.other_class_dict[tuple[0]] -= 1
                    if self.other_class_dict[tuple[0]] == 0:
                        self.other_class_dict.pop(tuple[0])
            
            drug_matrix_list = list()
            drug_x_axis_list = list()
            other_matrix_list = list()
            other_x_axis_list = list()
            for heatmap in self.megares_heatmap_list:
                # output_tuple = (drug_matrix, drug_x_axis, other_matrix, other_x_axis)
                output_tuple = heatmap.make_maps(
                    self.drug_bool_dict, self.other_bool_dict,
                    list(self.drug_class_dict.keys()),
                    list(self.other_class_dict.keys()))
                drug_matrix_list.append(output_tuple[0])
                drug_x_axis_list.append(output_tuple[1])
                other_matrix_list.append(output_tuple[2])
                other_x_axis_list.append(output_tuple[3])

            heatmap_maker(drug_matrix_list, drug_x_axis_list,
                        self.drug_class_dict, 'ARG_drugs', 'ARG - Drugs')
            heatmap_maker(other_matrix_list, other_x_axis_list,
                        self.other_class_dict, 'ARG_bio_metals', 'ARG - Biocides/Metals')
            
        if self.mge_analysis:
            for org_index, organism in enumerate(self.mge_bool_dict):
                temp_mge_type_dict = {k:v for k,v in self.mge_type_dict.items()}
                for tuple in self.mge_bool_dict[organism]:
                    if not(self.mge_bool_dict[organism][tuple]):
                        temp_mge_type_dict[tuple[0]] -= 1
                        if temp_mge_type_dict[tuple[0]] == 0:
                            temp_mge_type_dict.pop(tuple[0])

                # output_tuple = (mge_matrix, mge_x_axis)
                output_tuple = self.mge_heatmap_list[org_index].make_map(
                    self.mge_bool_dict[organism], list(temp_mge_type_dict.keys()))
                heatmap_maker(output_tuple[0], output_tuple[1],
                                temp_mge_type_dict, 'MGE_' + organism, 'MGE - ' + organism)

            

            