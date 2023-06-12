from tels_analysis.compos_richness_analysis.indiv_composition import indiv_composition
from tels_analysis import get_mge_annot_dict
from matplotlib import pyplot
from PIL import Image
import seaborn
import numpy
import os


class compos_richness_analyzer:
    filePath = lambda this, fileName, extension : this.source_prefix + fileName + this.source_suffix + extension

    def __init__(this, SOURCE_PREFIX, SOURCE_SUFFIX, SHORT_AMR_DIV, SHORT_MGE, MGE_CLASSIFICATION):
        this.source_prefix = SOURCE_PREFIX
        this.source_suffix = SOURCE_SUFFIX
        this.amr_reads = SHORT_AMR_DIV
        this.mge_reads = SHORT_MGE
        this.mge_dict = get_mge_annot_dict(MGE_CLASSIFICATION)
        this.compos_richness_list = []
        this.composition_data = []
        this.MGE_list = []

    def analyzeFile(this, fileName):
        fileComp = indiv_composition(fileName, this.mge_dict)
        
        fileComp.findAllData(this.filePath(fileName, this.amr_reads), this.filePath(fileName, this.mge_reads))
        this.composition_data.append(fileComp.getComposition())
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

    def makeBarCharts(this, argOutputFolder, mgeOutputFolder, INDIV_COMPOS_CHART):
        prime_numbers = [2, 3]
        for num in range(1, 19134):
            is_prime = True
            for prime_num in prime_numbers:
                if (6*num-1) % prime_num == 0:
                    is_prime = False
                    break
                if (prime_num*prime_num) > (6*num-1):
                    break
            if is_prime:
                prime_numbers.append((6*num-1))

            is_prime = True
            for prime_num in prime_numbers:
                if (6*num+1) % prime_num == 0:
                    is_prime = False
                    break
                if (prime_num*prime_num) > (6*num+1):
                    break
            if is_prime:
                prime_numbers.append((6*num+1))
        def factorization(num):
            factors = []
            remainder = num
            prime_numbers_index = 0
            while (remainder > prime_numbers[prime_numbers_index]):
                if remainder % prime_numbers[prime_numbers_index] == 0:
                    factors.append(prime_numbers[prime_numbers_index])
                    remainder = int(remainder / prime_numbers[prime_numbers_index])
                else:
                    prime_numbers_index += 1
            factors.append(remainder)
            return factors
        def least_common_multiple(sum_list):
            all_nums = []
            for value in sum_list:
                if value != 0 and value not in all_nums:
                    all_nums.append(value)
            if len(all_nums) == 1:
                return all_nums[0]
            prev_LCM = all_nums[0]
            for num in all_nums[1:]:
                prev_factorization = factorization(prev_LCM)
                curr_factorization = factorization(num)
                greatest_common_denominator = 1
                next_curr_index = 0
                for prev_factor in prev_factorization:
                    while (next_curr_index < len(curr_factorization) 
                           and curr_factorization[next_curr_index] < prev_factor):
                        next_curr_index += 1
                    if next_curr_index == len(curr_factorization):
                        break
                    if prev_factor == curr_factorization[next_curr_index]:
                        greatest_common_denominator *= prev_factor
                        next_curr_index += 1
                prev_LCM = prev_LCM * int(num / greatest_common_denominator)
            return prev_LCM
            
        tab10_colorblind = ['#006BA4', '#FF800E', '#ABABAB', '#595959', 
                            '#5F9ED1', '#C85200', '#898989', '#A2C8EC', 
                            '#FFBC79', '#CFCFCF']

        if not(os.path.exists(argOutputFolder)):
            os.makedirs(argOutputFolder)
        if not(os.path.exists(mgeOutputFolder)):
            os.makedirs(mgeOutputFolder)
        for index in range(0, len(this.composition_data), 3):
            arg_composition_data_sum = [sum(this.composition_data[index][1]),
                                        sum(this.composition_data[index+1][1]),
                                        sum(this.composition_data[index+2][1])]
            mge_composition_data_sum = [sum(this.composition_data[index][2]),
                                        sum(this.composition_data[index+1][2]),
                                        sum(this.composition_data[index+2][2])]

            filepath_arg_output = argOutputFolder + "/" + this.composition_data[index][0] + INDIV_COMPOS_CHART
            filepath_mge_output = mgeOutputFolder + "/" + this.composition_data[index][0] + INDIV_COMPOS_CHART

            if sum(arg_composition_data_sum) != 0:
                height = least_common_multiple(arg_composition_data_sum)
            else:
                height = 1
            new_list = list()
            for class_num in range(4):
                row = list()
                for sample_num in range(3):
                    if arg_composition_data_sum[sample_num] == 0:
                        row.append(0)
                    else:
                        proportion = int(height / arg_composition_data_sum[sample_num])
                        row.append(this.composition_data[index+sample_num][1][class_num]*proportion)
                new_list.append(row)
            
            bottom_array = [0,0,0]
            pyplot.figure(figsize=[16,30])
            for row, c_index in zip(new_list, range(4)):
                curr_array = numpy.array(row)
                pyplot.bar(['A','B','C'], 
                            width=0.9,
                            height=curr_array,
                            bottom=bottom_array,
                            #color=tab10_colorblind[c_index],
                            color=seaborn.color_palette("colorblind", n_colors=4)[c_index],
                            label=this.composition_data[index][0] + " ARG")
                bottom_array = curr_array + bottom_array
            pyplot.legend(['Drugs', 'Metal', 'Biocide', 'Multi-Compound'], 
                          bbox_to_anchor=(0.5, -0.03),
                          loc='lower center',
                          reverse=True,
                          fontsize=20,
                          ncol=4)
            pyplot.title(this.composition_data[index][0] + " ARG", fontsize=20)
            pyplot.ylim(0,height)
            pyplot.xticks([])
            pyplot.yticks([])
            pyplot.tight_layout()
            pyplot.savefig(filepath_arg_output)
            pyplot.close()

            if sum(mge_composition_data_sum) != 0:
                height = least_common_multiple(mge_composition_data_sum)
            else:
                height = 1
            new_list = list()
            for class_num in range(7):
                row = list()
                for sample_num in range(3):
                    if mge_composition_data_sum[sample_num] == 0:
                        row.append(0)
                    else:
                        proportion = int(height / mge_composition_data_sum[sample_num])
                        row.append(this.composition_data[index+sample_num][2][class_num]*proportion)
                new_list.append(row)
            
            bottom_array = [0,0,0]
            pyplot.figure(figsize=[16,30])
            for row, c_index in zip(new_list, range(7)):
                curr_array = numpy.array(row)
                pyplot.bar(['A','B','C'], 
                            width=0.9,
                            height=curr_array, 
                            bottom=bottom_array,
                            #color=tab10_colorblind[c_index],
                            color=seaborn.color_palette("colorblind", n_colors=13)[12-c_index],
                            label=this.composition_data[index][0] + " MGE")
                bottom_array = curr_array + bottom_array
            # pyplot.legend(['Plasmid', 'Phage', 'TE', 'IS', 'ICE', 'Virus', 'Unclassified'], 
            #               bbox_to_anchor=(0.5, -0.03),
            #               loc='lower center',
            #               reverse=True,
            #               fontsize=20,
            #               ncol=7)
            # pyplot.title(this.composition_data[index][0] + " MGE", fontsize=20)
            pyplot.ylim(0,height)
            pyplot.xticks([])
            pyplot.yticks([])
            pyplot.tight_layout()
            pyplot.savefig(filepath_mge_output)
            pyplot.close()

    def printAnalysis(this, outputFolder, COMPOS_RICHNESS_ANALYSIS):
        if not(os.path.exists(outputFolder)):
                os.makedirs(outputFolder)
        analysis = open(outputFolder + "/" + COMPOS_RICHNESS_ANALYSIS, "w")
        analysis.write("Sample,Sequencing platform,ARG composition,ARG Class richness,ARG Mechanism richness,ARG Group richness,MGE composition,MGE Accession richness\n")
        for c_r in this.compos_richness_list:
            analysis.write(c_r + "\n")
        analysis.close()

        analysis = open(outputFolder + "/MGE_list.csv", "w")
        for MGE in this.MGE_list:
            analysis.write(MGE + "\n")
        analysis.close()
