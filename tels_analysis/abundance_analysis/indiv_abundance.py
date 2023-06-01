from math import log10
import statistics

class indiv_abundance:
    def __init__(this, sample_type, x_axis):
        this.x_axis = x_axis
        this.sample_type = sample_type
        this.accession_to_group = dict()
        this.absolute_abundance = {}
        this.relative_abundance = {}

    def addToAbsolute(this, filepath):
        argFile = open(filepath, "r")
        lineNum = 0
        for line in argFile:
            lineNum += 1
            if lineNum < 20:
                continue
            line_list = line.split(',')
            arg_header = line_list[0].split("|")
            if arg_header[0] not in this.absolute_abundance:
                this.absolute_abundance[arg_header[0]] = 0
                this.accession_to_group[arg_header[0]] = arg_header[4]
            this.absolute_abundance[arg_header[0]] += int(line_list[1])
    def makeAbundanceRelative(this, file_size, gene_length):
        def aggregated_lengths():
            all_lengths_per_group = dict()
            for arg, length in gene_length.items():
                if arg not in this.accession_to_group: continue
                if this.accession_to_group[arg] not in all_lengths_per_group:
                    all_lengths_per_group[this.accession_to_group[arg]] = list()
                all_lengths_per_group[this.accession_to_group[arg]].append(length)
            return all_lengths_per_group
        def median_gene_length():
            all_lengths_per_group = aggregated_lengths()
            median_gene_length = dict()
            for group, length_list in all_lengths_per_group.items():
                median_gene_length[group] = statistics.median(length_list)
            return median_gene_length
        def average_gene_length():
            all_lengths_per_group = aggregated_lengths()
            average_gene_length = dict()
            for group, length_list in all_lengths_per_group.items():
                average_gene_length[group] = sum(length_list)/len(length_list)
            return average_gene_length
        
        #new_gene_length = median_gene_length()
        #new_gene_length = average_gene_length()

        not_logged_relative_abundance = dict()
        for arg, count in this.absolute_abundance.items():
            if this.accession_to_group[arg] not in not_logged_relative_abundance:
                not_logged_relative_abundance[this.accession_to_group[arg]] = 0

            # Required for summation
            not_logged_relative_abundance[this.accession_to_group[arg]] += (
                (100.0 * float(count))/(float(gene_length[arg]*file_size))
            )

            # Required for median or average
            # not_logged_relative_abundance[this.accession_to_group[arg]] += (
            #     (100.0 * float(count))/(float(new_gene_length[this.accession_to_group[arg]]*file_size))
            # )

            #this.relative_abundance[arg] = log10((100.0 * float(count))/(float(gene_length[arg]*file_size)))
        for group, relative_abundance in not_logged_relative_abundance.items():
            this.relative_abundance[group] = log10(relative_abundance)
    def getAbundance(this):
        return this.relative_abundance