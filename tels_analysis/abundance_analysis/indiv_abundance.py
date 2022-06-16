from math import log10
import numpy

class indiv_abundance:
    def __init__(this, sample_type, x_axis):
        this.x_axis = x_axis
        this.sample_type = sample_type
        this.absolute_abundance = {}
        this.relative_abundance = {}

    def addToAbsolute(this, filepath):
        argFile = open(filepath, "r")
        for line in argFile:
            blank_space = 0
            arg_header_list = line.split(',')[1:]
            for arg in arg_header_list:
                if arg == "":
                    blank_space += 1
                    if (blank_space % 2) == 0:
                        break
                    continue
                arg_header = arg.split("|")
                if arg_header[0] not in this.absolute_abundance:
                    this.absolute_abundance[arg_header[0]] = 0
                this.absolute_abundance[arg_header[0]] += 1
