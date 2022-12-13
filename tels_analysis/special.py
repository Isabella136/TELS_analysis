from re import M
from tels_analysis import mgeDict

class special:
    fiftyFile = lambda this, fileName, : this.source_prefix + fileName + this.source_suffix + this.extension
    eightyFile = lambda this, fileName : this.special_prefix + fileName + this.source_suffix + this.extension

    def __init__(this, SOURCE_PREFIX, SPECIAL_PREFIX, SOURCE_SUFFIX, SHORT_MGE, MGE_CLASSIFICATION):
        this.source_prefix = SOURCE_PREFIX
        this.special_prefix = SPECIAL_PREFIX
        this.source_suffix = SOURCE_SUFFIX
        this.extension = SHORT_MGE
        this.mge_dict = mgeDict(MGE_CLASSIFICATION)
        this.sample_list = []
        this.amr_mobilome = {}
        this.unknown_mobilome = {}

    def writeMobilomeInfo(this, output):
        file = open(output, "w")
        for sample in this.sample_list:
            file.write("," + sample + ",")
        file.write("\n")
        for i in range(len(this.sample_list)):
            file.write(",0.5,0.8")
        file.write("\n")
        for mge, info in this.amr_mobilome.items():
            file.write(mge)
            for sample in this.sample_list:
                if sample not in info:
                    file.write(",0,0")
                else:
                    file.write("," + str(info[sample][50]) + "," + str(info[sample][80]))
            file.write("\n")
        file.write("\n")
        for mge, info in this.unknown_mobilome.items():
            file.write(mge)
            for sample in this.sample_list:
                if sample not in info:
                    file.write(",0,0")
                else:
                    file.write("," + str(info[sample][50]) + "," + str(info[sample][80]))
            file.write("\n")
        file.close()

    def addToMobilomeInfo(this, sample):
        this.sample_list.append(sample)
        fifty = open(this.fiftyFile(sample), "r")
        i = 0
        for line in fifty:
            i += 1
            if i < 21:
                continue
            mge = line.split(',')[0]
            count = line.split(',')[1][:-1]
            if mge not in this.mge_dict:
                if mge not in this.unknown_mobilome:
                    this.unknown_mobilome.update({mge: dict()})
                this.unknown_mobilome[mge].update({sample: {50: count, 80: 0}})
            else:
                annot = this.mge_dict[mge]
                if annot == "UNKNOWN":
                    if mge not in this.unknown_mobilome:
                        this.unknown_mobilome.update({mge: dict()})
                    this.unknown_mobilome[mge].update({sample: {50: count, 80: 0}})
                elif annot == "AMR":
                    if mge not in this.amr_mobilome:
                        this.amr_mobilome.update({mge: dict()})
                    this.amr_mobilome[mge].update({sample: {50: count, 80: 0}})
        fifty.close()
        eighty = open(this.eightyFile(sample), "r")
        i = 0
        for line in eighty:
            i += 1
            if i < 21:
                continue
            mge = line.split(',')[0]
            count = line.split(',')[1][:-1]
            if mge not in this.mge_dict:
                if mge not in this.unknown_mobilome:
                    this.unknown_mobilome.update({mge: dict()})
                if sample not in this.unknown_mobilome[mge]:
                    this.unknown_mobilome[mge].update({sample: {50: 0, 80: 0}})
                this.unknown_mobilome[mge][sample][80] = count
            else:
                annot = this.mge_dict[mge]
                if annot == "UNKNOWN":
                    if mge not in this.unknown_mobilome:
                        this.unknown_mobilome.update({mge: dict()})
                    if sample not in this.unknown_mobilome[mge]:
                        this.unknown_mobilome[mge].update({sample: {50: 0, 80: 0}})
                    this.unknown_mobilome[mge][sample][80] = count
                elif annot == "AMR":
                    if mge not in this.amr_mobilome:
                        this.amr_mobilome.update({mge: dict()})
                    if sample not in this.amr_mobilome[mge]:
                        this.amr_mobilome[mge].update({sample: {50: 0, 80: 0}})
                    this.amr_mobilome[mge][sample][80] = count
        