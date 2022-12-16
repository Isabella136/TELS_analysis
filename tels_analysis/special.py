from re import M
from tels_analysis import mgeDict

class special:
    filename = lambda this, fileName, : this.source_prefix + fileName + this.source_suffix + this.extension

    def __init__(this, SOURCE_PREFIX, MGEAlignedToMegares, SOURCE_SUFFIX, SHORT_MGE, MGE_CLASSIFICATION):
        this.source_prefix = SOURCE_PREFIX
        this.mge_aligned_to_megares = []
        file = open(MGEAlignedToMegares, "r")
        for line in file:
            mge1, mge2 = tuple(line.split(","))
            mge2 = mge2[:-1]
            this.mge_aligned_to_megares.append(mge1)
            if mge2 != "":
                this.mge_aligned_to_megares.append(mge2)
        file.close()
        this.source_suffix = SOURCE_SUFFIX
        this.extension = SHORT_MGE
        this.mge_dict = mgeDict(MGE_CLASSIFICATION)
        this.sample_list = []
        this.classified_arg = {}
        this.classified_aligned_arg = {}
        this.aligned_arg = {}
        this.richness = {}

    def writeMobilomeInfo(this, output):
        file = open(output, "w")
        for sample in this.sample_list:
            file.write("," + sample)
        file.write("\n")
        for mge, info in this.classified_arg.items():
            file.write(mge)
            for sample in this.sample_list:
                if sample not in info:
                    file.write(",0")
                else:
                    file.write("," + str(info[sample]))
            file.write("\n")
        file.write("RICHNESS:")
        for sample in this.richness:
            file.write("," + str(this.richness[sample]["classified_arg"]))
        file.write("\n")
        file.write("\n")
        for mge, info in this.classified_aligned_arg.items():
            file.write(mge)
            for sample in this.sample_list:
                if sample not in info:
                    file.write(",0")
                else:
                    file.write("," + str(info[sample]))
            file.write("\n")
        file.write("RICHNESS:")
        for sample in this.richness:
            file.write("," + str(this.richness[sample]["classified_aligned_arg"]))
        file.write("\n")
        file.write("\n")
        for mge, info in this.aligned_arg.items():
            file.write(mge)
            for sample in this.sample_list:
                if sample not in info:
                    file.write(",0")
                else:
                    file.write("," + str(info[sample]))
            file.write("\n")
        file.write("RICHNESS:")
        for sample in this.richness:
            file.write("," + str(this.richness[sample]["aligned_arg"]))
        file.close()

    def addToMobilomeInfo(this, sample):
        this.richness.update({sample: {"classified_arg": 0, "classified_aligned_arg": 0, "aligned_arg": 0}})
        this.sample_list.append(sample)
        file = open(this.filename(sample), "r")
        i = 0
        for line in file:
            i += 1
            if i < 21:
                continue
            mge = line.split(',')[0]
            count = line.split(',')[1][:-1]
            annot = this.mge_dict[mge]
            if annot == "AMR":
                if mge not in this.mge_aligned_to_megares:
                    if mge not in this.classified_arg:
                        this.classified_arg.update({mge: dict()})
                    this.classified_arg[mge].update({sample: count})
                    this.richness[sample]["classified_arg"] += 1
                else:
                    if mge not in this.classified_aligned_arg:
                        this.classified_arg.update({mge: dict()})
                        this.classified_aligned_arg.update({mge: dict()})
                        this.aligned_arg.update({mge: dict()})
                    this.classified_arg[mge].update({sample: count})
                    this.classified_aligned_arg[mge].update({sample: count})
                    this.aligned_arg[mge].update({sample: count})
                    this.richness[sample]["classified_arg"] += 1
                    this.richness[sample]["classified_aligned_arg"] += 1
                    this.richness[sample]["aligned_arg"] += 1
            else:
                if mge in this.mge_aligned_to_megares:
                    if mge not in this.aligned_arg:
                        this.aligned_arg.update({mge: dict()})
                    this.aligned_arg[mge].update({sample: count})
                    this.richness[sample]["aligned_arg"] += 1
        file.close()
        