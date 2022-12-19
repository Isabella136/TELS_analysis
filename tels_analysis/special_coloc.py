class special_coloc:
    filename = lambda this, fileName: this.source_prefix + fileName + this.source_suffix + this.extension

    def __init__(this, SOURCE_PREFIX, MGEAlignedToMegares, SOURCE_SUFFIX, COMPOS_RICHNESS_ANALYSIS):
        this.source_prefix = SOURCE_PREFIX
        this.mge_classified_aligned_to_arg = []
        file = open(MGEAlignedToMegares, "r")
        for line in file:
            mge1, mge2 = tuple(line.split(","))
            this.mge_classified_aligned_to_arg.append(mge1)
        file.close()
        this.source_suffix = SOURCE_SUFFIX
        this.extension = COMPOS_RICHNESS_ANALYSIS
        this.sample_list = []
        this.coloc_info = {}

    def addColocInfo(this, sample):
        this.coloc_info.update({sample : 0})
        this.sample_list.append(sample)
        special_mge_count = 0