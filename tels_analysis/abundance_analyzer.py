from tels_analysis.abundance_analysis.indiv_abundance import indiv_abundance
from tels_analysis import x_axis
from tels_analysis import getGenesLength
from matplotlib import pyplot
import seaborn, numpy, gzip, shutil, os

class abundance_analyzer:
    filePath = lambda this, fileName, extension, prefix : prefix + fileName + this.source_suffix + extension

    def __init__(this, ARG_SAM_ANALYSIS_SOURCE_PREFIX, SOURCE_SUFFIX, INITIAL_SOURCE_PREFIX, ARG_SAM_ANALYSIS, MGE_SAM_ANALYSIS, MEGARES):
        this.source_prefix = ARG_SAM_ANALYSIS_SOURCE_PREFIX
        this.source_suffix = SOURCE_SUFFIX
        this.initial_source_prefix = INITIAL_SOURCE_PREFIX
        this.amr_reads = ARG_SAM_ANALYSIS
        this.mge_reads = MGE_SAM_ANALYSIS
        this.genes_length = getGenesLength(MEGARES)
        this.abundance_dict = {"Bovine":{}, "Human":{}, "Soil":{}, "Mock":{}}
        this.initial_source_size = {"Bovine":{}, "Human":{}, "Soil":{}, "Mock":{}}

    def findAbsoluteAbundance(this, fileName):
        sample, x_axis_name = x_axis(fileName)
        if x_axis_name not in list(this.abundance_dict[sample].keys()):
            this.abundance_dict[sample].update({x_axis_name:indiv_abundance(sample,x_axis_name)})
            this.initial_source_size[sample.update({x_axis_name:0})]
        this.abundance_dict[sample][x_axis_name].addToAbsolute(this.filePath(fileName, this.amr_reads, this.source_prefix))
        with gzip.open(this.filePath(fileName, "", this.initial_source_prefix), "rb") as input:
            with open ("temp_files/deduplicated_sequel-demultiplex.temp.ccs.fastq", "wb") as output:
                shutil.copyfileobj(input, output)
        file_size = os.stat("temp_files/deduplicated_sequel-demultiplex.temp.ccs.fastq").st_size
        this.initial_source_size[sample][x_axis_name] += file_size


    def makeAbundanceRelative(this):
        