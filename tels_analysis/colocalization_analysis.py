from tels_analysis import colocalization_analysis
from colocalization_analysis import colocalization

class colocalization_analysis:
    def __init__(this, filename):
        this.readsDict = {}
        colocalizationFile = open("deduplicated_sequel-demultiplex." + filename + "ccs.fastq.gz_colocalizations.csv", "r")
