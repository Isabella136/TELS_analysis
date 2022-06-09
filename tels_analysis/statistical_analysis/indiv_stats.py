from tels_analysis import fileDict

class indiv_stats:
    def __init__ (this, fileName):
        this.sample, this.seqPlatform = fileDict(fileName)
        this.raw_reads = 0
        this.deduplicated_reads = 0
        this.duplication = 0
        this.ARG_on_target = 0
        this.MGE_on_target = 0

    def getStats(this):
        toReturn = this.sample
        toReturn = toReturn + "," + this.seqPlatform
        toReturn = toReturn + "," + str(this.raw_reads)
        toReturn = toReturn + "," + str(this.deduplicated_reads)
        toReturn = toReturn + "," + str(this.duplication)
        toReturn = toReturn + "," + str(this.ARG_on_target)
        toReturn = toReturn + "," + str(this.MGE_on_target)
        return toReturn

    def findAllStats(this, filePath_stats, filePath_ARG, filePath_MGE):
        this.findReadStats(filePath_stats)
        this.findARGStats(filePath_ARG)
        this.findMGEStats(filePath_MGE)

    def findReadStats(this, filePath):
        statFile = open(filePath, "r")
        statFile.readline()
        this.raw_reads = int(statFile.readline().split(',')[1])
        this.deduplicated_reads = int(statFile.readline().split(',')[1])
        statFile.close()
        if this.deduplicated_reads == 0:
            this.deduplicated_reads = "__"
            this.duplication = "__"
        else:
            this.duplication = round(100 - ((this.deduplicated_reads/this.raw_reads) * 100), 1)
    def findARGStats(this, filePath):
        ARGstatFile = open(filePath, "r")
        ARGstatFile.readline()
        this.ARG_on_target = round((int(ARGstatFile.readline().split(',')[1]) / this.raw_reads) * 100,1)
        ARGstatFile.close()

    def findMGEStats(this, filePath):
        MGEstatFile = open(filePath, "r")
        MGEstatFile.readline()
        this.MGE_on_target = round((int(MGEstatFile.readline().split(',')[1]) / this.raw_reads) * 100, 1)
        MGEstatFile.close()