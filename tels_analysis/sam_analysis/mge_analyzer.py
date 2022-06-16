import pysam

class mge_analyzer:
    def __init__(this, fileName, SOURCE_PREFIX, SOURCE_SUFFIX, A_TO_MGES):
        this.fileName = fileName
        this.filePath = SOURCE_PREFIX + fileName + SOURCE_SUFFIX + A_TO_MGES

    def mge_genes_list(this, OUTPUT_PREFIX, OUTPUT_SUFFIX, SAM_ANALYSIS):
        MGEstatFile = open(OUTPUT_PREFIX + this.fileName + OUTPUT_SUFFIX + "_mge" + SAM_ANALYSIS, "w")
        readDict = {}
        pysam.sort("-o", "temp_files/temp_sorted_sam_file.sam", this.filePath)
        samfile = pysam.AlignmentFile("temp_files/temp_sorted_sam_file.sam", "r")
        iter = samfile.fetch()
        for read in iter:
            mge = read.reference_name
            readName = read.query_name
            if mge == None:
                break
            flagInt = read.flag
            flagBinary = []
            for i in range(0,12):
                flagBinary.append(int(flagInt/2**(11-i)))
                flagInt %= 2**(11-i)
            if flagBinary[9] == 1:
                continue
            if flagBinary[0] == 1:
                if readDict.get(readName, False) == False:
                    readDict.update({readName:([],[],[mge])})
                else:
                    readDict[readName][2].append(mge)
            elif flagBinary[3] == 1:
                if readDict.get(readName, False) == False:
                    readDict.update({readName:([],[mge],[])})
                else:
                    readDict[readName][1].append(mge)
            else:
                if readDict.get(readName, False) == False:
                    readDict.update({readName:([mge],[],[])})
                else:
                    readDict[readName][0].append(mge)
        samfile.close()
        for read, mgeTuple in readDict.items():
            MGEstatFile.write(read)
            for mgeList in mgeTuple:
                for mge in mgeList:
                    MGEstatFile.write("," + mge)
                MGEstatFile.write(",")
            MGEstatFile.write("\n")
        MGEstatFile.close()
