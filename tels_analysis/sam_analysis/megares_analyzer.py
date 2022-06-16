import pysam

class megares_analyzer:
    def __init__(this, fileName, SOURCE_PREFIX, SOURCE_SUFFIX, A_TO_MEGARES):
        this.fileName = fileName
        this.filePath = SOURCE_PREFIX + fileName + SOURCE_SUFFIX + A_TO_MEGARES

    def megares_genes_list(this, OUTPUT_PREFIX, OUTPUT_SUFFIX, SAM_ANALYSIS):
        ARGstatFile = open(OUTPUT_PREFIX + this.fileName + OUTPUT_SUFFIX + "_megares" + SAM_ANALYSIS, "w")
        readDict = {}
        pysam.sort("-o", "temp_files/temp_sorted_sam_file.sam", this.filePath)
        samfile = pysam.AlignmentFile("temp_files/temp_sorted_sam_file.sam", "r")
        iter = samfile.fetch()
        for read in iter:
            arg = read.reference_name
            readName = read.query_name
            if arg == None:
                break
            elif "RequiresSNPConfirmation" in arg:
                continue
            flagInt = read.flag
            flagBinary = []
            for i in range(0,12):
                flagBinary.append(int(flagInt/2**(11-i)))
                flagInt %= 2**(11-i)
            if flagBinary[9] == 1:
                continue
            if flagBinary[0] == 1:
                if readDict.get(readName, False) == False:
                    readDict.update({readName:([],[],[arg])})
                else:
                    readDict[readName][2].append(arg)
            elif flagBinary[3] == 1:
                if readDict.get(readName, False) == False:
                    readDict.update({readName:([],[arg],[])})
                else:
                    readDict[readName][1].append(arg)
            else:
                if readDict.get(readName, False) == False:
                    readDict.update({readName:([arg],[],[])})
                else:
                    readDict[readName][0].append(arg)
        samfile.close()
        for read, argTuple in readDict.items():
            ARGstatFile.write(read)
            for argList in argTuple:
                for arg in argList:
                    ARGstatFile.write("," + arg)
                ARGstatFile.write(",")
            ARGstatFile.write("\n")
        ARGstatFile.close()

