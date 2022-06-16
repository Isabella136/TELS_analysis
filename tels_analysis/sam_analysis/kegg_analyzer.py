import pysam

class kegg_analyzer:
    def __init__(this, fileName, SOURCE_PREFIX, SOURCE_SUFFIX, A_TO_KEGG):
        this.fileName = fileName
        this.filePath = SOURCE_PREFIX + fileName + SOURCE_SUFFIX + A_TO_KEGG

    def kegg_genes_list(this, OUTPUT_PREFIX, OUTPUT_SUFFIX, SAM_ANALYSIS):
        KEGGstatFile = open(OUTPUT_PREFIX + this.fileName + OUTPUT_SUFFIX + "_kegg" + SAM_ANALYSIS, "w")
        readDict = {}
        pysam.sort("-o", "temp_files/temp_sorted_sam_file.sam", this.filePath)
        samfile = pysam.AlignmentFile("temp_files/temp_sorted_sam_file.sam", "r")
        iter = samfile.fetch()
        for read in iter:
            kegg = read.reference_name
            readName = read.query_name
            if kegg == None:
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
                    readDict.update({readName:([],[],[kegg])})
                else:
                    readDict[readName][2].append(kegg)
            elif flagBinary[3] == 1:
                if readDict.get(readName, False) == False:
                    readDict.update({readName:([],[kegg],[])})
                else:
                    readDict[readName][1].append(kegg)
            else:
                if readDict.get(readName, False) == False:
                    readDict.update({readName:([kegg],[],[])})
                else:
                    readDict[readName][0].append(kegg)
        samfile.close()
        for read, keggTuple in readDict.items():
            KEGGstatFile.write(read)
            for keggList in keggTuple:
                for kegg in keggList:
                    KEGGstatFile.write("," + kegg)
                KEGGstatFile.write(",")
            KEGGstatFile.write("\n")
        KEGGstatFile.close()
