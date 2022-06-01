from tels_analysis.colocalization_analysis import colocalization
from tels_analysis.reads import reads

class colocalization_analyzer:
    def __init__(this, fileName, SOURCE_PREFIX, SOURCE_SUFFIX, COLOCALIZATIONS, READS_LENGTH, COLOCALIZATION_ANALYSIS):
        readList = reads.reads(fileName, SOURCE_PREFIX, SOURCE_SUFFIX, READS_LENGTH)
        this.readsDict = {}
        this.colocalizationList = []
        colocalizationFile = open(SOURCE_PREFIX + fileName + SOURCE_SUFFIX + COLOCALIZATIONS, "r")
        i = 0
        for line in colocalizationFile:
            i += 1
            if i > 1:
                params = line.split(',')
                read_length = readList.getLength(params[0])
                ARG = params[1]
                ARG_start = params[2][:params[2].find(":")]
                ARG_end = params[2][params[2].find(":")+1:]
                MGE_list = params[3].split(';')
                MGE_pos = params[4].split(';')
                MGE_start_list = []
                MGE_end_list = []
                for pos in MGE_pos:
                    MGE_start_list.append(pos[:pos.find(":")])
                    MGE_end_list.append(pos[pos.find(":")+1:])
                KEGG_list = params[5].split(';')
                KEGG_pos = params[6].split(';')
                KEGG_start_list = []
                KEGG_end_list = []
                for pos in KEGG_pos:
                    KEGG_start_list.append(pos[:pos.find(":")])
                    KEGG_end_list.append(pos[pos.find(":")+1:])
                this.colocalizationList.append(colocalization.colocalization(read_length, ARG, ARG_start, ARG_end, MGE_list, MGE_start_list, MGE_end_list, KEGG_list, KEGG_start_list, KEGG_end_list))
