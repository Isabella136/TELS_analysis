from tels_analysis.colocalization_analysis import colocalization
from tels_analysis.reads import reads
import rpy2
from rpy2 import robjects

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
    
    def makeChart(this, fileName, OUTPUT_PREFIX, OUTPUT_SUFFIX, COLOCALIZATION_ANALYSIS):
        colocInfoList = []
        for coloc in this.colocalizationList:
            colocInfoList.append(coloc)
        notCompleted = len(colocInfoList)
        index = 0
        while notCompleted != 0:
            x_val = []
            y_val = []
            for i in range(0,len(colocInfoList)):
                if index >= len(colocInfoList[i]):
                    x_val.extend(-1,-1,-1,-1,-1)
                    y_val.append("empty")
                else:
                    x_val.extend(0,colocInfoList[i][index][1], colocInfoList[i][index][1] + 1, colocInfoList[i][index][2], colocInfoList[i][index][4])
                    y_val.append(colocInfoList[i][index][0])
        #Sample R code to modify in the future
        robjects.r('''
        install.packages("tidyverse")
        library(ggplot2)
        myColors <- c("0" = "red3", "1" = "gold", "2" = "deepskyblue", "3" = "springgreen1", "4" = "green4", "5" = "darkgreen")
        mydata<-data.frame(x=c(0,80,120,1600,10000,0,500,560,6200,13000,-1,-1,-1,-1,-1),y=rep(c('A','B','C'),each = 5))
        png(file="D:/GitHub/TELS_analysis/test.png", width = 1000, height = 1000)
        ggplot(mydata,aes(x=y,y=x,fill=factor(c))) + geom_boxplot(coef = 500, outlier.shape = NA, fatten = NULL, alpha=0.3) + scale_y_continuous(limits = c(0,20000)) + theme(legend.position="none", panel.background = element_rect(fill = "transparent")) + scale_fill_manual(values=myColors,breaks = waiver()) + coord_flip()
        dev.off()
        ''')
