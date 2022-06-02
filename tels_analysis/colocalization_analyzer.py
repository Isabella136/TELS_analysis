from tels_analysis.colocalization_analysis import colocalization
from tels_analysis.reads import reads
import rpy2
from rpy2 import robjects
from PIL import Image

class colocalization_analyzer:
    def __init__(this, fileName, SOURCE_PREFIX, SOURCE_SUFFIX, COLOCALIZATIONS, READS_LENGTH):
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
        robjects.r('''
            install.packages("tidyverse")
            library(ggplot2)
            r_makeChart <- function(x_val, y_val, col_val, fileName) {
                x_val_mean <-(x_val_start+x_val_end)/2
                mylabels<-data.frame(labels=c('A','B','C'),pos=x_val_mean)
                myColors <- c("0" = "red3", "1" = "gold", "2" = "deepskyblue", "3" = "springgreen1", "4" = "green4", "5" = "darkgreen")
                mydata<-data.frame(x=x_val,y=rep(y_val,each=5),c=rep(col_val,each=5))
                png(file=fileName, width = 10000, height = 10000)
                ggplot(mydata,aes(x=y,y=x,fill=factor(c))) + geom_boxplot(coef = 500, outlier.shape = NA, fatten = NULL, alpha=0.75) + geom_text(data=mylabels,aes(label=labels,y=pos,x=labels),position = position_dodge(width = .75),inherit.aes = FALSE) + scale_y_continuous(limits = c(0,20000)) + theme(legend.position="none", panel.background = element_rect(fill = "transparent")) + scale_fill_manual(values=myColors,breaks = waiver()) + coord_flip()
                dev.off()
            }
            ''')
        img = Image.new(mode = "RGB", size = (10000, 10000), color = (255, 255, 255))
        colocInfoList = []
        for coloc in this.colocalizationList:
            colocInfoList.append(coloc.getColocInfo())
        notCompleted = len(colocInfoList)
        index = 0
        while notCompleted != 0:
            x_val = []
            y_val = []
            col_val = []
            for i in range(0,len(colocInfoList)):
                if index >= len(colocInfoList[i]):
                    x_val.extend(-1,-1,-1,-1,-1)
                    y_val.append("empty")
                    col_val.append(0)
                else:
                    x_val.extend(0,colocInfoList[i][index][1], colocInfoList[i][index][1] + 1, colocInfoList[i][index][2], colocInfoList[i][index][4])
                    y_val.append(colocInfoList[i][index][0])
                    col_val.append(colocInfoList[i][index][3])
                    if index == (len(colocInfoList[i]) - 1):
                        notCompleted -= 1
            index += 1
            r_x_val = robjects.IntVector(x_val)
            r_y_val = robjects.StrVector(y_val)
            r_col_val = robjects.IntVector(col_val)
            tempPNGfile = "temp_file/" + fileName + "_" + index + ".png"
            robjects.globalenv["r_makeChart"](r_x_val, r_y_val, r_col_val, tempPNGfile)
            tempIMG = Image.open(tempPNGfile)
            for w in range(0, 10000):
                for h in range(0, 10000):
                    tempPixel = tempIMG.getpixel((w,h))
                    if tempPixel != (255,255,255):
                        img.putpixel((w,h),tempPixel)
            tempIMG.close()
        img.save( OUTPUT_PREFIX + fileName + OUTPUT_SUFFIX + COLOCALIZATION_ANALYSIS )
