if (!require("rjson")) {
  install.packages("rjson")
}
if (!require("rgl")) {
  install.packages("rgl")
}
library(rjson)
library(rgl)
df=data.frame(Samples=c("BFV2AA","BFV2AB","BFV2AC","BFV2AMA","BFV2AMB","BFV2AMC","BFV2MA","BFV2MB","BFV2MC","BFV2NEGA","BFV2NEGAM","BFV2NEGM",
                        "BFXTAA","BFXTAB","BFXTAC","BFXTAMA","BFXTAMB","BFXTAMC","BFXTMA","BFXTMB","BFXTMC",'BFXTNEGA',"BFXTNEGAM","BFXTNEGM",
                        "HFV2AA", "HFV2AB", "HFV2AC","HFV2AMA","HFV2AMB","HFV2AMC","HFV2MA","HFV2MB","HFV2MC","HFV2NEGA","HFV2NEGAM","HFV2NEGM",
                        "HFXTAA", "HFXTAB", "HFXTAC","HFXTAMA","HFXTAMB","HFXTAMC","HFXTMA","HFXTMB","HFXTMC","HFXTNEGA","HFXTNEGAM","HFXTNEGM",
                        "MOV2AA", "MOV2AB", "MOV2AC","MOV2AMA","MOV2AMB","MOV2AMC","MOV2MA","MOV2MB","MOV2MC","MOV2NEGA","MOV2NEGAM","MOV2NEGM",
                        "MOXTAA", "MOXTAB", "MOXTAC","MOXTAMA","MOXTAMB","MOXTAMC","MOXTMA","MOXTMB","MOXTMC","MOXTNEGA",'MOXTNEGAM',"MOXTNEGM",
                        "SV2AA", "SV2AB", "SV2AC","SV2AMA","SV2AMB","SV2AMC","SV2MA","SV2MB","SV2MC","SV2NEGA","SV2NEGAM","SV2NEGM",
                        "SXTAA", "SXTAB", "SXTAC","SXTAMA","SXTAMB","SXTAMC","SXTMA","SXTMB","SXTMC","SXTNEGA","SXTNEGAM","SXTNEGM"),
              Unique_Coloc = integer(sample_count),
              XTvsV2 = integer(sample_count),
              SampleType = integer(sample_count),
              Read_Count = integer(sample_count),
              ARG = integer(sample_count),
              MGE = integer(sample_count),
              LogARG = integer(sample_count),
              LogMGE = integer(sample_count),
              Read_Length = integer(sample_count),
              Classified = integer(sample_count),
              Classified_Aligned = integer(sample_count),
              Aligned = integer(sample_count))
sample_count = 96
mobilome_info <- read.csv("~/Documents/GitHub/TELS_analysis/output/mobilome_info.csv")
for (i in 1:sample_count) {
  if ((i-1)%%24 > 11) {
    df$XTvsV2[i] <- 'XT'
  }
  else {
    df$XTvsV2[i] <- 'V2'
  }
  if ((i-1)%/%24 == 0) {
    df$SampleType[i] <- 'Bovine'
  }
  else if ((i-1)%/%24 == 1) {
    df$SampleType[i] <- 'Human'
  }
  else if ((i-1)%/%24 == 2) {
    df$SampleType[i] <- 'Mock'
  }
  else {
    df$SampleType[i] <- 'Soil'
  }
  if ((i-1)%%12 < 3) {
    df$Probe[i] <- 'ARG'
  }
  else if ((i-1)%%12 < 6) {
    df$Probe[i] <- 'ARG-MGE'
  }
  else if ((i-1)%%12 < 9) {
    df$Probe[i] <- 'MGE'
  }
  else {
    df$Probe[i] <- 'None'
  }
  colocalizations_richness <- read.csv(paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", ".ccs.fastq_deduplicated.fastq_colocalizations_richness.csv", sep=df$Samples[i]), header=FALSE, row.names=NULL, comment.char="#")
  stats <- read.csv(paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", ".ccs.fastq_deduplicated.fastq_stats.csv", sep=df$Samples[i]), header=FALSE, row.names=NULL)
  arg <- read.csv(paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", ".ccs.fastq_deduplicated.fastq_SHORT_amr_diversity.csv", sep=df$Samples[i]), row.names=NULL)
  mge <- read.csv(paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", ".ccs.fastq_deduplicated.fastq_SHORT_mobilome.csv", sep=df$Samples[i]), row.names=NULL)
  readlength <- fromJSON(file = paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", ".ccs.fastq_deduplicated.fastq_reads_length.json", sep=df$Samples[i]))
  df$Read_Count[i] <- stats$V2[2]
  df$ARG[i] <- arg$Statistics[1]
  df$MGE[i] <- mge$Statistics[1]
  df$LogARG[i] <- log2(1+as.numeric(arg$Statistics[1]))
  df$LogMGE[i] <- log2(1+as.numeric(mge$Statistics[1]))
  df$Read_Length[i] <- sum(unlist(readlength))/stats$V2[2]
  df$Unique_Coloc[i] = colocalizations_richness$V2[1]
  df$Classified[i] <- mobilome_info[105, i+1]
  df$Classified_Aligned[i] <- mobilome_info[151, i+1]
  df$Aligned[i] <- mobilome_info[283, i+1]
}
fit_tels_arg <- lm(data=df, Unique_Coloc ~ as.numeric(ARG))
summary(fit_tels_arg) # model summary
coefficients(fit_tels_arg) # model coefficients
fitted(fit_tels_arg) # predicted values
residuals(fit_tels_arg) # model residuals
anova(fit_tels_arg) # anova table, this is typically what we report
plot(fit_tels_arg)

fit_tels_arg_mge <- lm(data=df, Unique_Coloc ~ as.numeric(ARG) + as.numeric(MGE))
summary(fit_tels_arg_mge) # model summary
coefficients(fit_tels_arg_mge) # model coefficients
fitted(fit_tels_arg_mge) # predicted values
residuals(fit_tels_arg_mge) # model residuals
anova(fit_tels_arg_mge) # anova table, this is typically what we report
plot(fit_tels_arg_mge)

fit_tels_arg_classifiedalignedarg <- lm(data=df, Unique_Coloc ~ as.numeric(ARG) + Classified_Aligned)
summary(fit_tels_arg_classifiedalignedarg) # model summary
coefficients(fit_tels_arg_classifiedalignedarg) # model coefficients
fitted(fit_tels_arg_classifiedalignedarg) # predicted values
residuals(fit_tels_arg_classifiedalignedarg) # model residuals
anova(fit_tels_arg_classifiedalignedarg) # anova table, this is typically what we report
plot(fit_tels_arg_classifiedalignedarg)

fit_tels_classifiedalignedarg <- lm(data=df, Unique_Coloc ~ Classified_Aligned)
summary(fit_tels_classifiedalignedarg) # model summary
coefficients(fit_tels_classifiedalignedarg) # model coefficients
fitted(fit_tels_classifiedalignedarg) # predicted values
residuals(fit_tels_classifiedalignedarg) # model residuals
anova(fit_tels_classifiedalignedarg) # anova table, this is typically what we report
plot(fit_tels_classifiedalignedarg)

y_arg = integer(120001)
for (i in 1:120001) {
  y_arg[i] = fit_tels_arg$coefficients[1] + (i-1) * fit_tels_arg$coefficients[2]
}
plot(as.numeric(df$ARG), df$Unique_Coloc, main = "ARG Linear Regression", xlab = "ARG count", ylab = "Unique Colocalization")
lines(y_arg, type="l")

y_ca = integer(81)
for (i in 1:81) {
  y_ca[i] = fit_tels_classifiedalignedarg$coefficients[1] + (i-1) * fit_tels_classifiedalignedarg$coefficients[2]
}
plot(as.numeric(df$Classified_Aligned), df$Unique_Coloc, main = "MGE Classified as and Aligned to ARG Linear Regression", xlab = "MGE count", ylab = "Unique Colocalization")
lines(y_ca, type="l")

z_arg_ca = integer(81)
x_arg_ca = 0:80
y_arg_ca = seq(0, 120000, by=1500)
for (i in 1:81) {
  z_arg_ca[i] = fit_tels_arg_classifiedalignedarg$coefficients[1] + y_arg_ca[i] * fit_tels_arg_classifiedalignedarg$coefficients[2] + x_arg_ca[i] * fit_tels_arg_classifiedalignedarg$coefficients[3]
}
plot3d(x=as.numeric(df$Classified_Aligned), y=as.numeric(df$ARG), z=df$Unique_Coloc,
       main = "ARG vs MGE Classified as and Aligned to ARG Linear Regression", 
       xlab = "MGE count",
       ylab = "ARG count", 
       zlab = "Unique Colocalization")
lines3d(x=x_arg_ca, y=y_arg_ca, z=z_arg_ca)
