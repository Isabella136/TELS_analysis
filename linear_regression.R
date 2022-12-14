if (!require("rjson")) {
  install.packages("rjson")
}
library(rjson)
df=data.frame(Samples=c("BFV2AA","BFV2AB","BFV2AC","BFV2AMA","BFV2AMB","BFV2AMC","BFV2MA","BFV2MB","BFV2MC","BFV2NEGA","BFV2NEGAM","BFV2NEGM",
                        "BFXTAA","BFXTAB","BFXTAC","BFXTAMA","BFXTAMB","BFXTAMC","BFXTMA","BFXTMB","BFXTMC",'BFXTNEGA',"BFXTNEGAM","BFXTNEGM",
                        "HFV2AA", "HFV2AB", "HFV2AC","HFV2AMA","HFV2AMB","HFV2AMC","HFV2MA","HFV2MB","HFV2MC","HFV2NEGA","HFV2NEGAM","HFV2NEGM",
                        "HFXTAA", "HFXTAB", "HFXTAC","HFXTAMA","HFXTAMB","HFXTAMC","HFXTMA","HFXTMB","HFXTMC","HFXTNEGA","HFXTNEGAM","HFXTNEGM",
                        "MOV2AA", "MOV2AB", "MOV2AC","MOV2AMA","MOV2AMB","MOV2AMC","MOV2MA","MOV2MB","MOV2MC","MOV2NEGA","MOV2NEGAM","MOV2NEGM",
                        "MOXTAA", "MOXTAB", "MOXTAC","MOXTAMA","MOXTAMB","MOXTAMC","MOXTMA","MOXTMB","MOXTMC","MOXTNEGA",'MOXTNEGAM',"MOXTNEGM",
                        "SV2AA", "SV2AB", "SV2AC","SV2AMA","SV2AMB","SV2AMC","SV2MA","SV2MB","SV2MC","SV2NEGA","SV2NEGAM","SV2NEGM",
                        "SXTAA", "SXTAB", "SXTAC","SXTAMA","SXTAMB","SXTAMC","SXTMA","SXTMB","SXTMC","SXTNEGA","SXTNEGAM","SXTNEGM"))
sample_count = 96
df <- cbind(df, Unique_Coloc = integer(sample_count))
for (i in 1:sample_count) {
  colocalizations_richness <- read.csv(paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", ".ccs.fastq_deduplicated.fastq_colocalizations_richness.csv", sep=df$Samples[i]), header=FALSE, row.names=NULL, comment.char="#")
  df$Unique_Coloc[i] = colocalizations_richness$V2[1]
}
df <- cbind(df, XTvsV2 = integer(sample_count))
df <- cbind(df, SampleType = integer(sample_count))
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
}
fit_tels <- lm(data=df, Unique_Coloc ~ as.factor(SampleType) * as.factor(XTvsV2) * as.factor(Probe))
summary(fit_tels) # model summary
coefficients(fit_tels) # model coefficients
confint(fit_tels, level = 0.95) # CIs for model parameters
fitted(fit_tels) # predicted values
residuals(fit_tels) # model residuals
#vcov(fit_tels) # covariance matrix for model parameters
influence(fit_tels) # regression diagnostics
anova(fit_tels) # anova table, this is typically what we report
plot(fit_tels)
df <- cbind(df, Read_Count = integer(sample_count))
df <- cbind(df, ARG = integer(sample_count))
df <- cbind(df, MGE = integer(sample_count))
df <- cbind(df, Read_Length = integer(sample_count))
for (i in 1:sample_count) {
  stats <- read.csv(paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", ".ccs.fastq_deduplicated.fastq_stats.csv", sep=df$Samples[i]), header=FALSE, row.names=NULL)
  arg <- read.csv(paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", ".ccs.fastq_deduplicated.fastq_SHORT_amr_diversity.csv", sep=df$Samples[i]), row.names=NULL)
  mge <- read.csv(paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", ".ccs.fastq_deduplicated.fastq_SHORT_mobilome.csv", sep=df$Samples[i]), row.names=NULL)
  readlength <- fromJSON(file = paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", ".ccs.fastq_deduplicated.fastq_reads_length.json", sep=df$Samples[i]))
  df$Read_Count[i] <- stats$V2[2]
  df$ARG[i] <- arg$Statistics[1]
  df$MGE[i] <- mge$Statistics[1]
  df$Read_Length[i] <- sum(unlist(readlength))/stats$V2[2]
}
fit_tels <- lm(data=df, Unique_Coloc ~ Read_Count * Read_Length)
summary(fit_tels) # model summary
coefficients(fit_tels) # model coefficients
confint(fit_tels, level = 0.95) # CIs for model parameters
fitted(fit_tels) # predicted values
residuals(fit_tels) # model residuals
#vcov(fit_tels) # covariance matrix for model parameters
influence(fit_tels) # regression diagnostics
anova(fit_tels) # anova table, this is typically what we report
plot(fit_tels)