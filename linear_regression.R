if (!require("rjson")) {
  install.packages("rjson")
}
if (!require("rgl")) {
  install.packages("rgl")
}
if (!require("nlme")) {
  install.packages("nlme")
}
if (!require("lme4")) {
  install.packages("lme4")
}
if (!require("RLRsim")) {
  install.packages("RLRsim")
}
if (!require("MuMIn")) {
  install.packages('MuMIn')
}

if (!require("dplyr")) {
  install.packages('dplyr')
}

library(rjson)
library(rgl)
library(nlme)
library(lme4)
library(RLRsim)
library(MuMIn)
library(dplyr)

#Creating data frame of results
sample_count = 96
df=data.frame(Samples=c("BFV2AA","BFV2AB","BFV2AC","BFV2AMA","BFV2AMB","BFV2AMC","BFV2MA","BFV2MB","BFV2MC","BFV2NEGA","BFV2NEGAM","BFV2NEGM",
                        "BFXTAA","BFXTAB","BFXTAC","BFXTAMA","BFXTAMB","BFXTAMC","BFXTMA","BFXTMB","BFXTMC",'BFXTNEGA',"BFXTNEGAM","BFXTNEGM",
                        "HFV2AA","HFV2AB","HFV2AC","HFV2AMA","HFV2AMB","HFV2AMC","HFV2MA","HFV2MB","HFV2MC","HFV2NEGA","HFV2NEGAM","HFV2NEGM",
                        "HFXTAA","HFXTAB","HFXTAC","HFXTAMA","HFXTAMB","HFXTAMC","HFXTMA","HFXTMB","HFXTMC","HFXTNEGA","HFXTNEGAM","HFXTNEGM",
                        "MOV2AA","MOV2AB","MOV2AC","MOV2AMA","MOV2AMB","MOV2AMC","MOV2MA","MOV2MB","MOV2MC","MOV2NEGA","MOV2NEGAM","MOV2NEGM",
                        "MOXTAA","MOXTAB","MOXTAC","MOXTAMA","MOXTAMB","MOXTAMC","MOXTMA","MOXTMB","MOXTMC","MOXTNEGA",'MOXTNEGAM',"MOXTNEGM",
                        "SV2AA","SV2AB","SV2AC","SV2AMA","SV2AMB","SV2AMC","SV2MA","SV2MB","SV2MC","SV2NEGA","SV2NEGAM","SV2NEGM",
                        "SXTAA","SXTAB","SXTAC","SXTAMA","SXTAMB","SXTAMC","SXTMA","SXTMB","SXTMC","SXTNEGA","SXTNEGAM","SXTNEGM"),
              Chemistry = integer(sample_count),
              SampleType = integer(sample_count),
              Probe = integer(sample_count),
              Read_Count = integer(sample_count),
              Total_Read_Length = integer(sample_count),
              Unique_Coloc = integer(sample_count),
              ARG = integer(sample_count),
              MGE = integer(sample_count),
              LogARG = integer(sample_count),
              LogMGE = integer(sample_count),
              Random = integer(sample_count),
              Random2 = integer(sample_count),
              Random3 = integer(sample_count))
              

#Filling data frame
ref_length <- read.csv(
  "~/Documents/GitHub/TELS_analysis/output/references_length.csv", header=FALSE)
for (i in 1:sample_count) {
  if ((i-1)%%(sample_count/4) > 11) df$Chemistry[i] <- 'XT'
  else  df$Chemistry[i] <- 'V2'
  
  if ((i-1)%/%(sample_count/4) == 0) df$SampleType[i] <- 'Bovine'
  else if ((i-1)%/%(sample_count/4) == 1) df$SampleType[i] <- 'Human'
  else if ((i-1)%/%(sample_count/4) == 2) df$SampleType[i] <- 'Mock'
  else df$SampleType[i] <- 'Soil'
  
  if ((i-1)%%3 == 0) df$Random[i] <- 'A'
  else if ((i-1)%%3 == 1) df$Random[i] <- 'B'
  else df$Random[i] <- 'C'
  df$Random2[i] <- df$Random[i]
  df$Random3[i] <- df$Random[i]
  
  if ((i-1)%%(sample_count/8) < 3) df$Probe[i] <- 'ARG'
  else if ((i-1)%%(sample_count/8) < 6) df$Probe[i] <- 'ARG-MGE'
  else if ((i-1)%%(sample_count/8) < 9) df$Probe[i] <- 'MGE'
  else {
    df$Probe[i] <- 'None'
    df$Random[i] <- df$Chemistry[i]
    if ((i-1)%%12 == 9) {
      df$Random2[i] <- 'ARG'
      df$Random3[i] <- 'A'
    }
    else if ((i-1)%%12 == 10) {
      df$Random2[i] <- 'ARG-MGE'
      df$Random3[i] <- 'B'
    }
    else if ((i-1)%%12 == 11) {
      df$Random2[i] <- 'MGE'
      df$Random3[i] <- 'C'
    }
    df$Chemistry[i] <- 'None'
  }  
  
  colocalizations_richness <- read.csv(
    paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", 
          ".ccs.fastq_deduplicated.fastq_colocalizations_richness.csv", 
          sep=df$Samples[i]), 
    header=FALSE, row.names=NULL, comment.char="#")
  stats <- read.csv(
    paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", 
          ".ccs.fastq_deduplicated.fastq_stats.csv", 
          sep=df$Samples[i]), 
    header=FALSE, row.names=NULL)
  arg <- read.csv(
    paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", 
          ".ccs.fastq_deduplicated.fastq_SHORT_amr_diversity.csv", 
          sep=df$Samples[i]), 
    row.names=NULL)
  mge <- read.csv(
    paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", 
          ".ccs.fastq_deduplicated.fastq_SHORT_mobilome.csv", 
          sep=df$Samples[i]), 
    row.names=NULL)
  readlength <- fromJSON(
    file = paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", 
                 ".ccs.fastq_deduplicated.fastq_reads_length.json", 
                 sep=df$Samples[i]))
  
  df$Read_Count[i] <- stats$V2[2]
  df$ARG[i] <- as.integer(arg$Statistics[1])
  df$MGE[i] <- as.integer(mge$Statistics[1])
  df$LogARG[i] <- log2(1+as.numeric(arg$Statistics[1]))
  df$LogMGE[i] <- log2(1+as.numeric(mge$Statistics[1]))
  df$Total_Read_Length[i] <- sum(unlist(readlength))
  df$Unique_Coloc[i] = as.integer(colocalizations_richness$V2[1])
}

lmer.step <- function(object, steps = 1000, current_step = 0){
  if (current_step == steps)
    return(object)
  fit_change <- drop1(object)
  fit_change <- arrange(fit_change, AIC)
  smallest <- rownames(fit_change)[1]
  if (smallest == "<none>")
    return(object)
  new <- update(object, paste(".~. -" , smallest))
  return(lmer.step(new, current_step = current_step+1))
}


#For control, use assigned chemistry as random variable
lmer.coloc_grand_fit <- lmer(Unique_Coloc ~ Probe 
                            * SampleType 
                            * Chemistry 
                            * Total_Read_Length
                            + (Probe:SampleType | Random),
                            data=df)
lmer.coloc_best_fit <- lmer.step(lmer.coloc_grand_fit)

r.squaredGLMM(lmer.coloc_grand_fit)
r.squaredGLMM(lmer.coloc_best_fit)
extractAIC(lmer.coloc_grand_fit)
extractAIC(lmer.coloc_best_fit)

summary(lmer.coloc_best_fit)
residuals(lmer.coloc_best_fit)[abs(residuals(lmer.coloc_best_fit)) >= .5]
ranef(lmer.coloc_best_fit)
plot(x = df$Unique_Coloc, 
     y = round(fitted(lmer.coloc_best_fit)),
     main = "Unique Colocalization",
     xlab = "Original Values",
     ylab = "Fitted Values")

lmer.read_grand_fit <- lmer(Read_Count ~ Probe 
                            * SampleType 
                            * Chemistry 
                            + (Probe:SampleType | Random),
                            data=df)
lmer.read_best_fit <- lmer.step(lmer.read_grand_fit)

r.squaredGLMM(lmer.read_grand_fit)
r.squaredGLMM(lmer.read_best_fit)
extractAIC(lmer.read_grand_fit)
extractAIC(lmer.read_best_fit)

summary(lmer.read_best_fit)
residuals(lmer.read_best_fit)[abs(residuals(lmer.read_best_fit)) >= 50000]
ranef(lmer.read_best_fit)
plot(x = df$Read_Count, 
     y = round(fitted(lmer.read_best_fit)),
     main = "Total Reads",
     xlab = "Original Values",
     ylab = "Fitted Values")

lmer.ARG_grand_fit <- lmer(ARG ~ Probe 
                            * SampleType 
                            * Chemistry 
                            * Total_Read_Length
                            + (Probe:SampleType | Random),
                            data=df)
lmer.ARG_best_fit <- lmer.step(lmer.ARG_grand_fit)

r.squaredGLMM(lmer.ARG_grand_fit)
r.squaredGLMM(lmer.ARG_best_fit)
extractAIC(lmer.ARG_grand_fit)
extractAIC(lmer.ARG_best_fit)

summary(lmer.ARG_best_fit)
residuals(lmer.ARG_best_fit)[abs(residuals(lmer.ARG_best_fit)) >= 5]
ranef(lmer.ARG_best_fit)
plot(x = df$ARG, 
     y = round(fitted(lmer.ARG_best_fit)),
     main = "ARG Richness",
     xlab = "Original Values",
     ylab = "Fitted Values")

lmer.MGE_grand_fit <- lmer(MGE ~ Probe 
                           * SampleType 
                           * Chemistry 
                           * Total_Read_Length
                           + (Probe:SampleType | Random),
                           data=df)
lmer.MGE_best_fit <- lmer.step(lmer.MGE_grand_fit)

r.squaredGLMM(lmer.MGE_grand_fit)
r.squaredGLMM(lmer.MGE_best_fit)
extractAIC(lmer.MGE_grand_fit)
extractAIC(lmer.MGE_best_fit)

summary(lmer.MGE_best_fit)
residuals(lmer.MGE_best_fit)[abs(residuals(lmer.MGE_best_fit)) >= 50]
ranef(lmer.MGE_best_fit)
plot(x = df$MGE, 
     y = round(fitted(lmer.MGE_best_fit)),
     main = "MGE Richness",
     xlab = "Original Values",
     ylab = "Fitted Values")

#For control, use assigned probe as random variable
lmer.coloc_grand_fit_2 <- lmer(Unique_Coloc ~ Probe 
                             * SampleType 
                             * Chemistry 
                             * Total_Read_Length
                             + (SampleType:Chemistry | Random2),
                             data=df)
lmer.coloc_best_fit_2 <- lmer.step(lmer.coloc_grand_fit_2)

r.squaredGLMM(lmer.coloc_grand_fit_2)
r.squaredGLMM(lmer.coloc_best_fit_2)
extractAIC(lmer.coloc_grand_fit_2)
extractAIC(lmer.coloc_best_fit_2)

summary(lmer.coloc_best_fit_2)
residuals(lmer.coloc_best_fit_2)[abs(residuals(lmer.coloc_best_fit_2)) >= .5]
ranef(lmer.coloc_best_fit_2)
plot(x = df$Unique_Coloc, 
     y = round(fitted(lmer.coloc_best_fit_2)),
     main = "Unique Colocalization",
     xlab = "Original Values",
     ylab = "Fitted Values")

lmer.read_grand_fit_2 <- lmer(Read_Count ~ Probe 
                            * SampleType 
                            * Chemistry 
                            + (Probe:SampleType | Random),
                            data=df)
lmer.read_best_fit_2 <- lmer.step(lmer.read_grand_fit_2)

r.squaredGLMM(lmer.read_grand_fit_2)
r.squaredGLMM(lmer.read_best_fit_2)
extractAIC(lmer.read_grand_fit_2)
extractAIC(lmer.read_best_fit_2)

summary(lmer.read_best_fit_2)
residuals(lmer.read_best_fit_2)[abs(residuals(lmer.read_best_fit_2)) >= 50000]
ranef(lmer.read_best_fit_2)
plot(x = df$Read_Count, 
     y = round(fitted(lmer.read_best_fit_2)),
     main = "Total Reads",
     xlab = "Original Values",
     ylab = "Fitted Values")

lmer.ARG_grand_fit_2 <- lmer(ARG ~ Probe 
                           * SampleType 
                           * Chemistry 
                           * Total_Read_Length
                           + (Probe:SampleType | Random),
                           data=df)
lmer.ARG_best_fit_2 <- lmer.step(lmer.ARG_grand_fit_2)

r.squaredGLMM(lmer.ARG_grand_fit_2)
r.squaredGLMM(lmer.ARG_best_fit_2)
extractAIC(lmer.ARG_grand_fit_2)
extractAIC(lmer.ARG_best_fit_2)

summary(lmer.ARG_best_fit_2)
residuals(lmer.ARG_best_fit_2)[abs(residuals(lmer.ARG_best_fit_2)) >= 5]
ranef(lmer.ARG_best_fit_2)
plot(x = df$ARG, 
     y = round(fitted(lmer.ARG_best_fit_2)),
     main = "ARG Richness",
     xlab = "Original Values",
     ylab = "Fitted Values")

lmer.MGE_grand_fit_2 <- lmer(MGE ~ Probe 
                           * SampleType 
                           * Chemistry 
                           * Total_Read_Length
                           + (Probe:SampleType | Random),
                           data=df)
lmer.MGE_best_fit_2 <- lmer.step(lmer.MGE_grand_fit_2)

r.squaredGLMM(lmer.MGE_grand_fit_2)
r.squaredGLMM(lmer.MGE_best_fit_2)
extractAIC(lmer.MGE_grand_fit_2)
extractAIC(lmer.MGE_best_fit_2)

summary(lmer.MGE_best_fit_2)
residuals(lmer.MGE_best_fit_2)[abs(residuals(lmer.MGE_best_fit_2)) >= 50]
ranef(lmer.MGE_best_fit_2)
plot(x = df$MGE, 
     y = round(fitted(lmer.MGE_best_fit_2)),
     main = "MGE Richness",
     xlab = "Original Values",
     ylab = "Fitted Values")
