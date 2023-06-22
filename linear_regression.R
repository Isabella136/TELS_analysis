if (!require("rjson")) {
  install.packages("rjson")
}
if (!require("lme4")) {
  install.packages("lme4")
}
if (!require("MuMIn")) {
  install.packages('MuMIn')
}
if (!require("dplyr")) {
  install.packages('dplyr')
}

library(rjson)
library(lme4)
library(MuMIn)
library(dplyr)

#Creating data frame of results
sample_count = 72
df=data.frame(Samples=c("BFV2AA","BFV2AB","BFV2AC","BFV2AMA","BFV2AMB","BFV2AMC","BFV2MA","BFV2MB","BFV2MC",
                        "BFXTAA","BFXTAB","BFXTAC","BFXTAMA","BFXTAMB","BFXTAMC","BFXTMA","BFXTMB","BFXTMC",
                        "HFV2AA","HFV2AB","HFV2AC","HFV2AMA","HFV2AMB","HFV2AMC","HFV2MA","HFV2MB","HFV2MC",
                        "HFXTAA","HFXTAB","HFXTAC","HFXTAMA","HFXTAMB","HFXTAMC","HFXTMA","HFXTMB","HFXTMC",
                        "MOV2AA","MOV2AB","MOV2AC","MOV2AMA","MOV2AMB","MOV2AMC","MOV2MA","MOV2MB","MOV2MC",
                        "MOXTAA","MOXTAB","MOXTAC","MOXTAMA","MOXTAMB","MOXTAMC","MOXTMA","MOXTMB","MOXTMC",
                        "SV2AA","SV2AB","SV2AC","SV2AMA","SV2AMB","SV2AMC","SV2MA","SV2MB","SV2MC",
                        "SXTAA","SXTAB","SXTAC","SXTAMA","SXTAMB","SXTAMC","SXTMA","SXTMB","SXTMC"),
              SampleID=c("BFV2A","BFV2A","BFV2A","BFV2AM","BFV2AM","BFV2AM","BFV2M","BFV2M","BFV2M",
                         "BFXTA","BFXTA","BFXTA","BFXTAM","BFXTAM","BFXTAM","BFXTM","BFXTM","BFXTM",
                         "HFV2A","HFV2A","HFV2A","HFV2AM","HFV2AM","HFV2AM","HFV2M","HFV2M","HFV2M",
                         "HFXTA","HFXTA","HFXTA","HFXTAM","HFXTAM","HFXTAM","HFXTM","HFXTM","HFXTM",
                         "MOV2A","MOV2A","MOV2A","MOV2AM","MOV2AM","MOV2AM","MOV2M","MOV2M","MOV2M",
                         "MOXTA","MOXTA","MOXTA","MOXTAM","MOXTAM","MOXTAM","MOXTM","MOXTM","MOXTM",
                         "SV2A","SV2A","SV2A","SV2AM","SV2AM","SV2AM","SV2M","SV2M","SV2M",
                         "SXTA","SXTA","SXTA","SXTAM","SXTAM","SXTAM","SXTM","SXTM","SXTM"),
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
              Random = integer(sample_count)
)
              

#Filling data frame
ref_length <- read.csv(
  "./output/references_length.csv", header=FALSE)
for (i in 1:sample_count) {
  if ((i-1)%%(sample_count/4) > 8) df$Chemistry[i] <- 'XT'
  else  df$Chemistry[i] <- 'V2'
  
  if ((i-1)%/%(sample_count/4) == 0) df$SampleType[i] <- 'Bovine'
  else if ((i-1)%/%(sample_count/4) == 1) df$SampleType[i] <- 'Human'
  else if ((i-1)%/%(sample_count/4) == 2) df$SampleType[i] <- 'Mock'
  else df$SampleType[i] <- 'Soil'
  
  if ((i-1)%%3 == 0) df$Random[i] <- 'A'
  else if ((i-1)%%3 == 1) df$Random[i] <- 'B'
  else df$Random[i] <- 'C'
  
  if ((i-1)%%(sample_count/8) < 3) df$Probe[i] <- 'ARG'
  else if ((i-1)%%(sample_count/8) < 6) df$Probe[i] <- 'ARG-MGE'
  else if ((i-1)%%(sample_count/8) < 9) df$Probe[i] <- 'MGE'
  
  colocalizations_richness <- read.csv(
    paste("./TELS_output/sequel-demultiplex.", 
          ".ccs.fastq_deduplicated.fastq_colocalizations_richness.csv", 
          sep=df$Samples[i]), 
    header=FALSE, row.names=NULL, comment.char="#")
  stats <- read.csv(
    paste("./TELS_output/sequel-demultiplex.", 
          ".ccs.fastq_deduplicated.fastq_stats.csv", 
          sep=df$Samples[i]), 
    header=FALSE, row.names=NULL)
  arg <- read.csv(
    paste("./TELS_output/sequel-demultiplex.", 
          ".ccs.fastq_deduplicated.fastq_SHORT_amr_diversity.csv", 
          sep=df$Samples[i]), 
    row.names=NULL)
  mge <- read.csv(
    paste("./TELS_output/sequel-demultiplex.", 
          ".ccs.fastq_deduplicated.fastq_SHORT_mobilome.csv", 
          sep=df$Samples[i]), 
    row.names=NULL)
  readlength <- fromJSON(
    file = paste("./TELS_output/sequel-demultiplex.", 
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

if (file.exists('./output/regression_results') == FALSE) {
  dir.create('./output/regression_results')
}

#For control, use assigned chemistry as random variable
lmer.coloc_grand_fit <- lmer(Unique_Coloc ~ Probe 
                            * SampleType 
                            * Chemistry 
                            * Total_Read_Length
                            + (1 | SampleID),
                            data=df)
lmer.coloc_best_fit <- lmer.step(lmer.coloc_grand_fit)

lmer.coloc_best_fit.rsquared <- r.squaredGLMM(lmer.coloc_best_fit)
lmer.coloc_best_fit.summary <- summary(lmer.coloc_best_fit)
lmer.coloc_best_fit.random_effects <- ranef(lmer.coloc_best_fit)
lmer.coloc_best_fit.coefficients <- lmer.coloc_best_fit.summary[["coefficients"]]
lmer.coloc_best_fit.confidence_intervals <- confint(lmer.coloc_best_fit, method='Wald')
coloc_coeff_amt = length(lmer.coloc_best_fit.coefficients)/3
df.lmer.coloc <- data_frame(
  row.names=row.names(lmer.coloc_best_fit.coefficients), 
  'Estimate' = lmer.coloc_best_fit.coefficients[
    c(1:coloc_coeff_amt)],
  'Standard Error' = lmer.coloc_best_fit.coefficients[
    c((coloc_coeff_amt+1):(2*coloc_coeff_amt))],
  't-Value' = lmer.coloc_best_fit.coefficients[
    c((2*coloc_coeff_amt+1):(coloc_coeff_amt*3))],
  '2.5%' = lmer.coloc_best_fit.confidence_intervals[
    c(3:(2+coloc_coeff_amt))],
  '97.5%' = lmer.coloc_best_fit.confidence_intervals[
    c((5+coloc_coeff_amt):(2*(coloc_coeff_amt+2)))])
write.csv(df.lmer.coloc,file='./output/regression_results/colocalization_coefficients.csv')
write.csv(lmer.coloc_best_fit.random_effects,file='./output/regression_results/coloc_rand_effects.csv')

png(file='./output/regression_results/colocalization_model_fit.png',
    width=600, height=600)
op <- par(mar = c(5,6,4,2) + 0.1)
plot(x = df$Unique_Coloc, 
     y = round(fitted(lmer.coloc_best_fit)),
     main = "Unique Colocalization",
     xlab = "Original Values",
     ylab = " ",
     xlim = c(-5, 50),
     ylim = c(-5, 50),
     las = 1)
title(ylab = "Fitted Values", cex.lab = 1, line = 4.5)
par(op)
dev.off()

lmer.read_grand_fit <- lmer(Read_Count ~ Probe 
                            * SampleType 
                            * Chemistry 
                            + (1 | SampleID),
                            data=df)
lmer.read_best_fit <- lmer.step(lmer.read_grand_fit)

lmer.read_best_fit.rsquared <- r.squaredGLMM(lmer.read_best_fit)
lmer.read_best_fit.summary <- summary(lmer.read_best_fit)
lmer.read_best_fit.random_effects <- ranef(lmer.read_best_fit)
lmer.read_best_fit.coefficients <- lmer.read_best_fit.summary[["coefficients"]]
lmer.read_best_fit.confidence_intervals <- confint(lmer.read_best_fit, method='Wald')
read_coeff_amt = length(lmer.read_best_fit.coefficients)/3
df.lmer.read <- data_frame(
  row.names=row.names(lmer.read_best_fit.coefficients), 
  'Estimate' = lmer.read_best_fit.coefficients[
    c(1:read_coeff_amt)],
  'Standard Error' = lmer.read_best_fit.coefficients[
    c((read_coeff_amt+1):(2*read_coeff_amt))],
  't-Value' = lmer.read_best_fit.coefficients[
    c((2*read_coeff_amt+1):(read_coeff_amt*3))],
  '2.5%' = lmer.read_best_fit.confidence_intervals[
    c(3:(2+read_coeff_amt))],
  '97.5%' = lmer.read_best_fit.confidence_intervals[
    c((5+read_coeff_amt):(2*(read_coeff_amt+2)))])
write.csv(df.lmer.read,file='./output/regression_results/read_counts_coefficients.csv')
write.csv(lmer.read_best_fit.random_effects,file='./output/regression_results/read_rand_effects.csv')

png(file='./output/regression_results/read_counts_model_fit.png',
    width=600, height=600)
op <- par(mar = c(5,6,4,2) + 0.1)
plot(x = df$Read_Count, 
     y = round(fitted(lmer.read_best_fit)),
     main = "Total Reads",
     xlab = "Original Values",
     ylab = " ",
     xlim = c(0, 800000),
     ylim = c(0, 800000),
     las = 1)
title(ylab = "Fitted Values", cex.lab = 1, line = 4.5)
par(op)
dev.off()

lmer.ARG_grand_fit <- lmer(ARG ~ Probe 
                            * SampleType 
                            * Chemistry 
                            * Total_Read_Length
                            + (1 | SampleID),
                            data=df)
lmer.ARG_best_fit <- lmer.step(lmer.ARG_grand_fit)

lmer.ARG_best_fit.rsquared <- r.squaredGLMM(lmer.ARG_best_fit)
lmer.ARG_best_fit.summary <- summary(lmer.ARG_best_fit)
lmer.ARG_best_fit.random_effects <- ranef(lmer.ARG_best_fit)
lmer.ARG_best_fit.coefficients <- lmer.ARG_best_fit.summary[["coefficients"]]
lmer.ARG_best_fit.confidence_intervals <- confint(lmer.ARG_best_fit, method='Wald')
ARG_coeff_amt = length(lmer.ARG_best_fit.coefficients)/3
df.lmer.ARG <- data_frame(
  row.names=row.names(lmer.ARG_best_fit.coefficients), 
  'Estimate' = lmer.ARG_best_fit.coefficients[
    c(1:ARG_coeff_amt)],
  'Standard Error' = lmer.ARG_best_fit.coefficients[
    c((ARG_coeff_amt+1):(2*ARG_coeff_amt))],
  't-Value' = lmer.ARG_best_fit.coefficients[
    c((2*ARG_coeff_amt+1):(ARG_coeff_amt*3))],
  '2.5%' = lmer.ARG_best_fit.confidence_intervals[
    c(3:(2+ARG_coeff_amt))],
  '97.5%' = lmer.ARG_best_fit.confidence_intervals[
    c((5+ARG_coeff_amt):(2*(ARG_coeff_amt+2)))])
write.csv(df.lmer.ARG,file='./output/regression_results/ARG_richness_coefficients.csv')
write.csv(lmer.ARG_best_fit.random_effects,file='./output/regression_results/ARG_rand_effects.csv')

png(file='./output/regression_results/ARG_richness_model_fit.png',
    width=600, height=600)
op <- par(mar = c(5,6,4,2) + 0.1)
plot(x = df$ARG, 
     y = round(fitted(lmer.ARG_best_fit)),
     main = "ARG Richness",
     xlab = "Original Values",
     ylab = " ",
     xlim = c(0, 120000),
     ylim = c(0, 120000),
     las = 1)
title(ylab = "Fitted Values", cex.lab = 1, line = 4.5)
par(op)
dev.off()

lmer.MGE_grand_fit <- lmer(MGE ~ Probe 
                           * SampleType 
                           * Chemistry 
                           * Total_Read_Length
                           + (1 | SampleID),
                           data=df)
lmer.MGE_best_fit <- lmer.step(lmer.MGE_grand_fit)

lmer.MGE_best_fit.rsquared <- r.squaredGLMM(lmer.MGE_best_fit)
lmer.MGE_best_fit.summary <- summary(lmer.MGE_best_fit)
lmer.MGE_best_fit.random_effects <- ranef(lmer.MGE_best_fit)
lmer.MGE_best_fit.coefficients <- lmer.MGE_best_fit.summary[["coefficients"]]
lmer.MGE_best_fit.confidence_intervals <- confint(lmer.MGE_best_fit, method='Wald')
MGE_coeff_amt = length(lmer.MGE_best_fit.coefficients)/3
df.lmer.MGE <- data_frame(
  row.names=row.names(lmer.MGE_best_fit.coefficients), 
  'Estimate' = lmer.MGE_best_fit.coefficients[
    c(1:MGE_coeff_amt)],
  'Standard Error' = lmer.MGE_best_fit.coefficients[
    c((MGE_coeff_amt+1):(2*MGE_coeff_amt))],
  't-Value' = lmer.MGE_best_fit.coefficients[
    c((2*MGE_coeff_amt+1):(MGE_coeff_amt*3))],
  '2.5%' = lmer.MGE_best_fit.confidence_intervals[
    c(3:(2+MGE_coeff_amt))],
  '97.5%' = lmer.MGE_best_fit.confidence_intervals[
    c((5+MGE_coeff_amt):(2*(MGE_coeff_amt+2)))])
write.csv(df.lmer.MGE,file='./output/regression_results/MGE_richness_coefficients.csv')
write.csv(lmer.MGE_best_fit.random_effects,file='./output/regression_results/MGE_rand_effects.csv')

png(file='./output/regression_results/MGE_richness_model_fit.png',
    width=600, height=600)
op <- par(mar = c(5,6,4,2) + 0.1)
plot(x = df$MGE, 
     y = round(fitted(lmer.MGE_best_fit)),
     main = "MGE Richness",
     xlab = "Original Values",
     ylab = " ",
     xlim = c(0, 120000),
     ylim = c(0, 120000),
     las = 1)
title(ylab = "Fitted Values", cex.lab = 1, line = 4.5)
par(op)
dev.off()
