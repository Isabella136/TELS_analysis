if (!require("rjson"))
  install.packages("rjson")

if (!require("lme4"))
  install.packages("lme4")

if (!require("MuMIn"))
  install.packages('MuMIn')

if (!require("dplyr"))
  install.packages('dplyr')

library(rjson)
library(lme4)
library(MuMIn)
library(dplyr)

# Filling data frame function
fill_data_frame <- function(df, count, df_organisms_count, organism=FALSE){
  for (index in 1:count) {
    # Chemistry
    if ((index-1)%%(count/df_organisms_count) > 8)
      df$Chemistry[index] <- 'XT'
    else
      df$Chemistry[index] <- 'V2'
    
    # Organisms
    if (organism) {
      if ((index-1)%/%(count/df_organisms_count) == 0)
        df$SampleType[index] <- 'Bovine'
      else if ((index-1)%/%(count/df_organisms_count) == 1)
        df$SampleType[index] <- 'Human'
      else if (((index-1)%/%(count/df_organisms_count) == 2) & df_organisms_count > 3)
        df$SampleType[index] <- 'Mock'
      else
        df$SampleType[index] <- 'Soil'
    }
    
    # Probe Specificity
    if ((index-1)%%(count/(2*df_organisms_count)) < 3)
      df$Probe[index] <- 'ARG'
    else if ((index-1)%%(count/(2*df_organisms_count)) < 6)
      df$Probe[index] <- 'ARG-MGE'
    else
      df$Probe[index] <- 'MGE'
    
    # Colocalization Richness
    colocalizations_richness <- read.csv(
      paste("./TELS_output/sequel-demultiplex.", 
            ".ccs.fastq_deduplicated.fastq_colocalizations_richness.csv", 
            sep=df$Samples[index]), 
      header=FALSE, row.names=NULL, comment.char="#")
    df$Unique_Coloc[index] = as.integer(colocalizations_richness$V2[1])
    
    # Dedup Read Count
    stats <- read.csv(
      paste("./TELS_output/sequel-demultiplex.", 
            ".ccs.fastq_deduplicated.fastq_stats.csv", 
            sep=df$Samples[index]), 
      header=FALSE, row.names=NULL)
    df$Read_Count[index] <- stats$V2[2]
    
    # ARG Richness
    arg <- read.csv(
      paste("./TELS_output/sequel-demultiplex.", 
            ".ccs.fastq_deduplicated.fastq_SHORT_amr_diversity.csv", 
            sep=df$Samples[index]), 
      row.names=NULL)
    df$ARG[index] <- as.integer(arg$Statistics[1])
    
    # MGE Richness
    mge <- read.csv(
      paste("./TELS_output/sequel-demultiplex.", 
            ".ccs.fastq_deduplicated.fastq_SHORT_mobilome.csv", 
            sep=df$Samples[index]), 
      row.names=NULL)
    df$MGE[index] <- as.integer(mge$Statistics[1])
    
    # Total Read Length
    readlength <- fromJSON(
      file = paste("./TELS_output/sequel-demultiplex.", 
                   ".ccs.fastq_deduplicated.fastq_reads_length.json", 
                   sep=df$Samples[index]))
    df$Total_Read_Length[index] <- sum(unlist(readlength))
    
  }
  return(df)
}

# Step-wise regression
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

# Linear Mixed-Effect Regression for colocalization
lmer_coloc <- function(df, organism = FALSE) {
  lmer.coloc_grand_fit <- {
    if (organism) {
      lmer(
        Unique_Coloc ~  
          Probe + SampleType + Chemistry + Total_Read_Length 
        + Probe:SampleType + Probe:Chemistry + Probe:Total_Read_Length 
        + SampleType:Chemistry + SampleType:Total_Read_Length 
        + Total_Read_Length:Chemistry + (1 | SampleID),
        data=df)
    }
    else {
      lmer(
        Unique_Coloc ~ 
          Probe + Chemistry + Total_Read_Length + Probe:Total_Read_Length 
        + Chemistry:Total_Read_Length + (1 | SampleID),
        data=df)
    }
  }
  return(lmer.step(lmer.coloc_grand_fit))
}

# Linear Mixed-Effect Regression for ARG
lmer_arg <- function(df, organism = FALSE) {
  lmer.arg_grand_fit <- {
    if (organism) {
      lmer(
        ARG ~ 
          Probe + SampleType + Chemistry + Total_Read_Length 
        + Probe:SampleType + Probe:Chemistry + Probe:Total_Read_Length 
        + SampleType:Chemistry + SampleType:Total_Read_Length 
        + Total_Read_Length:Chemistry + (1 | SampleID),
        data=df)
    }
    else {
      lmer(
        ARG ~ 
          Probe + Chemistry + Total_Read_Length + Probe:Total_Read_Length 
        + Chemistry:Total_Read_Length + (1 | SampleID),
        data=df)
    }
  }
  return(lmer.step(lmer.arg_grand_fit))
}

# Linear Mixed-Effect Regression for MGE
lmer_mge <- function(df, organism = FALSE) {
  lmer.mge_grand_fit <- {
    if (organism) {
      lmer(
        MGE ~  
          Probe + SampleType + Chemistry + Total_Read_Length 
        + Probe:SampleType + Probe:Chemistry + Probe:Total_Read_Length 
        + SampleType:Chemistry + SampleType:Total_Read_Length 
        + Total_Read_Length:Chemistry + (1 | SampleID),
        data=df)
    }
    else {
      lmer(
        MGE ~ 
          Probe + Chemistry + Total_Read_Length + Probe:Total_Read_Length 
        + Chemistry:Total_Read_Length + (1 | SampleID),
        data=df)
    }
  }
  return(lmer.step(lmer.mge_grand_fit))
}

# Linear Mixed-Effect Regression for read counts
lmer_read <- function(df, organism = FALSE) {
  lmer.read_grand_fit <- {
    if (organism) {
      lmer(
        Read_Count ~  
          Probe + SampleType + Chemistry + Probe:SampleType 
        + Probe:Chemistry + SampleType:Chemistry + (1 | SampleID),
        data=df)
    }
    else {
      lmer(
        Read_Count ~ Probe * Chemistry + (1 | SampleID),
        data=df)
    }
  }
  return(lmer.step(lmer.read_grand_fit))
}

plot_coloc <- function(df, lmer) {
  plot(x = df$Unique_Coloc, 
       y = round(fitted(lmer)),
       main = "Unique Colocalization",
       xlab = "Original Values",
       ylab = " ",
       xlim = c(-5, 50),
       ylim = c(-5, 50),
       las = 1)
}

plot_arg  <- function(df, lmer) {
  plot(x = df$ARG, 
       y = round(fitted(lmer)),
       main = "ARG Richness",
       xlab = "Original Values",
       ylab = " ",
       xlim = c(0, 120000),
       ylim = c(0, 120000),
       las = 1)
}

plot_mge <- function(df, lmer) {
  plot(x = df$MGE, 
       y = round(fitted(lmer)),
       main = "MGE Richness",
       xlab = "Original Values",
       ylab = " ",
       xlim = c(0, 120000),
       ylim = c(0, 120000),
       las = 1)
}

plot_read <- function(df, lmer) {
  plot(x = df$Read_Count, 
       y = round(fitted(lmer)),
       main = "Total Reads",
       xlab = "Original Values",
       ylab = " ",
       xlim = c(0, 800000),
       ylim = c(0, 800000),
       las = 1)
}

# Write LMER information
write_lmer <- function(lmer_function, lmer_filename, lmer_org_info, 
                       plot_function, df, organism = FALSE) {
  lmer <- lmer_function(df, organism)
  lmer.summary <- summary(lmer)
  lmer.random_effects <- ranef(lmer)
  lmer.rsquared <- r.squaredGLMM(lmer)
  lmer.coefficients <- lmer.summary[["coefficients"]]
  lmer.confidence_intervals <- confint(lmer, method='Wald')
  coeff_amt = length(lmer.coefficients) / 3
  df.lmer <- data_frame(
    row.names = row.names(lmer.coefficients),
    'Estimate' = lmer.coefficients[c(1:coeff_amt)],
    'Standard Error' = lmer.coefficients[c((coeff_amt+1):(2*coeff_amt))],
    't-Value' = lmer.coefficients[c((2*coeff_amt+1):(coeff_amt*3))],
    '2.5%' = lmer.confidence_intervals[c(3:(2+coeff_amt))],
    '97.5%' = lmer.confidence_intervals[c((5+coeff_amt):(2*(coeff_amt+2)))])
  write.csv(df.lmer, file=paste(
    output_folder, lmer_org_info, lmer_filename, 'coefficients.csv'))
  write.csv(lmer.random_effects, file=paste(
    output_folder, lmer_org_info, lmer_filename, 'rand_effects.csv'))
  write.csv(lmer.rsquared, file=paste(
    output_folder, lmer_org_info, lmer_filename, 'rsquared.csv'))
  
  png(file=paste(output_folder, lmer_org_info, lmer_filename, 'model_fit.png'),
      width=600, height=600)
  op <- par(mar = c(5,6,4,2) + 0.1)
  plot_function(df, lmer)
  title(ylab = "Fitted Values", cex.lab = 1, line = 4.5)
  par(op)
  dev.off()
  data = data.frame(row.names = paste(lmer_org_info, lmer_filename),
                    'AIC' = extractAIC(lmer)[2],
                    'edf' = extractAIC(lmer)[1])
  browser()
  write.table(data, file=paste(output_folder, 'AIC.csv'), append=TRUE)
  return(lmer.summary[["AICtab"]][["REML"]])
}

# Creating data frame of results
mock_count <- 18
mock_df_organisms <- 1
mock_df=data.frame(
  Samples=c(
    "MOV2AA","MOV2AB","MOV2AC","MOV2AMA","MOV2AMB","MOV2AMC","MOV2MA","MOV2MB","MOV2MC",
    "MOXTAA","MOXTAB","MOXTAC","MOXTAMA","MOXTAMB","MOXTAMC","MOXTMA","MOXTMB","MOXTMC"),
  SampleID=c(
    "MOV2A","MOV2A","MOV2A","MOV2AM","MOV2AM","MOV2AM","MOV2M","MOV2M","MOV2M",
    "MOXTA","MOXTA","MOXTA","MOXTAM","MOXTAM","MOXTAM","MOXTM","MOXTM","MOXTM"),
  Chemistry = integer(mock_count),
  Probe = integer(mock_count),
  Read_Count = integer(mock_count),
  Total_Read_Length = integer(mock_count),
  Unique_Coloc = integer(mock_count),
  ARG = integer(mock_count),
  MGE = integer(mock_count)
)
non_mock_count <- 54
non_mock_df_organisms <- 3
non_mock_df=data.frame(
  Samples=c(
    "BFV2AA","BFV2AB","BFV2AC","BFV2AMA","BFV2AMB","BFV2AMC","BFV2MA","BFV2MB","BFV2MC",
    "BFXTAA","BFXTAB","BFXTAC","BFXTAMA","BFXTAMB","BFXTAMC","BFXTMA","BFXTMB","BFXTMC",
    "HFV2AA","HFV2AB","HFV2AC","HFV2AMA","HFV2AMB","HFV2AMC","HFV2MA","HFV2MB","HFV2MC",
    "HFXTAA","HFXTAB","HFXTAC","HFXTAMA","HFXTAMB","HFXTAMC","HFXTMA","HFXTMB","HFXTMC",
    "SV2AA","SV2AB","SV2AC","SV2AMA","SV2AMB","SV2AMC","SV2MA","SV2MB","SV2MC",
    "SXTAA","SXTAB","SXTAC","SXTAMA","SXTAMB","SXTAMC","SXTMA","SXTMB","SXTMC"),
  SampleID=c(
    "BFV2A","BFV2A","BFV2A","BFV2AM","BFV2AM","BFV2AM","BFV2M","BFV2M","BFV2M",
    "BFXTA","BFXTA","BFXTA","BFXTAM","BFXTAM","BFXTAM","BFXTM","BFXTM","BFXTM",
    "HFV2A","HFV2A","HFV2A","HFV2AM","HFV2AM","HFV2AM","HFV2M","HFV2M","HFV2M",
    "HFXTA","HFXTA","HFXTA","HFXTAM","HFXTAM","HFXTAM","HFXTM","HFXTM","HFXTM",
    "SV2A","SV2A","SV2A","SV2AM","SV2AM","SV2AM","SV2M","SV2M","SV2M",
    "SXTA","SXTA","SXTA","SXTAM","SXTAM","SXTAM","SXTM","SXTM","SXTM"),
  Chemistry = integer(non_mock_count),
  SampleType = integer(non_mock_count),
  Probe = integer(non_mock_count),
  Read_Count = integer(non_mock_count),
  Total_Read_Length = integer(non_mock_count),
  Unique_Coloc = integer(non_mock_count),
  ARG = integer(non_mock_count),
  MGE = integer(non_mock_count)
)

all_count <- 72
all_df_organisms <- 4
all_df=data.frame(
  Samples=c(
    "BFV2AA","BFV2AB","BFV2AC","BFV2AMA","BFV2AMB","BFV2AMC","BFV2MA","BFV2MB","BFV2MC",
    "BFXTAA","BFXTAB","BFXTAC","BFXTAMA","BFXTAMB","BFXTAMC","BFXTMA","BFXTMB","BFXTMC",
    "HFV2AA","HFV2AB","HFV2AC","HFV2AMA","HFV2AMB","HFV2AMC","HFV2MA","HFV2MB","HFV2MC",
    "HFXTAA","HFXTAB","HFXTAC","HFXTAMA","HFXTAMB","HFXTAMC","HFXTMA","HFXTMB","HFXTMC",
    "MOV2AA","MOV2AB","MOV2AC","MOV2AMA","MOV2AMB","MOV2AMC","MOV2MA","MOV2MB","MOV2MC",
    "MOXTAA","MOXTAB","MOXTAC","MOXTAMA","MOXTAMB","MOXTAMC","MOXTMA","MOXTMB","MOXTMC",
    "SV2AA","SV2AB","SV2AC","SV2AMA","SV2AMB","SV2AMC","SV2MA","SV2MB","SV2MC",
    "SXTAA","SXTAB","SXTAC","SXTAMA","SXTAMB","SXTAMC","SXTMA","SXTMB","SXTMC"),
  SampleID=c(
    "BFV2A","BFV2A","BFV2A","BFV2AM","BFV2AM","BFV2AM","BFV2M","BFV2M","BFV2M",
    "BFXTA","BFXTA","BFXTA","BFXTAM","BFXTAM","BFXTAM","BFXTM","BFXTM","BFXTM",
    "HFV2A","HFV2A","HFV2A","HFV2AM","HFV2AM","HFV2AM","HFV2M","HFV2M","HFV2M",
    "HFXTA","HFXTA","HFXTA","HFXTAM","HFXTAM","HFXTAM","HFXTM","HFXTM","HFXTM",
    "MOV2A","MOV2A","MOV2A","MOV2AM","MOV2AM","MOV2AM","MOV2M","MOV2M","MOV2M",
    "MOXTA","MOXTA","MOXTA","MOXTAM","MOXTAM","MOXTAM","MOXTM","MOXTM","MOXTM",
    "SV2A","SV2A","SV2A","SV2AM","SV2AM","SV2AM","SV2M","SV2M","SV2M",
    "SXTA","SXTA","SXTA","SXTAM","SXTAM","SXTAM","SXTM","SXTM","SXTM"),
  Chemistry = integer(all_count),
  SampleType = integer(all_count),
  Probe = integer(all_count),
  Read_Count = integer(all_count),
  Total_Read_Length = integer(all_count),
  Unique_Coloc = integer(all_count),
  ARG = integer(all_count),
  MGE = integer(all_count)
)
              
non_mock_df <- fill_data_frame(
  non_mock_df, non_mock_count, non_mock_df_organisms, TRUE)

mock_df <- fill_data_frame(
  mock_df, mock_count, mock_df_organisms)

all_df <- fill_data_frame(
  all_df, all_count, all_df_organisms, TRUE)

output_folder <- './output/regression_results/restricted/'

if (file.exists(output_folder) == FALSE)
  dir.create(output_folder)

lmer_function_list <- c(lmer_coloc, lmer_arg, lmer_mge, lmer_read)
plot_function_list <- c(plot_coloc, plot_arg, plot_mge, plot_read)
lmer_filename_list <- c("coloc_", "arg_", "mge_", "read_counts_")
lmer_org_info_list <- c("re_mock_", "re_nonmock_", "re_all_")
df_list <- list(mock_df, non_mock_df, all_df)

summary <- list()
name <- list()

for (i in 1:4) {
  for (j in 1:3) {
    summary <- append(summary, write_lmer(
      lmer_function_list[[i]], lmer_filename_list[i], lmer_org_info_list[j], 
      plot_function_list[[i]], df_list[[j]], j >= 2))
    name <- append(name, paste(lmer_org_info_list[j], lmer_filename_list[i]))
  }
}
df.reml <- data_frame(
  row.names = array(unlist(name)),
  'REML' = array(unlist(summary))
)

write.csv(df.reml, file=paste(output_folder, 're_REML.csv'))