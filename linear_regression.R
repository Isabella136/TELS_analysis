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

library(rjson)
library(rgl)
library(nlme)
library(lme4)
library(RLRsim)
library(MuMIn)

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
              Random2 = integer(sample_count))
              

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
  
  if ((i-1)%%(sample_count/8) < 3) df$Probe[i] <- 'ARG'
  else if ((i-1)%%(sample_count/8) < 6) df$Probe[i] <- 'ARG-MGE'
  else if ((i-1)%%(sample_count/8) < 9) df$Probe[i] <- 'MGE'
  else {
    df$Probe[i] <- 'None'
    df$Random[i] <- df$Chemistry[i]
    if ((i-1)%%12 == 9) df$Random2[i] <- 'ARG'
    else if ((i-1)%%12 == 10) df$Random2[i] <- 'ARG-MGE'
    else if ((i-1)%%12 == 11) df$Random2[i] <- 'MGE'
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

bestfit <- function(dependent, indep, df, deciding_factor, random_var){
  returns_best <- function(first, second) {
    if (deciding_factor == "BIC") {
      bic_df <- BIC(first, second)
      if (bic_df$BIC[1] > bic_df$BIC[2]) return(second)
      return(first)
    }
    else if (deciding_factor == "AIC") {
      aic_df <- AIC(first, second)
      if (aic_df$AIC[1] > aic_df$AIC[2]) return(second)
      return(first)
    }
    else if (deciding_factor == "marginal_r") {
      r.sq1 <- r.squaredGLMM(first)
      r.sq2 <- r.squaredGLMM(second)
      if (r.sq1$R2m[1] < r.sq2$R2m[1]) return(second)
      return(first)
    }
    else if (deciding_factor == "conditional_r") {
      r.sq1 <- r.squaredGLMM(first)
      r.sq2 <- r.squaredGLMM(second)
      if (r.sq1$R2c[1] < r.sq2$R2c[1]) return(second)
      return(first)
    }
  }
  
  find_all_random_models <- function(rhs) {
    toReturn <- list()
    if (grepl(paste(random_var, " )"), as.character(rhs[2])) ) {
      toReturn<-find_all_random_models(as.list(rhs[[2]]))
    }
    if (grepl(paste(random_var, " )"), as.character(rhs[3])) ) {
      toReturn <- c(toReturn,as.character(rhs[[3]][[2]])[2])
    }
    return(toReturn)
  }
  
  recursiveBestFit <- function(base_model, fixed_to_add_list){
    best <- base_model
    best_var <- NA
    for (fixed_to_add in fixed_to_add_list) {
      fixed_added <- update(base_model, .~. + str2lang(fixed_to_add))
      fixed_added_in_effect <- update(fixed_added, .~. + (str2lang(fixed_to_add)|random_var))
      
      fixed_added_with_inter_plus0 <- fixed_added
      fixed_added_with_inter_plus1 <- fixed_added
      fixed_added_with_inter_plus2 <- fixed_added
      fixed_added_in_effect_plus0 <- fixed_added_in_effect
      fixed_added_in_effect_plus1 <- fixed_added_in_effect
      fixed_added_in_effect_plus2 <- fixed_added_in_effect
      fixed_added_with_mult_plus0 <- fixed_added
      fixed_added_with_mult_plus1 <- fixed_added
      fixed_added_with_mult_plus2 <- fixed_added
      
      
      for (fixed_in_effect in find_all_random_models(as.list(base_model@call[["formula"]][[3]]))) {
        if (grepl(":", fixed_in_effect) ||  grepl("*", fixed_in_effect)) {
          for (single_fixed in as.character(as.list(str2lang(fixed_in_effect))[2:3])) {
            
            plus1.fixed_and_interaction_added <- update(fixed_added, .~. + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
            plus1.fixed_and_interaction_added_in_effect <- update(fixed_added_in_effect, .~. + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
            plus1.fixed_and_mult_added <- update(fixed_added, .~. + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))
            
            plus1.updated_fixed_and_interaction_added <- update(fixed_added_with_inter_plus0, .~. + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
            plus1.updated_fixed_and_interaction_added_in_effect <- update(fixed_added_in_effect_plus0, .~. + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
            plus1.updated_fixed_and_mult_added <- update(fixed_added_with_mult_plus0, .~. + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))
            
            plus2.fixed_and_interaction_added <- update(fixed_added_with_inter_plus1, .~. + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))                                     
            plus2.fixed_and_interaction_added_in_effect <- update(fixed_added_in_effect_plus1, .~. + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))                                     
            plus2.fixed_and_mult_added <- update(fixed_added_with_mult_plus1, .~. + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))                                     
            
            
            plus0.fixed_and_interaction_added_with_replacement <- update(fixed_added, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
            plus0.fixed_and_interaction_added_in_effect_with_replacement <- update(fixed_added_in_effect, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
            plus0.fixed_and_mult_added_with_replacement <- update(fixed_added, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))
            
            plus0.updated_fixed_and_interaction_added_with_replacement <- update(fixed_added_with_inter_plus0, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
            plus0.updated_fixed_and_interaction_added_in_effect_with_replacement <- update(fixed_added_in_effect_plus0, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
            plus0.updated_fixed_and_mult_added_with_replacement <- update(fixed_added_with_mult_plus0, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))
            
            
            plus1.fixed_and_interaction_added_with_replacement <- update(fixed_added_with_inter_plus1, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
            plus1.fixed_and_interaction_added_in_effect_with_replacement <- update(fixed_added_in_effect_plus1, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
            plus1.fixed_and_mult_added_with_replacement <- update(fixed_added_with_mult_plus1, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))
            
            plus2.fixed_and_interaction_added_with_replacement <- update(fixed_added_with_inter_plus2, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
            plus2.fixed_and_interaction_added_in_effect_with_replacement <- update(fixed_added_in_effect_plus2, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
            plus2.fixed_and_mult_added_with_replacement <- update(fixed_added_with_mult_plus2, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))
            
            fixed_added_with_inter_plus0 <- returns_best(fixed_added_with_inter_plus0, plus0.fixed_and_interaction_added_with_replacement)
            fixed_added_with_inter_plus0 <- returns_best(fixed_added_with_inter_plus0, plus0.updated_fixed_and_interaction_added_with_replacement)
            fixed_added_with_inter_plus1 <- returns_best(fixed_added_with_inter_plus1, plus1.fixed_and_interaction_added)
            fixed_added_with_inter_plus1 <- returns_best(fixed_added_with_inter_plus1, plus1.updated_fixed_and_interaction_added)
            fixed_added_with_inter_plus1 <- returns_best(fixed_added_with_inter_plus1, plus1.fixed_and_interaction_added_with_replacement)
            fixed_added_with_inter_plus2 <- returns_best(fixed_added_with_inter_plus2, plus2.fixed_and_interaction_added)
            fixed_added_with_inter_plus2 <- returns_best(fixed_added_with_inter_plus2, plus2.fixed_and_interaction_added_with_replacement)
            
            fixed_added_in_effect_plus0 <- returns_best(fixed_added_in_effect_plus0, plus0.fixed_and_interaction_added_in_effect_with_replacement)
            fixed_added_in_effect_plus0 <- returns_best(fixed_added_in_effect_plus0, plus0.updated_fixed_and_interaction_added_in_effect_with_replacement)
            fixed_added_in_effect_plus1 <- returns_best(fixed_added_in_effect_plus1, plus1.fixed_and_interaction_added_in_effect)
            fixed_added_in_effect_plus1 <- returns_best(fixed_added_in_effect_plus1, plus1.updated_fixed_and_interaction_added_in_effect)
            fixed_added_in_effect_plus1 <- returns_best(fixed_added_in_effect_plus1, plus1.fixed_and_interaction_added_in_effect_with_replacement)
            fixed_added_in_effect_plus2 <- returns_best(fixed_added_in_effect_plus2, plus2.fixed_and_interaction_added_in_effect)
            fixed_added_in_effect_plus2 <- returns_best(fixed_added_in_effect_plus2, plus2.fixed_and_interaction_added_in_effect_with_replacement)
            
            fixed_added_with_mult_plus0 <- returns_best(fixed_added_with_mult_plus0, plus0.fixed_and_mult_added_with_replacement)
            fixed_added_with_mult_plus0 <- returns_best(fixed_added_with_mult_plus0, plus0.updated_fixed_and_mult_added_with_replacement)
            fixed_added_with_mult_plus1 <- returns_best(fixed_added_with_mult_plus1, plus1.fixed_and_mult_added)
            fixed_added_with_mult_plus1 <- returns_best(fixed_added_with_mult_plus1, plus1.updated_fixed_and_mult_added)
            fixed_added_with_mult_plus1 <- returns_best(fixed_added_with_mult_plus1, plus1.fixed_and_mult_added_with_replacement)
            fixed_added_with_mult_plus2 <- returns_best(fixed_added_with_mult_plus2, plus2.fixed_and_mult_added)
            fixed_added_with_mult_plus2 <- returns_best(fixed_added_with_mult_plus2, plus2.fixed_and_mult_added_with_replacement)
            
          }
            
        }
        
        else {
          plus1.fixed_and_interaction_added <- update(fixed_added, .~. + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
          plus1.fixed_and_interaction_added_in_effect <- update(fixed_added_in_effect, .~. + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
          plus1.fixed_and_mult_added <- update(fixed_added, .~. + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))
          
          plus1.updated_fixed_and_interaction_added <- update(fixed_added_with_inter_plus0, .~. + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
          plus1.updated_fixed_and_interaction_added_in_effect <- update(fixed_added_in_effect_plus0, .~. + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
          plus1.updated_fixed_and_mult_added <- update(fixed_added_with_mult_plus0, .~. + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))
          
          plus2.fixed_and_interaction_added <- update(fixed_added_with_inter_plus1, .~. + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))                                     
          plus2.fixed_and_interaction_added_in_effect <- update(fixed_added_in_effect_plus1, .~. + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))                                     
          plus2.fixed_and_mult_added <- update(fixed_added_with_mult_plus1, .~. + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))                                     
          
          
          plus0.fixed_and_interaction_added_with_replacement <- update(fixed_added, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
          plus0.fixed_and_interaction_added_in_effect_with_replacement <- update(fixed_added_in_effect, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
          plus0.fixed_and_mult_added_with_replacement <- update(fixed_added, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))
          
          plus0.updated_fixed_and_interaction_added_with_replacement <- update(fixed_added_with_inter_plus0, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
          plus0.updated_fixed_and_interaction_added_in_effect_with_replacement <- update(fixed_added_in_effect_plus0, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
          plus0.updated_fixed_and_mult_added_with_replacement <- update(fixed_added_with_mult_plus0, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))
          
          
          plus1.fixed_and_interaction_added_with_replacement <- update(fixed_added_with_inter_plus1, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
          plus1.fixed_and_interaction_added_in_effect_with_replacement <- update(fixed_added_in_effect_plus1, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
          plus1.fixed_and_mult_added_with_replacement <- update(fixed_added_with_mult_plus1, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))
          
          plus2.fixed_and_interaction_added_with_replacement <- update(fixed_added_with_inter_plus2, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
          plus2.fixed_and_interaction_added_in_effect_with_replacement <- update(fixed_added_in_effect_plus2, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add):str2lang(single_fixed)|random_var))
          plus2.fixed_and_mult_added_with_replacement <- update(fixed_added_with_mult_plus2, .~. - str2lang(fixed_in_effect) + (str2lang(fixed_to_add)*str2lang(single_fixed)|random_var))
          
          fixed_added_with_inter_plus0 <- returns_best(fixed_added_with_inter_plus0, plus0.fixed_and_interaction_added_with_replacement)
          fixed_added_with_inter_plus0 <- returns_best(fixed_added_with_inter_plus0, plus0.updated_fixed_and_interaction_added_with_replacement)
          fixed_added_with_inter_plus1 <- returns_best(fixed_added_with_inter_plus1, plus1.fixed_and_interaction_added)
          fixed_added_with_inter_plus1 <- returns_best(fixed_added_with_inter_plus1, plus1.updated_fixed_and_interaction_added)
          fixed_added_with_inter_plus1 <- returns_best(fixed_added_with_inter_plus1, plus1.fixed_and_interaction_added_with_replacement)
          fixed_added_with_inter_plus2 <- returns_best(fixed_added_with_inter_plus2, plus2.fixed_and_interaction_added)
          fixed_added_with_inter_plus2 <- returns_best(fixed_added_with_inter_plus2, plus2.fixed_and_interaction_added_with_replacement)
          
          fixed_added_in_effect_plus0 <- returns_best(fixed_added_in_effect_plus0, plus0.fixed_and_interaction_added_in_effect_with_replacement)
          fixed_added_in_effect_plus0 <- returns_best(fixed_added_in_effect_plus0, plus0.updated_fixed_and_interaction_added_in_effect_with_replacement)
          fixed_added_in_effect_plus1 <- returns_best(fixed_added_in_effect_plus1, plus1.fixed_and_interaction_added_in_effect)
          fixed_added_in_effect_plus1 <- returns_best(fixed_added_in_effect_plus1, plus1.updated_fixed_and_interaction_added_in_effect)
          fixed_added_in_effect_plus1 <- returns_best(fixed_added_in_effect_plus1, plus1.fixed_and_interaction_added_in_effect_with_replacement)
          fixed_added_in_effect_plus2 <- returns_best(fixed_added_in_effect_plus2, plus2.fixed_and_interaction_added_in_effect)
          fixed_added_in_effect_plus2 <- returns_best(fixed_added_in_effect_plus2, plus2.fixed_and_interaction_added_in_effect_with_replacement)
          
          fixed_added_with_mult_plus0 <- returns_best(fixed_added_with_mult_plus0, plus0.fixed_and_mult_added_with_replacement)
          fixed_added_with_mult_plus0 <- returns_best(fixed_added_with_mult_plus0, plus0.updated_fixed_and_mult_added_with_replacement)
          fixed_added_with_mult_plus1 <- returns_best(fixed_added_with_mult_plus1, plus1.fixed_and_mult_added)
          fixed_added_with_mult_plus1 <- returns_best(fixed_added_with_mult_plus1, plus1.updated_fixed_and_mult_added)
          fixed_added_with_mult_plus1 <- returns_best(fixed_added_with_mult_plus1, plus1.fixed_and_mult_added_with_replacement)
          fixed_added_with_mult_plus2 <- returns_best(fixed_added_with_mult_plus2, plus2.fixed_and_mult_added)
          fixed_added_with_mult_plus2 <- returns_best(fixed_added_with_mult_plus2, plus2.fixed_and_mult_added_with_replacement)
          
        }
        
      }
    
      temp_best <- best
      best <- returns_best(best, fixed_added_with_inter_plus0)
      best <- returns_best(best, fixed_added_with_inter_plus1)
      best <- returns_best(best, fixed_added_with_inter_plus2)
      best <- returns_best(best, fixed_added_in_effect_plus0)
      best <- returns_best(best, fixed_added_in_effect_plus1)
      best <- returns_best(best, fixed_added_in_effect_plus2)
      best <- returns_best(best, fixed_added_with_mult_plus0)
      best <- returns_best(best, fixed_added_with_mult_plus1)
      best <- returns_best(best, fixed_added_with_mult_plus2)
      
      if (best != temp_best) best_var <- fixed_to_add
      
    }
    if (best == base_model)
      return(best)
  
  }
  
  
  

}




lmer.coloc_grand_fit <- lmer(Unique_Coloc ~ Probe 
                             * SampleType 
                             * Chemistry 
                             * Total_Read_Length
                             + (SampleType:Probe | Random),
                             data=df)
r.squaredGLMM(lmer.coloc_grand_fit)
AIC(lmer.coloc_grand_fit)

lmer.coloc_grand_fit2 <- update(lmer.coloc_grand_fit, .~.-(SampleType:Probe | Random) + (SampleType*Probe | Random))
lmer.coloc_grand_fit2 <- update(lmer.coloc_grand_fit, .~. + (Probe | Random))

lmer.coloc_grand_fit2$call

test <- r.squaredGLMM(lmer.coloc_grand_fit2)
AIC(lmer.coloc_grand_fit,lmer.coloc_grand_fit2)
test <- BIC(lmer.coloc_grand_fit,lmer.coloc_grand_fit2)

summary(lmer.coloc_grand_fit)
residuals(lmer.coloc_grand_fit)[abs(residuals(lmer.coloc_grand_fit)) >= .5]
ranef(lmer.coloc_grand_fit)
plot(x = df$Unique_Coloc, 
     y = round(fitted(lmer.coloc_grand_fit)),
     main = "Unique Colocalization",
     xlab = "Original Values",
     ylab = "Fitted Values")

lmer.no_probe <- lmer(Unique_Coloc ~ SampleType
                      * Chemistry 
                      * Total_Read_Length 
                      + (SampleType | Random),
                      data=df)
r.squaredGLMM(lmer.no_probe)

summary(lmer.no_probe)
residuals(lmer.no_probe)[abs(residuals(lmer.no_probe)) >= 5]
ranef(lmer.no_probe)

lmer.no_chemistry <- lmer(Unique_Coloc ~ Probe
                          * SampleType 
                          * Total_Read_Length 
                          + (SampleType:Probe | Random),
                          data=df)
r.squaredGLMM(lmer.no_chemistry)

summary(lmer.no_chemistry)
residuals(lmer.no_chemistry)[abs(residuals(lmer.no_chemistry)) >= 5]
ranef(lmer.no_chemistry)

lmer.no_length <- lmer(Unique_Coloc ~ Probe 
                       * SampleType 
                       * Chemistry 
                       + (SampleType:Probe | Random),
                       data=df)
r.squaredGLMM(lmer.no_length)

summary(lmer.no_length)
residuals(lmer.no_length)[abs(residuals(lmer.no_length)) >= 5]
ranef(lmer.no_length)







lmer.read_grand_fit <- lmer(Read_Count ~ Probe 
                            * SampleType 
                            * Chemistry 
                            * Total_Read_Length
                            + (Probe:SampleType | Random),
                            data=df)
r.squaredGLMM(lmer.read_grand_fit)

summary(lmer.read_grand_fit)
residuals(lmer.read_grand_fit)[abs(residuals(lmer.read_grand_fit)) >= 5000]
ranef(lmer.read_grand_fit)
plot(x = df$Read_Count, 
     y = round(fitted(lmer.read_grand_fit)),
     main = "Unique Colocalization",
     xlab = "Original Values",
     ylab = "Fitted Values")


