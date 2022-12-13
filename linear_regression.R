df=data.frame(Samples=c("BFV2AA","BFV2AB","BFV2AC","BFV2AMA","BFV2AMB","BFV2AMC","BFV2MA","BFV2MB","BFV2MC","BFV2NEGA","BFV2NEGAM","BFV2NEGM",
                        "BFXTAA","BFXTAB","BFXTAC","BFXTAMA","BFXTAMB","BFXTAMC","BFXTMA","BFXTMB","BFXTMC",'BFXTNEGA',"BFXTNEGAM","BFXTNEGM",
                        "HFV2AA", "HFV2AB", "HFV2AC","HFV2AMA","HFV2AMB","HFV2AMC","HFV2MA","HFV2MB","HFV2MC","HFV2NEGA","HFV2NEGAM","HFV2NEGM",
                        "HFXTAA", "HFXTAB", "HFXTAC","HFXTAMA","HFXTAMB","HFXTAMC","HFXTMA","HFXTMB","HFXTMC","HFXTNEGA","HFXTNEGAM","HFXTNEGM",
                        "MOV2AA", "MOV2AB", "MOV2AC","MOV2AMA","MOV2AMB","MOV2AMC","MOV2MA","MOV2MB","MOV2MC","MOV2NEGA","MOV2NEGAM","MOV2NEGM",
                        "MOXTAA", "MOXTAB", "MOXTAC","MOXTAMA","MOXTAMB","MOXTAMC","MOXTMA","MOXTMB","MOXTMC","MOXTNEGA",'MOXTNEGAM',"MOXTNEGM",
                        "SV2AA", "SV2AB", "SV2AC","SV2AMA","SV2AMB","SV2AMC","SV2MA","SV2MB","SV2MC","SV2NEGA","SV2NEGAM","SV2NEGM",
                        "SXTAA", "SXTAB", "SXTAC","SXTAMA","SXTAMB","SXTAMC","SXTMA","SXTMB","SXTMC","SXTNEGA","SXTNEGAM","SXTNEGM"))
df <- cbind(df, Unique_Coloc = integer(96))
for (i in 1:96) {
  colocalizations_richness <- read.csv(paste("~/Documents/GitHub/TELS_analysis/TELS_output/sequel-demultiplex.", ".ccs.fastq_deduplicated.fastq_colocalizations_richness.csv", sep=df$Samples[i]), header=FALSE, row.names=NULL, comment.char="#")
  df$Unique_Coloc[i] = colocalizations_richness$V2[1]
}
