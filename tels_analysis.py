import configparser
#from tels_analysis.sam_analyzer import sam_analyzer
#from tels_analysis.colocalization_analyzer import colocalization_analyzer
from tels_analysis.heatmap_analyzer import heatmap_analyzer
from tels_analysis.statistical_analyzer import statistical_analyzer
from tels_analysis.compos_richness_analyzer import compos_richness_analyzer
from tels_analysis.abundance_analyzer import abundance_analyzer
import sys, getopt

fileList = ["BFV2AA",		"BFV2AB",		"BFV2AC",		"BFV2AMA",		"BFV2AMB",		"BFV2AMC",
			"BFV2MA",		"BFV2MB",		"BFV2MC",		"BFV2NEGA",		"BFV2NEGAM",	"BFV2NEGM",		
			"BFXTAA",		"BFXTAB",		"BFXTAC",		"BFXTAMA",		"BFXTAMB",		"BFXTAMC",
			"BFXTMA",		"BFXTMB",		"BFXTMC",		"BFXTNEGA",		"BFXTNEGAM",	"BFXTNEGM", 
			"HFV2AA",		"HFV2AB",		"HFV2AC",		"HFV2AMA",		"HFV2AMB",		"HFV2AMC", 
			"HFV2MA",		"HFV2MB",		"HFV2MC",		"HFV2NEGA",		"HFV2NEGAM",	"HFV2NEGM",
			"HFXTAA",		"HFXTAB",		"HFXTAC",		"HFXTAMA",		"HFXTAMB",		"HFXTAMC", 
			"HFXTMA",		"HFXTMB",		"HFXTMC",		"HFXTNEGA",		"HFXTNEGAM",	"HFXTNEGM", 
			"MOV2AA",		"MOV2AB",		"MOV2AC",		"MOV2AMB",		"MOV2MA",		"MOV2MB",
			"MOV2MC",		"MOV2NEGA",		"MOV2NEGAM",	"MOV2NEGM",		"MOXTAA",		"MOXTAB", 
			"MOXTAC",		"MOXTAMA",		"MOXTAMB",		"MOXTAMC",		"MOXTMA",		"MOXTMB",
			"MOXTMC",		"MOXTNEGA",		"MOXTNEGAM",	"MOXTNEGM",		"SV2AA",		"SV2AB",
			"SV2AC",		"SV2AMA",		"SV2AMB",		"SV2AMC",		"SV2MA",		"SV2MB", 
			"SV2MC",		"SV2NEGA",		"SV2NEGAM",		"SV2NEGM",		"SXTAA",		"SXTAB", 
			"SXTAC",		"SXTAMA",		"SXTAMB",		"SXTAMC",		"SXTMA",		"SXTMB",
			"SXTMC",		"SXTNEGA",		"SXTNEGAM",		"SXTNEGM"]

outputFolder = "output"
configFile = "config.ini"
try:
    options, args = getopt.getopt(sys.argv[1:], "hc:o:")
except getopt.GetoptError:
    print("tels_analysis.py -c <configFile> -o <outputFolder>")
    sys.exit(-1)
for opt, arg in options:
    if opt == "-h":
        print("List of arguments:\n\n\n-c: config file\n-h: help\n-o: output folder\n")
        sys.exit()
    elif opt == "-c":
        configFile = arg
    elif opt == "-o":
        outputFolder = arg

config = configparser.ConfigParser()
config.read(configFile)

#if config.getboolean("STEPS", "SAM_ANALYSIS"):
#	for fileName in fileList:
#		samAnalyzer = sam_analyzer(fileName, config.get("SOURCE_FILE", "SOURCE_PREFIX"), 
#							 config.get("SOURCE_FILE", "SOURCE_SUFFIX"), 
#							 config.get("SOURCE_EXTENSION", "A_TO_MEGARES"), 
#							 config.get("SOURCE_EXTENSION", "A_TO_MGES"), 
#							 config.get("SOURCE_EXTENSION", "A_TO_KEGG"))
#		samAnalyzer.all_genes_list(outputFolder, config.get("OUTPUT_FILE", "OUTPUT_PREFIX"),
#							 config.get("OUTPUT_FILE", "OUTPUT_SUFFIX"),
#							 config.get("OUTPUT_EXTENSION", "SAM_ANALYSIS"))

#if config.getboolean("STEPS", "COLOCALIZATION_ANALYSIS"):
#	for fileName in fileList:
#		colocalizationAnalyzer = colocalization_analyzer(fileName, config.get("SOURCE_FILE", "SOURCE_PREFIX"), 
#													   config.get("SOURCE_FILE", "SOURCE_SUFFIX"), 
#													   config.get("SOURCE_EXTENSION", "COLOCALIZATIONS"), 
#													   config.get("SOURCE_EXTENSION", "READS_LENGTH"))
#		colocalizationAnalyzer.makeChart( fileName, outputFolder + "/" + 
#									  config.get("OUTPUT_FILE", "OUTPUT_PREFIX"), 
#									  config.get("OUTPUT_FILE", "OUTPUT_SUFFIX"), 
#									  config.get("OUTPUT_EXTENSION", "COLOCALIZATION_ANALYSIS"))

if config.getboolean("STEPS", "STATISTICAL_ANALYSIS"):
	statisticalAnalyzer = statistical_analyzer(config.get("SOURCE_FILE", "SOURCE_PREFIX"), 
													   config.get("SOURCE_FILE", "SOURCE_SUFFIX"), 
													   config.get("SOURCE_EXTENSION", "ARG_SAM_ANALYSIS"), 
													   config.get("SOURCE_EXTENSION", "MGE_SAM_ANALYSIS"),
													   config.get("SOURCE_EXTENSION", "STATS"))
	for fileName in fileList:
		statisticalAnalyzer.analyzeFile(fileName)
	statisticalAnalyzer.printAnalysis(outputFolder, config.get("OUTPUT_EXTENSION", "STATISTICAL_ANALYSIS"))

if config.getboolean("STEPS", "COMPOS_RICHNESS_ANALYSIS"):
	crAnalyzer = compos_richness_analyzer(config.get("SOURCE_FILE", "SOURCE_PREFIX"), 
													   config.get("SOURCE_FILE", "SOURCE_SUFFIX"), 
													   config.get("SOURCE_EXTENSION", "ARG_SAM_ANALYSIS"), 
													   config.get("SOURCE_EXTENSION", "MGE_SAM_ANALYSIS"))
	for fileName in fileList:
		crAnalyzer.analyzeFile(fileName, outputFolder + "/ARG_composition", config.get("OUTPUT_EXTENSION", "INDIV_COMPOS_CHART"))
	crAnalyzer.printAnalysis(outputFolder, config.get("OUTPUT_EXTENSION", "COMPOS_RICHNESS_ANALYSIS"))
if config.getboolean("STEPS", "HEATMAP"):
	heatmapAnalyzer = heatmap_analyzer(config.get("SOURCE_FILE", "SOURCE_PREFIX"), 
													   config.get("SOURCE_FILE", "SOURCE_SUFFIX"), 
													   config.get("SOURCE_EXTENSION", "ARG_SAM_ANALYSIS"), 
													   config.get("SOURCE_EXTENSION", "MGE_SAM_ANALYSIS"),
													   config.get("SOURCE_FILE", "MEGARES"))
	for fileName in fileList:
		heatmapAnalyzer.addToMaps(fileName)
	heatmapAnalyzer.makeMaps(outputFolder + "/ARG_heatmap", config.get("OUTPUT_EXTENSION", "HEATMAP"))
if config.getboolean("STEPS", "VIOLIN"):
	abundanceAnalyzer = abundance_analyzer(config.get("SOURCE_FILE", "ARG_SAM_ANALYSIS_SOURCE_PREFIX"), 
													   config.get("SOURCE_FILE", "SOURCE_SUFFIX"), 
													   config.get("SOURCE_FILE", "INITIAL_SOURCE_PREFIX"),
													   config.get("SOURCE_EXTENSION", "ARG_SAM_ANALYSIS"), 
													   config.get("SOURCE_EXTENSION", "MGE_SAM_ANALYSIS"),
													   config.get("SOURCE_FILE", "MEGARES"))
	for fileName in fileList:
		abundanceAnalyzer.findAbsoluteAbundance(fileName)