import configparser
from tels_analysis.colocalization_analyzer import colocalization_analyzer
from tels_analysis.statistical_analyzer import statistical_analyzer
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

outputFolder = ""
configFile = ""
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
statisticalAnalyzer = statistical_analyzer(config.get("SOURCE_FILE", "SOURCE_PREFIX"), 
												   config.get("SOURCE_FILE", "SOURCE_SUFFIX"), 
												   config.get("SOURCE_EXTENSION", "SHORT_AMR_DIV"), 
												   config.get("SOURCE_EXTENSION", "SHORT_MGE"),
												   config.get("SOURCE_EXTENSION", "STATS"))
for fileName in fileList:
	colocalizationAnalyzer = colocalization_analyzer(fileName, config.get("SOURCE_FILE", "SOURCE_PREFIX"), 
												   config.get("SOURCE_FILE", "SOURCE_SUFFIX"), 
												   config.get("SOURCE_EXTENSION", "COLOCALIZATIONS"), 
												   config.get("SOURCE_EXTENSION", "READS_LENGTH"))
	colocalizationAnalyzer.makeChart( fileName, outputFolder + "/" + 
								  config.get("OUTPUT_FILE", "OUTPUT_PREFIX"), 
								  config.get("OUTPUT_FILE", "OUTPUT_SUFFIX"), 
								  config.get("OUTPUT_EXTENSION", "COLOCALIZATION_ANALYSIS"))
	statisticalAnalyzer.analyzeFile(fileName)

statisticalAnalyzer.printAnalysis(outputFolder, config.get("OUTPUT_EXTENSION", "STATISTICAL_ANALYSIS"))
