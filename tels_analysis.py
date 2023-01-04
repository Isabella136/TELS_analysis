import configparser
#from tels_analysis.colocalization_analyzer import colocalization_analyzer
from tels_analysis.heatmap_analyzer import heatmap_analyzer
from tels_analysis.statistical_analyzer import statistical_analyzer
from tels_analysis.compos_richness_analyzer import compos_richness_analyzer
from tels_analysis.abundance_analyzer import abundance_analyzer
from tels_analysis.stacked_abundance_analyzer import stacked_abundance_analyzer
from tels_analysis.venn_analyzer import venn_analyzer
from tels_analysis.special import special
from tels_analysis.special_coloc import special_coloc
import sys, getopt, gzip, shutil, os

fileList = ["BFV2AA",		"BFV2AB",		"BFV2AC",		"BFV2AMA",		"BFV2AMB",		"BFV2AMC",
			"BFV2MA",		"BFV2MB",		"BFV2MC",		"BFV2NEGA",		"BFV2NEGAM",	"BFV2NEGM",		
			"BFXTAA",		"BFXTAB",		"BFXTAC",		"BFXTAMA",		"BFXTAMB",		"BFXTAMC",
			"BFXTMA",		"BFXTMB",		"BFXTMC",		"BFXTNEGA",		"BFXTNEGAM",	"BFXTNEGM", 
			"HFV2AA",		"HFV2AB",		"HFV2AC",		"HFV2AMA",		"HFV2AMB",		"HFV2AMC", 
			"HFV2MA",		"HFV2MB",		"HFV2MC",		"HFV2NEGA",		"HFV2NEGAM",	"HFV2NEGM",
			"HFXTAA",		"HFXTAB",		"HFXTAC",		"HFXTAMA",		"HFXTAMB",		"HFXTAMC", 
			"HFXTMA",		"HFXTMB",		"HFXTMC",		"HFXTNEGA",		"HFXTNEGAM",	"HFXTNEGM", 
			"MOV2AA",		"MOV2AB",		"MOV2AC",		"MOV2AMA",		"MOV2AMB",		"MOV2AMC",
			"MOV2MA",		"MOV2MB",		"MOV2MC",		"MOV2NEGA",		"MOV2NEGAM",	"MOV2NEGM",
			"MOXTAA",		"MOXTAB",		"MOXTAC",		"MOXTAMA",		"MOXTAMB",		"MOXTAMC",
			"MOXTMA",		"MOXTMB",		"MOXTMC",		"MOXTNEGA",		"MOXTNEGAM",	"MOXTNEGM",
			"SV2AA",		"SV2AB",		"SV2AC",		"SV2AMA",		"SV2AMB",		"SV2AMC",
			"SV2MA",		"SV2MB", 		"SV2MC",		"SV2NEGA",		"SV2NEGAM",		"SV2NEGM",
			"SXTAA",		"SXTAB", 		"SXTAC",		"SXTAMA",		"SXTAMB",		"SXTAMC",
			"SXTMA",		"SXTMB",		"SXTMC",		"SXTNEGA",		"SXTNEGAM",		"SXTNEGM"]

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

if not os.path.exists(outputFolder + "/"):
	os.makedirs(outputFolder + "/")

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

if config.getboolean("STEPS", "SPECIAL"):
	mgeInfo = special(config.get("SOURCE_FILE", "SOURCE_PREFIX"),
									"AlignedToMegares.csv",
									config.get("SOURCE_FILE", "SOURCE_SUFFIX"),
									config.get("SOURCE_EXTENSION", "SHORT_MGE"),
									config.get("SOURCE_FILE", "MGE_CLASSIFICATION"))
	colocInfo = special_coloc(config.get("SOURCE_FILE", "SOURCE_PREFIX"),
									"AlignedToMegares.csv",
									config.get("SOURCE_FILE", "SOURCE_SUFFIX"),
									config.get("SOURCE_EXTENSION", "COLOCALIZATIONS_RICHNESS"))
	for fileName in fileList:
		mgeInfo.addToMobilomeInfo(fileName)
		colocInfo.addColocInfo(fileName)
	mgeInfo.writeMobilomeInfo(outputFolder + "/mobilome_info.csv")
	mgeInfo.addComparisonInfo()
	mgeInfo.writeComparisonMobilomeInfo(outputFolder + "/comparison_mobilome_info.csv")
	colocInfo.writeColocInfo(outputFolder + "/coloc_info.csv")

if config.getboolean("STEPS", "STATS"):
	for fileName in fileList:
		statsFileOutput = config.get("SOURCE_FILE", "SOURCE_PREFIX") + fileName + config.get("SOURCE_FILE", "SOURCE_SUFFIX") + config.get("OUTPUT_EXTENSION", "STATS")
		if not(os.path.exists(statsFileOutput)):
			readListSource = config.get("SOURCE_FILE", "SOURCE_PREFIX") + fileName + config.get("SOURCE_FILE", "SOURCE_SUFFIX")
			inputFile = open(readListSource, "r")
			lineNum = 0
			for line in inputFile:
				lineNum += 1
			inputFile.close()
			readListSource = config.get("SOURCE_FILE", "SOURCE_PREFIX") + fileName + ".ccs.fastq"
			dupFile = open(readListSource, "r")
			dupLineNum = 0
			for line in dupFile:
				dupLineNum += 1
			dupFile.close()
			outputFile = open(statsFileOutput, "w")
			outputFile.write("DUPLICATED_STATS_NUM_OF_READS," + str(int(dupLineNum/4)))
			outputFile.write("\n")
			outputFile.write("DEDUPLICATED_STATS_NUM_OF_READS," + str(int(lineNum/4)))
			outputFile.close()

if config.getboolean("STEPS", "STATISTICAL_ANALYSIS"):
	statisticalAnalyzer = statistical_analyzer(config.get("SOURCE_FILE", "SOURCE_PREFIX"), 
													   config.get("SOURCE_FILE", "SOURCE_SUFFIX"), 
													   config.get("SOURCE_EXTENSION", "SHORT_AMR_DIV"), 
													   config.get("SOURCE_EXTENSION", "SHORT_MGE"),
													   config.get("OUTPUT_EXTENSION", "STATS"))
	for fileName in fileList:
		statisticalAnalyzer.analyzeFile(fileName)
	statisticalAnalyzer.printAnalysis(outputFolder, config.get("OUTPUT_EXTENSION", "STATISTICAL_ANALYSIS"))

if config.getboolean("STEPS", "COMPOS_RICHNESS_ANALYSIS"):
	crAnalyzer = compos_richness_analyzer(config.get("SOURCE_FILE", "SOURCE_PREFIX"), 
													   config.get("SOURCE_FILE", "SOURCE_SUFFIX"), 
													   config.get("SOURCE_EXTENSION", "SHORT_AMR_DIV"), 
													   config.get("SOURCE_EXTENSION", "SHORT_MGE"),
													   config.get("SOURCE_FILE", "MGE_CLASSIFICATION"))
	for fileName in fileList:
		crAnalyzer.analyzeFile(fileName, outputFolder + "/ARG_composition", outputFolder + "/MGE_composition", config.get("OUTPUT_EXTENSION", "INDIV_COMPOS_CHART"))
	crAnalyzer.printAnalysis(outputFolder, config.get("OUTPUT_EXTENSION", "COMPOS_RICHNESS_ANALYSIS"))

if config.getboolean("STEPS", "HEATMAP"):
	heatmapAnalyzer = heatmap_analyzer(config.get("SOURCE_FILE", "SOURCE_PREFIX"), 
													   config.get("SOURCE_FILE", "SOURCE_SUFFIX"), 
													   config.get("SOURCE_EXTENSION", "SHORT_AMR_DIV"), 
													   config.get("SOURCE_EXTENSION", "SHORT_MGE"),
													   config.get("SOURCE_FILE", "MEGARES"))
	for fileName in fileList:
		heatmapAnalyzer.addToMaps(fileName)
	heatmapAnalyzer.makeMaps(outputFolder + "/ARG_heatmap", config.get("OUTPUT_EXTENSION", "HEATMAP"))

if config.getboolean("STEPS", "FILE_SIZE"):
	fileOfSizes = open(outputFolder + "/" + config.get("OUTPUT_EXTENSION", "FILE_SIZE"), "w")
	for fileName in fileList:
		with gzip.open(config.get("SOURCE_FILE", "INITIAL_SOURCE_PREFIX") + fileName + config.get("SOURCE_FILE", "SOURCE_SUFFIX"), "rb") as input:
			with open ("temp_files/deduplicated_sequel-demultiplex.temp.ccs.fastq", "wb") as output:
				shutil.copyfileobj(input, output)
		fileOfSizes.write(fileName + "," + str(os.stat("temp_files/deduplicated_sequel-demultiplex.temp.ccs.fastq").st_size) + ",\n")
	fileOfSizes.close()

if config.getboolean("STEPS", "VIOLIN"):
	abundanceAnalyzer = abundance_analyzer(config.get("SOURCE_FILE", "SOURCE_PREFIX"), 
													   config.get("SOURCE_FILE", "SOURCE_SUFFIX"), 
													   outputFolder + "/" + config.get("OUTPUT_EXTENSION", "FILE_SIZE"),
													   config.get("SOURCE_EXTENSION", "SHORT_AMR_DIV"), 
													   config.get("SOURCE_EXTENSION", "SHORT_MGE"),
													   config.get("SOURCE_FILE", "MEGARES_FASTA"))
	for fileName in fileList:
		abundanceAnalyzer.findAbsoluteAbundance(fileName)
	abundanceAnalyzer.makeViolinPlot(outputFolder, config.get("OUTPUT_EXTENSION", "VIOLIN"))

if config.getboolean("STEPS", "STACKED"):
	stackedAnalyzer = stacked_abundance_analyzer(config.get("SOURCE_FILE", "SOURCE_PREFIX"), 
													   config.get("SOURCE_FILE", "SOURCE_SUFFIX"),
													   config.get("SOURCE_EXTENSION", "SHORT_AMR_DIV"), 
													   config.get("SOURCE_EXTENSION", "SHORT_MGE"),
													   config.get("OUTPUT_EXTENSION", "STATS"))
	for fileName in fileList:
		stackedAnalyzer.findAbsoluteAbundance(fileName)
	stackedAnalyzer.makeStack(outputFolder, config.get("OUTPUT_EXTENSION", "STACKED"))

if config.getboolean("STEPS", "VENN"):
	vennAnalyzer = venn_analyzer(config.get("SOURCE_FILE", "SOURCE_PREFIX"), 
													   config.get("SOURCE_FILE", "SOURCE_SUFFIX"),
													   config.get("SOURCE_EXTENSION", "SHORT_AMR_DIV"), 
													   config.get("SOURCE_EXTENSION", "SHORT_MGE"))
	for fileName in fileList:
		vennAnalyzer.addToCount(fileName)
	vennAnalyzer.makeVenn(outputFolder, config.get("OUTPUT_EXTENSION", "VENN"))
