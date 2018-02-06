# Before running the script load the modules present in MODULES file
# This is modified from the 'script for combined MRNA and MiRNA analysis'
# Needs gene counts matrix files for MiRNA RNA-seq samples  
# and samples information (Samples.txt) file,
# in the current directory or specified by their complete path

# Start fresh
rm(list=ls())

# Working directory
# Set the working directory, which has all function files
WD <- "/home/smalkaram/Projects/travis20180117151534/TRIM/MIRDEEP_MAPPER/MIRDEEP2_run3/DESEQINT"
setwd(WD)

# Start Log
sink()
sink(file=paste(format(Sys.Date(), format="%m-%d-%Y"), format(Sys.time(), format="%H-%M-%S"), "log", sep="."), split=TRUE)

# Load Libraries & Functions
source("PACKAGES")
source("FUNCTIONS.R")

# Set the necessary Inputs
# MRNA_COUNTS_FILE = Counts file for all samples of mRNA generated using HISAT3 and STRINGTIE pipeline (See PRE_MRNA_ANALYSIS subdirectory)
# MiRNA_COUNTS_FILE = Counts file for all samples of miRNA generated using miRDeep2 pipeline (See PRE_MIRNA_ANALYSIS subdirectory)
# MiRNA_MAP_FILE = A mapping file that maps the sample names and the identifiers output from miRDeep2
# PVAL = Adjusted pvalue cutoff
# TYPE = can be either gene or transcript depending on your counts files
# GTF = file generated from the HISAT2 and STRINGTIE pipeline
# refGTF = reference GTF file for the genome assembly used for alignments
# MAX_TARGETS = maximum number of miRNA targets to output in DEresult tables
# CUTOFF_READ = cutoff # for minimum reads for any gene in any sample
# TYPE = "gene" # necessary for compatibility with pipeline scripts

TYPE <- "gene"
PVAL <- 0.05
refGTF <- "/data/db/GenCode/GRCm38p5/gencode.vM14.primary_assembly.annotation.gtf"
#GTF <- "Merged.gtf"
#MRNA_COUNTS_FILE <- "stringtie.gene.count.matrix"
MiRNA_COUNTS_FILE <- "/home/smalkaram/Projects/travis20180117151534/TRIM/MIRDEEP_MAPPER/MIRDEEP2_run3/quantify_all_mature_novel/MATURE/miRNAs_expressed_all_samples_20180122084334.csv"
MiRNA_MAP_FILE <- "/home/smalkaram/Projects/travis20180117151534/TRIM/MIRDEEP_MAPPER/mapping_file3.txt" 
Samples_FILE<- "/home/smalkaram/Projects/travis20180117151534/TRIM/MIRDEEP_MAPPER/MIRDEEP2_run3/DESEQINT/Samples.txt"
CUTOFF_READ <- 10
MAX_TARGETS<-20

PC(paste("Working directory", WD, sep=":"))

#----------------miRNA----------------------------------
PC("Analyzing miRNA")

# Reads MRNA Counts file
PC("Reading Counts file", 1)
CountsOrig <- read.delim(MiRNA_COUNTS_FILE, sep = "\t", header = TRUE)
CountsOrig <- CountsOrig[CountsOrig$total != 0, ]

# Clean up
PC("Cleanup gene names", 1)
CountsOrig$refname <- gsub("chr.*-", "", CountsOrig[, 1])
CountsOrig$reftotal <- paste(CountsOrig$refname, CountsOrig$total, 
		sep = "_")
CountsOrig <- takeFirstCommonRow(CountsOrig, "reftotal")
CountsOrig <- CountsOrig[, -grep("reftotal", colnames(CountsOrig))]
CountsOrig <- CountsOrig[, -grep("read_count|total|norm", colnames(CountsOrig), 
				perl = TRUE, ignore.case = TRUE)]
CountsOrig <- addByColFactor(CountsOrig, "refname")

# Read MiRNA map file
PC("Read miRDeep2 id to samples map", 1)
MAP <- read.table(MiRNA_MAP_FILE, header = FALSE)
NAMES <- MAP[, 2]
#-------This line needs to be checked as it is specific for each analysis
MAP <- gsub("_trimmed.fq.gz", "", MAP[, 1])
names(MAP) <- NAMES

# Change column names of original matrix to that in map, and 
# get new matrix with the total corresponding to samples represented in Samples.txt
PC("Change miRDeep2 id to sample names", 1)
colnames(CountsOrig) <- as.vector(MAP[colnames(CountsOrig)])
#-------This line needs to be checked as it is specific for each analysis
SAMPLES <- unique(gsub("_S[0-9]+_L[0-9]_[0-9]$", "", colnames(CountsOrig)))
Counts <- sapply(SAMPLES, function(x) GETSUM(x, CountsOrig))

Counts <- Counts[rowSums(Counts) > CUTOFF_READ, ]

# Read sample information file
PC("Reading sample information file", 1)
ColData<-read.table(Samples_FILE, header=TRUE, sep="\t", row.names=1)
fat<-as.factor(ColData$fat)
rep<-as.factor(ColData$rep)
lib<-as.factor(ColData$lib)
fatrep<-as.factor(ColData$fatrep)
fatlib<-as.factor(ColData$fatlib)
fatreplib<-as.factor(ColData$fatreplib)

PC("Compute DE", 1)
# Compute DE for all samples (for comparisons between fat DM, FF, LF)
PC("Computing DE between fat DM, FF, LF", 2)
Counts <- Counts[, rownames(ColData)]
PC("Samples", 2)
print(colnames(Counts))
dds <- DESeqDataSetFromMatrix(countData = Counts, colData = ColData, 
		design = formula(~ lib + fat ))
dds$fat<-relevel(dds$fat, ref="DM")
dds <- DESeq(dds)

# Comparison between Fat levels
PC("fat", 2)
FF_DM<-results(dds, contrast=c("fat", "FF", "DM"))
LF_DM<-results(dds, contrast=c("fat", "LF", "DM"))
FF_LF<-results(dds, contrast=c("fat", "FF", "LF"))


#Collect all above results into this list
RESULTS <-c('FF_DM', 'LF_DM', 'FF_LF')
PC("MiRNA result variables", 1)
PC(paste(RESULTS, collapse= " "), 2)

# Clean up results tables
# Replace NAs in pvalue and padj column with 1. Add foldchange and absolute foldchange column
# prefix the results with "MiRNA"
PC("Cleaning up results", 1)
for (VAR in RESULTS) {
	PC(VAR, 2)
	assign(VAR, replaceDFColNAs(get(VAR), "pvalue", 1))
	assign(VAR, replaceDFColNAs(get(VAR), "padj", 1))
	assign(VAR, addFC(get(VAR), logFCCOL = "log2FoldChange"))
	assign(VAR, addAbsFC(get(VAR), FCCOL = "foldChange"))
	assign(paste("MiRNA", VAR, sep="_"), addUpDown(get(VAR), log2FCCOL = "log2FoldChange"))
}

#-------------------------------
# Combine the results and find miRNA targets
PC("Combine DE results")

PC("Get a common list of MiRNA present in all results", 1)
# Get the common list of all MiRNA names from all DE results. 
MiRNA_QUERY <- unique(unlist(lapply(RESULTS, function(X) { rownames(get(paste("MiRNA", X, sep = "_"))) })))
#Get targets for this list
PC("Get targets", 2)
MiRNA_TARGET_LIST <- get.multimir(org = "mmu", mirna = MiRNA_QUERY, summary = FALSE, table = "all", add.link = FALSE, 
    use.tibble = FALSE, predicted.site = "conserved", predicted.cutoff = 10000, predicted.cutoff.type = "n")

#------------------------------
PC("Writing final result tables")
# Analyze the results of DESeq2 results for both MRNA and MiRNA
# And the targets 
# Provide a variable 'VAR', which has DESeq2 results
# Writes DE tables with extra columns as for both MRNA and MiRNA 
# at 0.05, 0.1 and 1(all) adjusted pvalue cutoffs sorting on absolute FoldChange 

# If this is on compute cluster with torque resource manager, 
# enable PARALLEL=TRUE only if its requirements, as below are met, otherwise use PARALLEL=FALSE



#----------------MRNA----------------------------------
# Read DESeq2 MRNA results from csv file
# The function ReadJimsDeseqMRNAresults is used to read the csv files  provided for this project
# The function reads and processes the specific coulumns and information contained therein, required for this analysis

# Read the three csv files required for three comparisons. 
# To be compatible with the pipeline, the result objects are named similar to the result objects of miRNA results, and prefixed with "MRNA"

MRNA_FF_DM<-ReadJimsDeseqMRNAresults("/home/smalkaram/Projects/travis20180117151534/mRNA_Results_from_Jim/significant_Fat_DMEM.csv")
MRNA_FF_LF<-ReadJimsDeseqMRNAresults("/home/smalkaram/Projects/travis20180117151534/mRNA_Results_from_Jim/significant_Fat_Lean.csv")
MRNA_LF_DM<-ReadJimsDeseqMRNAresults("/home/smalkaram/Projects/travis20180117151534/mRNA_Results_from_Jim/significant_Lean_DMEM.csv")
head(MRNA_LF_DM)



# Get Ensembl-MGI map
PC("Creating Ensembl-MGI map", 1)
library(biomaRt)
MART <- ""
while (class(MART) != "Mart") {
	try(MART <- useDataset(as.character("mmusculus_gene_ensembl"), mart = useMart("ENSEMBL_MART_ENSEMBL")))
}
GENE_MGI <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "mgi_symbol"),
		mart = MART)
GENE_MGI <- ReplaceEmptyMappings(GENE_MGI, "hgnc_symbol", "mgi_symbol")
GENE_MGI <- ReplaceEmptyMappings(GENE_MGI, "ensembl_gene_id", "mgi_symbol")
# Changed separation from , to ;
GENE_MGI <- unlist(lapply(Map2List(GENE_MGI, "ensembl_gene_id", "mgi_symbol"),
				function(x) paste(x, collapse = ";")))



#-----------------FINAL DE TABLES---------------------------------

# SERIAL PROCESSING EACH RESULT-------------------
# All results are stored in the current directory
PC("Writing results in serial mode", 1)

for (VAR in RESULTS) {
	PC(VAR, 1)
	MiRNA_VAR <- paste("MiRNA", VAR, sep = "_")
	MiRNA_DE <- get(MiRNA_VAR)

	# If you have MRNA----------------------------
	MRNA_VAR <- paste("MRNA", VAR, sep = "_")
	MRNA_DE <- get(MRNA_VAR)

	# IF there are novel genes in MRNA like the stringtie novel ones, remove them
	if(length(grep("^MSTRG", rownames(MRNA_DE))) > 0) {
		MRNA_DE <- MRNA_DE[-grep("^MSTRG", rownames(MRNA_DE)), ]
	}

	rownames(MRNA_DE) <- gsub("[.][0-9]*$", "", rownames(MRNA_DE))

	 
	# Don't need this as we are not writing mRNA tables
	#MRNA_QUERY <- rownames(MRNA_DE)

	MiRNA_QUERY <- rownames(MiRNA_DE)

	#MiRNA------------------------------------------------------
	PC('MiRNA', 2)
	PREV_COLUMNS <- colnames(MiRNA_DE)


	# If you DONT have MRNA DE data ------------------------------
	#ADD_COLUMNS <- c("validated", "predicted", "validated_freq", "predicted.freq")
	#TEMP <- cbind(as.data.frame(MiRNA_DE), t(sapply(MiRNA_QUERY, function(QUERY) {
	#	EXTRACT_TARGETS_NO_TARGET_DE(QUERY, MiRNA_TARGET_LIST, "mature_mirna_id", "target_ensembl",
	#		MAX_TARGETS) })))

	# If you have MRNA DE data ------------------------------
	ADD_COLUMNS <- c("val.sig.up", "val.sig.down", "val.nsig.up", "val.nsig.down", "all.other.val",
        "val.sig.up.freq", "val.sig.down.freq", "val.nsig.up.freq", "val.nsig.down.freq", "all.other.val.freq",
        "pred.sig.up", "pred.sig.down", "pred.nsig.up", "pred.nsig.down", "all.other.pred",
        "pred.sig.up.freq", "pred.sig.down.freq", "pred.nsig.up.freq", "pred.nsig.down.freq", "all.other.pred.freq")

	TEMP <- cbind(as.data.frame(MiRNA_DE), t(sapply(MiRNA_QUERY, function(QUERY) {
		EXTRACT_TARGETS(QUERY, MiRNA_TARGET_LIST, "mature_mirna_id", "target_ensembl",
			MAX_TARGETS, MRNA_DE) })))

	colnames(TEMP) <- c(PREV_COLUMNS, ADD_COLUMNS)
	WRITE_SIG_TABLES(DF = TEMP, NAME = MiRNA_VAR, PVAL = "padj", ALPHA = 0.05,
			FC = "absFoldChange")
	WRITE_SIG_TABLES(DF = TEMP, NAME = MiRNA_VAR, PVAL = "padj", ALPHA = 0.1,
			FC = "absFoldChange")
	WRITE_SIG_TABLES(DF = TEMP, NAME = MiRNA_VAR, PVAL = "padj", ALPHA = 1,
			FC = "absFoldChange")

}
# ------------------------------------




sink()
# Sridhar A Malkaram (smalkaram@wvstateu.edu)
# Last modified on: 12/30/2017
# Sridhar A Malkaram (smalkaram@wvstateu.edu)
# Last modified on: 12/30/2017
