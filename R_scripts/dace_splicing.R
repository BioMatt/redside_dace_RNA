# A script to run differential exon usage (DEU) analysis on redside dace RNAseq data, following: https://github.com/Oshlack/Lace/wiki/Example:-Differential-Transcript-Usage-on-a-non-model-organism

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("statmod")
library(edgeR)
library(tidyverse)
library(enrichR)
library(extrafont)
#font_import()
loadfonts()

#Read in data
counts <- read.table("dace_counts.txt",header=TRUE,sep="\t")

# Define groups
treatment <- c(rep("1", 10), rep("2", 10), rep("3", 10))
sample <- c('CTmax1', 'CTmax10', 'CTmax2', 'CTmax3', 'CTmax4', 'CTmax5', 'CTmax6', 'CTmax7', 'CTmax8', 'CTmax9', 'Handle1', 'Handle10', 'Handle2', 'Handle3', 'Handle4', 'Handle5', 'Handle6', 'Handle7', 'Handle8', 'Handle9', 'Wild1', 'Wild10', 'Wild2', 'Wild3', 'Wild4', 'Wild5', 'Wild6', 'Wild7', 'Wild8', 'Wild9')

# Read in a metadata file that includes RIN scores
metadata <- read_delim("dace_metadata.txt", delim = "\t")
RIN <- metadata$RIN

#Make DGElist and normalise
dx <- DGEList(counts[,c(7:36)])
dx<-calcNormFactors(dx,group=metadata$treatment)

#Make exon id
eid = paste0(counts$Chr,":",counts$Start)

#Define design matrix
design <- model.matrix(~treatment+RIN)

#Voom transform
vx <- voom(dx,design)
#Fit with limma
fx <- lmFit(vx,design)

ex <- diffSplice(fx,geneid=counts$Chr,exonid=eid)
res<-topSplice(ex,number="20")
plotSplice(ex) #Plots the top gene

############################################################################
# Many superTranscripts have way too many transcripts within them, so we try filtering out the superTranscripts with that greater number and see 
dup_counts <- counts
dup_counts <- dup_counts %>% 
  separate(Chr, c("SuperTranscript", "transcript"), sep = "\\.")

# Based on this tally, a single SuperTranscript, Cluster-41845 has 200,212 transcripts. The next highest one has only 744 transcripts. 
supertranscript_tally <- dup_counts %>% 
  group_by(SuperTranscript) %>% 
  tally()

# Filtering out Cluster-41845. Removing this cluster cut down observations from 401541 to 201329
# Also removing all supertranscripts with >=500 transcripts within them
supertranscript_tally <- dplyr::filter(supertranscript_tally, n < 500)
dup_counts <- dplyr::filter(dup_counts, SuperTranscript %in% c(supertranscript_tally$SuperTranscript))
# Reverting the count table to the same number of columns as Counts for subsetting down the road
dup_counts <- dplyr::select(dup_counts, -transcript)

#Make exon id
eid = paste0(dup_counts$Geneid,":",dup_counts$Start)

########################################################
# Breaking apart the overall model into 3 pairwise comparisons
# CTmax versus Handling, first
CvH_counts <- DGEList(dup_counts[,c(7:26)])
CvH_treatment <- c(rep("1", 10), rep("2", 10))
CvH_RIN <- RIN[1:20]
CvH_dx<-calcNormFactors(CvH_counts,group=CvH_treatment)

CvH_design <- model.matrix(~CvH_treatment + CvH_RIN)

#Voom transform
CvH_vx <- voom(CvH_dx,CvH_design)
#Fit with limma
CvH_fx <- lmFit(CvH_vx,CvH_design)

CvH_ex <- diffSplice(CvH_fx,geneid=dup_counts$Geneid,exonid=eid, robust = TRUE)
plotSplice(CvH_ex)
CvH_res <- topSplice(CvH_ex, FDR = 0.05)

# Now, Wild vs CTmax
WvC_counts <- DGEList(dup_counts[,c(27:36,7:16)])
WvC_treatment <- c(rep("1", 10), rep("2", 10))
WvC_RIN <- RIN[c(21:30, 1:10)]
WvC_dx<-calcNormFactors(WvC_counts,group=WvC_treatment)

WvC_design <- model.matrix(~WvC_treatment+WvC_RIN)

#Voom transform
WvC_vx <- voom(WvC_dx,WvC_design)
#Fit with limma
WvC_fx <- lmFit(WvC_vx,WvC_design)

WvC_ex <- diffSplice(WvC_fx,geneid=dup_counts$Geneid,exonid=eid, robust = TRUE)
plotSplice(WvC_ex)
WvC_res <- topSplice(WvC_ex, FDR = 0.05)

# Finally, Wild v Handling
WvH_counts <- DGEList(dup_counts[,c(27:36,17:26)])
WvH_treatment <- c(rep("1", 10), rep("2", 10))
WvH_RIN <- RIN[c(21:30, 11:20)]
WvH_dx<-calcNormFactors(WvH_counts,group=WvH_treatment)

WvH_design <- model.matrix(~WvH_treatment+WvH_RIN)

#Voom transform
WvH_vx <- voom(WvH_dx,WvH_design)
#Fit with limma
WvH_fx <- lmFit(WvH_vx,WvH_design)

WvH_ex <- diffSplice(WvH_fx,geneid=dup_counts$Geneid,exonid=eid, robust = TRUE)
plotSplice(WvH_ex)
WvH_res <- topSplice(WvH_ex, FDR = 0.05)

# Write out the 3 results files
#write_delim(CvH_res, "CvH_DEU_res.txt", delim = "\t")
#write_delim(WvC_res, "WvC_DEU_res.txt", delim = "\t")
#write_delim(WvH_res, "WvH_DEU_res.txt", delim = "\t")

#############################################################################################################
# Trying DEXseq because it gives significance on an exon-by-exon basis
library(DEXSeq)
library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
library(doParallel)
cores<-detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
library(GenomicFeatures)

BPPARAM = SnowParam(13)

metadata <- read.delim("dace_metadata.txt", sep = "\t")
metadata$condition <- as.factor(metadata$treatment)
# Center and scale RIN scores to lower standard deviation and help with GLM convergence
metadata$scaledRIN <- scale(metadata$RIN)

counts_dex <- read.table("dace_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = F)

clusters <- counts_dex[,2]
counts_only <- counts_dex[,7:36]

# Convert data to integer
cc <- data.matrix(counts_only)
cc = round(cc)
colnames(cc) <- metadata$sample

# Unique exon names, using gene names and the starting base pair for the exon
exon_ids <- paste0(counts_dex[,1], "_", counts_dex[,3])

# Create the dataset and design formula. Scaled RIN is included to take into account RNA integrity
dxd <- DEXSeqDataSet(cc, design =~sample + scaledRIN + exon + condition:exon, featureID = as.factor(exon_ids), groupID = as.factor(clusters), sampleData = metadata, featureRanges = myGranges)

##Estimate size factors and dispersions (DEXseq does this based on a negative bionmial distribution
dxd <- DEXSeq::estimateSizeFactors(dxd)
dxd <- DEXSeq::estimateDispersions(dxd)

#Test for DEU
dxd <- testForDEU(dxd, BPPARAM=SerialParam())


# Estimating exon fold changes
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")

# Importing my GFF for sending plot information to DEXSeq

#gffRangedData<-import.gff("dexseq.gff")
#myGranges<-as(gffRangedData, "GRanges")
#test <- dxd
#test@rowRanges <- myGranges
#dex_res$genomicData <- myGranges

txdb <- makeTxDbFromGFF("dexseq.gff", format = "gff3")

#Extract the results
dex_res <- DEXSeqResults(dxd)
table ( dex_res$padj < 0.05 )

results <- as_tibble(dex_res)
# Check how many clusters show any DEU
results_filtered <- filter(results, padj < 0.05)
length(unique(results_filtered$groupID))
# An exploratory MA plot
#DEXSeq::plotMA(dex_res)

DEXSeq::plotDEXSeq(object = dex_res, geneID = "Cluster-43916.5", fitExpToVar = "condition", norCounts = T, legend = T)
mtext("HSP 70", line = 55, cex = 2.0)

DEXSeq::plotDEXSeq(object = dex_res, geneID = "Cluster-80266.0", fitExpToVar = "condition", norCounts = T, legend = T)
mtext("Apolipoprotein L6", line = 55, cex = 2.0)

DEXSeq::plotDEXSeq(object = dex_res, geneID = "Cluster-65835.1", fitExpToVar = "condition", norCounts = T, legend = T)
mtext("Palladin", line = 55, cex = 2.0)


# Srek1
DEXSeq::plotDEXSeq(object = dex_res, geneID = "Cluster-50970.0",  fitExpToVar = "condition", norCounts = FALSE, legend = F, splicing = T, expression = F, displayTranscripts=TRUE, color = c("#4575b4", "#fdae61", "#d73027"), family = "Arial")

DEXSeq::plotDEXSeq(object = dex_res, geneID = "Cluster-36495.0",  fitExpToVar = "condition", norCounts = FALSE, legend = F, splicing = T, expression = F, displayTranscripts=TRUE, color = c("#4575b4", "#fdae61", "#d73027"), family = "Arial")
mtext("RCC1-like G exchanging factor-like protein", line = 45, cex = 2.0)

DEXSeq::plotDEXSeq(object = dex_res, geneID = "Cluster-78257.2", fitExpToVar = "condition", norCounts = T, legend = T, fitExpToVar = "condition", norCounts = FALSE, legend = F, splicing = T, expression = F, displayTranscripts=TRUE, color = c("#4575b4", "#fdae61", "#d73027"), family = "Arial")
mtext("RNA-binding protein 39", line = 47, cex = 2.0)

DEXSeq::plotDEXSeq(object = dex_res, geneID = "Cluster-68486.1", fitExpToVar = "condition", norCounts = FALSE, legend = F, splicing = T, expression = F, displayTranscripts=TRUE, color = c("#4575b4", "#fdae61", "#d73027"), family = "Arial")
mtext("Pinin", line = 50, cex = 2.0)

DEXSeq::plotDEXSeq(object = dex_res, geneID = "Cluster-5910.1", fitExpToVar = "condition", norCounts = FALSE, legend = F, splicing = T, expression = F, displayTranscripts=TRUE, color = c("#4575b4", "#fdae61", "#d73027"), family = "Arial")
mtext("RNA-binding protein 25", line = 50, cex = 2.0)




dev.off()
plotDispEsts(dxd)
plotMA(dxd)

###############################################################################################################################
# Doing a functional analysis of the DEU results
# Read in the corset clusters, relating supertranscript to transcript from Trinity, and add informative column names.
corset_clusters <- read_delim("dace-clusters.txt", delim = "\t", col_names = FALSE)
corset_clusters <- rename(corset_clusters, transcript_id = X1, corset_cluster = X2)

# Read in the transcript annotations
annotations <- read_delim("dace_transcriptome_annotations.txt", delim = "\t")

# First, combine the corset clusters and annotation sheets by Trinity transcript ID
combined_data <- full_join(corset_clusters, annotations)
# Add the data to the combined spreadsheet.
combined_data <- combined_data %>% 
  right_join(results_filtered, by = c("corset_cluster" = "groupID"))


#################################################################################################################################
# Using EnrichR to look at GO terms
# Taking a look at what databases are available in EnrichR
listEnrichrDbs()
# A function to take an EnrichrR table and split the GO term and definition for input into Revigo
split_go <- function(x) {
  # Taking apart the GO descriptions and GO ID terms
  x <- x %>%
    select(-starts_with("Old")) %>%
    separate(Term, c("GO_term", "GO_ID"), sep = "GO")
  
  # Changing all the floating :#### GO terms to be GO:####
  x <- mutate_if(x, is.character, str_replace_all, pattern = ":", replacement = "GO:")
  
  # Removing the last parentheses in the GO ID 
  x$GO_ID <- x$GO_ID %>% 
    str_replace("\\)", "")
  
  # Removing the last parentheses in the GO term to have clean looking cells 
  x$GO_term <- x$GO_term %>% 
    str_replace("\\($", "")
  return(x)
}
# Run enrichR, searching the biological process and molecular function databases
pos_ctmax_deu <- dplyr::filter(combined_data, log2fold_Wild_CT_max < 0 & log2fold_Handle_CT_max < 0)
length(unique(pos_ctmax_deu$corset_cluster))

length(unique(pos_ctmax_deu$gene_name))
pos_CTmax_enriched <- enrichr(pos_ctmax_deu$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

pos_CTmax_enriched_bioGO <- filter(split_go(pos_CTmax_enriched$GO_Biological_Process_2018), Adjusted.P.value < 0.05)
pos_CTmax_enriched_molGO <- filter(split_go(pos_CTmax_enriched$GO_Molecular_Function_2018), Adjusted.P.value < 0.05)
pos_CTmax_enriched_cellGO <- filter(split_go(pos_CTmax_enriched$GO_Cellular_Component_2018), Adjusted.P.value < 0.05)


#################################################################################################################################
neg_ctmax_deu <- dplyr::filter(combined_data, log2fold_Wild_CT_max > 0 & log2fold_Handle_CT_max > 0)
length(unique(neg_ctmax_deu$gene_name))
length(unique(neg_ctmax_deu$corset_cluster))

# Run enrichR, searching the biological process and molecular function databases
neg_CTmax_enriched <- enrichr(neg_ctmax_deu$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

neg_CTmax_enriched_bioGO <- filter(split_go(neg_CTmax_enriched$GO_Biological_Process_2018), Adjusted.P.value < 0.05)
neg_CTmax_enriched_molGO <- filter(split_go(neg_CTmax_enriched$GO_Molecular_Function_2018), Adjusted.P.value < 0.05)
neg_CTmax_enriched_cellGO <- filter(split_go(neg_CTmax_enriched$GO_Cellular_Component_2018), Adjusted.P.value < 0.05)

#################################################################################################################################
overall_enriched <- enrichr(combined_data$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))
overall_enriched_bioGO <- filter(split_go(overall_enriched$GO_Biological_Process_2018), Adjusted.P.value < 0.05)
overall_enriched_molGO <- filter(split_go(overall_enriched$GO_Molecular_Function_2018), Adjusted.P.value < 0.05)
overall_enriched_cellGO <- filter(split_go(overall_enriched$GO_Cellular_Component_2018), Adjusted.P.value < 0.05)

###############################################################################################################################
write_delim(distinct(pos_ctmax_deu, gene_name, .keep_all = TRUE), "pos_CTmax_unique_DEU.txt", delim = "\t")
write_delim(distinct(neg_ctmax_deu, gene_name, .keep_all = TRUE), "neg_CTmax_unique_DEU.txt", delim = "\t")
write_delim(distinct(combined_data, gene_name, .keep_all = TRUE), "overall_unique_DEU.txt", delim = "\t")

write_delim(pos_CTmax_enriched_bioGO, "pos_CTmax_enriched_bioGO.txt", delim = "\t")
write_delim(pos_CTmax_enriched_molGO, "pos_CTmax_enriched_molGO.txt", delim = "\t")
write_delim(pos_CTmax_enriched_cellGO, "pos_CTmax_enriched_cellGO.txt", delim = "\t")

write_delim(neg_CTmax_enriched_bioGO, "neg_CTmax_enriched_bioGO.txt", delim = "\t")
write_delim(neg_CTmax_enriched_molGO, "neg_CTmax_enriched_molGO.txt", delim = "\t")
write_delim(neg_CTmax_enriched_cellGO, "neg_CTmax_enriched_cellGO.txt", delim = "\t")

write_delim(overall_enriched_bioGO, "overall_enriched_bioGO.txt", delim = "\t")
write_delim(overall_enriched_molGO, "overall_enriched_molGO.txt", delim = "\t")
write_delim(overall_enriched_cellGO, "overall_enriched_cellGO.txt", delim = "\t")
