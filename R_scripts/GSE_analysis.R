# This script functionally analyzes differential gene expression found using EdgeR, by combining significant clusters (transcripts) with annotations found using trinotate
library(tidyverse)
library(enrichR)

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

# Read in the corset clusters, relating supertranscript to transcript from Trinity, and add informative column names.
corset_clusters <- read_delim("dace-clusters.txt", delim = "\t", col_names = FALSE)
corset_clusters <- rename(corset_clusters, transcript_id = X1, corset_cluster = X2)

# Read in the transcript annotations
annotations <- read_delim("dace_transcriptome_annotations.txt", delim = "\t")

# Read in the Wild versus CTmax differential gene expression (DGE) results
WvC <- read_delim("dace_WvC_results.txt", delim = "\t")
# Read the Wild versus Handling DGE results
WvH <- read_delim("dace_WvH_results.txt", delim = "\t")
# Last, the CTmax versus Handling results
CvH <- read_delim("dace_CvH_results.txt", delim = "\t")

#Pulling the WvC columns we're interested in- corset cluster, log fold change, and false discovery rate. Rename to include the Wild versus CTmax line (WvC) so the column headers make sense when combined with other data
# Also, filtering all these data for FDR < 0.05 so only significant results get included in the combined data
WvC <- dplyr::select(WvC, corset_cluster, WvC.logFC = logFC, WvC.FDR = FDR)
WvC <- dplyr::filter(WvC, WvC.FDR < 0.05)
# Follow the same process for Wild versus Handling data
WvH <- dplyr::select(WvH, corset_cluster, WvH.logFC = logFC, WvH.FDR = FDR)
WvH <- dplyr::filter(WvH, WvH.FDR < 0.05)
# Doing this a final time with the CTmax versus Handling data
CvH <- dplyr::select(CvH, corset_cluster, CvH.logFC = logFC, CvH.FDR = FDR)
CvH <- dplyr::filter(CvH, CvH.FDR < 0.05)


# Now, read in the clusters (transcripts) that only show DGE under CTmax conditions. First, the ones showing positive DGE
# Adding a column of 1's next to each cluster. When these tables are joined to the combined data, the clusters showing positive or negative DGE under CTmax will have a 1 in the new column, allowing for filtering for gene names downstream.
pos_CTmax <- read_delim("dace_pos_CTmax_clusters.txt", delim = "\t")
pos_CTmax <- mutate(pos_CTmax, pos_CTmax = 1)
# Now the clusters that showed negative DGE only during CTmax, not during handling or in the wild
neg_CTmax <- read_delim("dace_neg_CTmax_clusters.txt", delim = "\t")
neg_CTmax <- mutate(neg_CTmax, neg_CTmax = 1)

########################################################################
# With all data read in, we start combining data sets
# First, combine the corset clusters and annotation sheets by Trinity transcript ID
combined_data <- full_join(corset_clusters, annotations)

# Add the data to the combined spreadsheet.
combined_data <- combined_data %>% 
  left_join(WvC) %>% 
  left_join(WvH) %>% 
  left_join(CvH)

# Add the positive or negative DGE under CTMax only conditions
combined_data <- combined_data %>% 
  left_join(pos_CTmax) %>% 
  left_join(neg_CTmax)

#################################################################################
# Filtering out clusters for which no information exists. That is, filtering out those without annotations and no DGE in any comparison
combined_data <- dplyr::filter(combined_data, !is.na(gene_name) | !is.na(transcript_description) | !is.na(WvC.logFC) | !is.na(WvH.logFC) | !is.na(CvH.logFC))

# Tally can be used to count up observations within groups
combined_data %>%
  group_by(pos_CTmax) %>% 
  tally()


# Counting the number of rows under different conditions
dplyr::filter(combined_data, !is.na(WvC.logFC), !is.na(gene_name)) %>% 
  distinct(gene_name) %>% 
  nrow()
dplyr::filter(combined_data, !is.na(WvH.logFC), !is.na(gene_name)) %>% 
  distinct(gene_name) %>% 
  nrow()
dplyr::filter(combined_data, !is.na(CvH.logFC), !is.na(gene_name)) %>% 
  distinct(gene_name) %>% 
  nrow()
####################################################################################
# Writing out CSV files showing differential gene expression under the three different comparisons
write_csv(dplyr::filter(combined_data, !is.na(WvC.logFC)), "WildvCTmax_dge.csv", col_names = TRUE)
write_csv(dplyr::filter(combined_data, !is.na(WvH.logFC)), "WildvHandling_dge.csv", col_names = TRUE)
write_csv(dplyr::filter(combined_data, !is.na(CvH.logFC)), "CTmaxvHandling_dge.csv", col_names = TRUE)

# Also writing positive and negative CTmax genes
write_csv(dplyr::filter(combined_data, pos_CTmax == 1 & !is.na(gene_name)), "PosCTmax_dge.csv", col_names = TRUE)
write_csv(dplyr::filter(combined_data, neg_CTmax == 1 & !is.na(gene_name)), "NegCTmax_dge.csv", col_names = TRUE)

# Creating simplified results tables with just the unique super transcripts (corset_clusters), because statistical DGE results were found from those. But when multiple clusters have annotation info, keep just one annotation
simplified_WvH <- dplyr::filter(combined_data, !is.na(WvH.logFC)) %>% 
  group_by(corset_cluster) %>% 
  filter(rank(uniprot_id, ties.method="first")==1)
write_csv(simplified_WvH, "WildvHandling_simplified_dge.csv", col_names = TRUE)

simplified_WvC <- dplyr::filter(combined_data, !is.na(WvC.logFC)) %>% 
  group_by(corset_cluster) %>% 
  filter(rank(uniprot_id, ties.method="first")==1)
write_csv(simplified_WvC, "WildvCTmax_simplified_dge.csv", col_names = TRUE)

simplified_CvH <- dplyr::filter(combined_data, !is.na(CvH.logFC)) %>% 
  group_by(corset_cluster) %>% 
  filter(rank(uniprot_id, ties.method="first")==1)
write_csv(simplified_CvH, "CTmaxvHandling_simplified_dge.csv", col_names = TRUE)

#####################################################################################
# With data combined and filtered, we use enrichR to check for gene set enrichment
# First, taking the unique known genes that showed positive DGE in CTmax conditions
pos_genes <- dplyr::filter(combined_data, pos_CTmax == 1 & !is.na(gene_name))
pos_genes <- distinct(pos_genes, gene_name)
pos_genes <- as.vector(pos_genes)

# Taking a look at what databases are available in EnrichR
listEnrichrDbs()
# Run enrichR, searching the biological process and molecular function databases
pos_enriched <- enrichr(pos_genes$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
pos_enrich_bio <- as_tibble(pos_enriched[["GO_Biological_Process_2018"]])
pos_enrich_bio <- dplyr::filter(pos_enrich_bio, Adjusted.P.value < 0.05)
print(pos_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
pos_enrich_bio <- split_go(pos_enrich_bio)

pos_enrich_mol <- as_tibble(pos_enriched[["GO_Molecular_Function_2018"]])
pos_enrich_mol <- dplyr::filter(pos_enrich_mol, Adjusted.P.value < 0.05)
print(pos_enrich_mol$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
pos_enrich_mol <- split_go(pos_enrich_mol)

pos_enrich_comp <- as_tibble(pos_enriched[["GO_Cellular_Component_2018"]])
pos_enrich_comp <- dplyr::filter(pos_enrich_comp, Adjusted.P.value < 0.05)
print(pos_enrich_comp$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
pos_enrich_comp <- split_go(pos_enrich_comp)


# Writing just the GO terms to a table for input into Revigo
write.table(as.data.frame(rbind(cbind(pos_enrich_bio$GO_ID, pos_enrich_bio$Adjusted.P.value), cbind(pos_enrich_mol$GO_ID, pos_enrich_mol$Adjusted.P.value), cbind(pos_enrich_comp$GO_ID, pos_enrich_comp$Adjusted.P.value))), file = "pos_CTmax_GO_terms.txt", sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

# Writing out the whole GO term spreadsheets into their own files
write_delim(as.data.frame(rbind(pos_enrich_bio, pos_enrich_mol, pos_enrich_comp)), "pos_CTmax_GO_terms_full.txt", delim = "\t", col_names = TRUE)

View(combined_data %>% 
  group_by(corset_cluster, gene_id) %>% 
  tally()
)

View(combined_data %>% 
       group_by(gene_id, corset_cluster) %>% 
       tally()
)

# Following the same process for genes showing negative DGE in response to CTmax conditions
neg_genes <- dplyr::filter(combined_data, neg_CTmax == 1 & !is.na(gene_name))
neg_genes <- dplyr::distinct(neg_genes, gene_name)
neg_genes <- as.vector(neg_genes)

# Run enrichR, searching the biological process and molecular function databases
neg_enriched <- enrichr(neg_genes$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
neg_enrich_bio <- as_tibble(neg_enriched[["GO_Biological_Process_2018"]])
neg_enrich_bio <- dplyr::filter(neg_enrich_bio, Adjusted.P.value < 0.05)
print(neg_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
neg_enrich_bio <- split_go(neg_enrich_bio)

neg_enrich_mol <- as_tibble(neg_enriched[["GO_Molecular_Function_2018"]])
neg_enrich_mol <- dplyr::filter(neg_enrich_mol, Adjusted.P.value < 0.05)
print(neg_enrich_mol$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
neg_enrich_mol <- split_go(neg_enrich_mol)

neg_enrich_comp <- as_tibble(neg_enriched[["GO_Cellular_Component_2018"]])
neg_enrich_comp <- dplyr::filter(neg_enrich_comp, Adjusted.P.value < 0.05)
print(neg_enrich_comp$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
neg_enrich_comp <- split_go(neg_enrich_comp)

# Following the same process for genes showing DGE in response to CTmax conditions, positive or negative
CTmax_genes <- dplyr::filter(combined_data, (neg_CTmax == 1 | pos_CTmax == 1) & !is.na(gene_name))
CTmax_genes <- dplyr::distinct(CTmax_genes, gene_name)
CTmax_genes <- as.vector(CTmax_genes)

# Run enrichR, searching the biological process and molecular function databases
CTmax_enriched <- enrichr(CTmax_genes$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
CTmax_enrich_bio <- as_tibble(CTmax_enriched[["GO_Biological_Process_2018"]])
CTmax_enrich_bio <- dplyr::filter(CTmax_enrich_bio, Adjusted.P.value < 0.05)
print(CTmax_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
CTmax_enrich_bio <- split_go(CTmax_enrich_bio)

CTmax_enrich_mol <- as_tibble(CTmax_enriched[["GO_Molecular_Function_2018"]])
CTmax_enrich_mol <- dplyr::filter(CTmax_enrich_mol, Adjusted.P.value < 0.05)
print(CTmax_enrich_mol$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
CTmax_enrich_mol <- split_go(CTmax_enrich_mol)

CTmax_enrich_comp <- as_tibble(CTmax_enriched[["GO_Cellular_Component_2018"]])
CTmax_enrich_comp <- dplyr::filter(CTmax_enrich_comp, Adjusted.P.value < 0.05)
print(CTmax_enrich_comp$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
CTmax_enrich_comp <- split_go(CTmax_enrich_comp)


###################################################################################################################
# Pulling GO terms for the 3 treatments
# First, CTmax vs Handling genes
# Ones showing positive DGE here
CvH_pos_genes <- dplyr::filter(combined_data, CvH.logFC > 0 & !is.na(gene_name))
CvH_pos_genes <- dplyr::distinct(CvH_pos_genes, gene_name)
CvH_pos_genes <- as.vector(CvH_pos_genes)

# Run enrichR, searching the biological process and molecular function databases
CvH_pos_enriched <- enrichr(CvH_pos_genes$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))
nrow(filter(CvH_pos_enriched$GO_Biological_Process_2018, Adjusted.P.value<0.05))
nrow(filter(CvH_pos_enriched$GO_Molecular_Function_2018, Adjusted.P.value<0.05))
nrow(filter(CvH_pos_enriched$GO_Cellular_Component_2018, Adjusted.P.value<0.05))

# Pull the results from one of the databases, biological process first
CvH_pos_enrich_bio <- as_tibble(CvH_pos_enriched[["GO_Biological_Process_2018"]])
CvH_pos_enrich_bio <- dplyr::filter(CvH_pos_enrich_bio, Adjusted.P.value < 0.05)
print(CvH_pos_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
CvH_pos_enrich_bio <- split_go(CvH_pos_enrich_bio)

CvH_pos_enrich_mol <- as_tibble(CvH_pos_enriched[["GO_Molecular_Function_2018"]])
CvH_pos_enrich_mol <- dplyr::filter(CvH_pos_enrich_mol, Adjusted.P.value < 0.05)
CvH_pos_enrich_mol <- split_go(CvH_pos_enrich_mol)

CvH_pos_enrich_comp <- as_tibble(CvH_pos_enriched[["GO_Cellular_Component_2018"]])
CvH_pos_enrich_comp <- dplyr::filter(CvH_pos_enrich_comp, Adjusted.P.value < 0.05)
CvH_pos_enrich_comp <- split_go(CvH_pos_enrich_comp)

# Now genes showing negative DGE under CTmax vs Handling conditions
CvH_neg_genes <- dplyr::filter(combined_data, CvH.logFC < 0 & !is.na(gene_name))
CvH_neg_genes <- dplyr::distinct(CvH_neg_genes, gene_name)
CvH_neg_genes <- as.vector(CvH_neg_genes)

# Run enrichR, searching the biological process and molecular function databases
CvH_neg_enriched <- enrichr(CvH_neg_genes$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))
nrow(filter(CvH_neg_enriched$GO_Biological_Process_2018, Adjusted.P.value<0.05))
nrow(filter(CvH_neg_enriched$GO_Molecular_Function_2018, Adjusted.P.value<0.05))
nrow(filter(CvH_neg_enriched$GO_Cellular_Component_2018, Adjusted.P.value<0.05))
# Pull the results from one of the databases, biological process first
CvH_neg_enrich_bio <- as_tibble(CvH_neg_enriched[["GO_Biological_Process_2018"]])
CvH_neg_enrich_bio <- dplyr::filter(CvH_neg_enrich_bio, Adjusted.P.value < 0.05)
print(CvH_neg_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
CvH_neg_enrich_bio <- split_go(CvH_neg_enrich_bio)

CvH_neg_enrich_mol <- as_tibble(CvH_neg_enriched[["GO_Molecular_Function_2018"]])
CvH_neg_enrich_mol <- dplyr::filter(CvH_neg_enrich_mol, Adjusted.P.value < 0.05)
CvH_neg_enrich_mol <- split_go(CvH_neg_enrich_mol)

CvH_neg_enrich_comp <- as_tibble(CvH_neg_enriched[["GO_Cellular_Component_2018"]])
CvH_neg_enrich_comp <- dplyr::filter(CvH_neg_enrich_comp, Adjusted.P.value < 0.05)
CvH_neg_enrich_comp <- split_go(CvH_neg_enrich_comp)

# Following the same process for genes showing positive and negative DGE, but for the Wild vs Handling condition
WvH_pos_genes <- dplyr::filter(combined_data, WvH.logFC > 0 & !is.na(gene_name))
WvH_pos_genes <- dplyr::distinct(WvH_pos_genes, gene_name)
WvH_pos_genes <- as.vector(WvH_pos_genes)

# Run enrichR, searching the biological process and molecular function databases
WvH_pos_enriched <- enrichr(WvH_pos_genes$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
WvH_pos_enrich_bio <- as_tibble(WvH_pos_enriched[["GO_Biological_Process_2018"]])
WvH_pos_enrich_bio <- dplyr::filter(WvH_pos_enrich_bio, Adjusted.P.value < 0.05)
print(WvH_pos_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvH_pos_enrich_bio <- split_go(WvH_pos_enrich_bio)

WvH_pos_enrich_mol <- as_tibble(WvH_pos_enriched[["GO_Molecular_Function_2018"]])
WvH_pos_enrich_mol <- dplyr::filter(WvH_pos_enrich_mol, Adjusted.P.value < 0.05)
print(WvH_pos_enrich_mol$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvH_pos_enrich_mol <- split_go(WvH_pos_enrich_mol)

WvH_pos_enrich_comp <- as_tibble(WvH_pos_enriched[["GO_Cellular_Component_2018"]])
WvH_pos_enrich_comp <- dplyr::filter(WvH_pos_enrich_comp, Adjusted.P.value < 0.05)
print(WvH_pos_enrich_comp$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvH_pos_enrich_comp <- split_go(WvH_pos_enrich_comp)

# Now genes showing negative DGE under Wild vs Handling conditions
WvH_neg_genes <- dplyr::filter(combined_data, WvH.logFC < 0 & !is.na(gene_name))
WvH_neg_genes <- dplyr::distinct(WvH_neg_genes, gene_name)
WvH_neg_genes <- as.vector(WvH_neg_genes)

# Run enrichR, searching the biological process and molecular function databases
WvH_neg_enriched <- enrichr(WvH_neg_genes$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
WvH_neg_enrich_bio <- as_tibble(WvH_neg_enriched[["GO_Biological_Process_2018"]])
WvH_neg_enrich_bio <- dplyr::filter(WvH_neg_enrich_bio, Adjusted.P.value < 0.05)
print(WvH_neg_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvH_neg_enrich_bio <- split_go(WvH_neg_enrich_bio)

WvH_neg_enrich_mol <- as_tibble(WvH_neg_enriched[["GO_Molecular_Function_2018"]])
WvH_neg_enrich_mol <- dplyr::filter(WvH_neg_enrich_mol, Adjusted.P.value < 0.05)
print(WvH_neg_enrich_mol$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvH_neg_enrich_mol <- split_go(WvH_neg_enrich_mol)

WvH_neg_enrich_comp <- as_tibble(WvH_neg_enriched[["GO_Cellular_Component_2018"]])
WvH_neg_enrich_comp <- dplyr::filter(WvH_neg_enrich_comp, Adjusted.P.value < 0.05)
print(WvH_neg_enrich_comp$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvH_neg_enrich_comp <- split_go(WvH_neg_enrich_comp)

# Following the same process for genes showing positive and negative DGE, but for the Wild vs CTmax condition
WvC_pos_genes <- dplyr::filter(combined_data, WvC.logFC > 0 & !is.na(gene_name))
WvC_pos_genes <- dplyr::distinct(WvC_pos_genes, gene_name)
WvC_pos_genes <- as.vector(WvC_pos_genes)
length(WvC_pos_genes$gene_name)
# Run enrichR, searching the biological process and molecular function databases
WvC_pos_enriched <- enrichr(WvC_pos_genes$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
WvC_pos_enrich_bio <- as_tibble(WvC_pos_enriched[["GO_Biological_Process_2018"]])
WvC_pos_enrich_bio <- dplyr::filter(WvC_pos_enrich_bio, Adjusted.P.value < 0.05)
print(WvC_pos_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvC_pos_enrich_bio <- split_go(WvC_pos_enrich_bio)

WvC_pos_enrich_mol <- as_tibble(WvC_pos_enriched[["GO_Molecular_Function_2018"]])
WvC_pos_enrich_mol <- dplyr::filter(WvC_pos_enrich_mol, Adjusted.P.value < 0.05)
print(WvC_pos_enrich_mol$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvC_pos_enrich_mol <- split_go(WvC_pos_enrich_mol)

WvC_pos_enrich_comp <- as_tibble(WvC_pos_enriched[["GO_Cellular_Component_2018"]])
WvC_pos_enrich_comp <- dplyr::filter(WvC_pos_enrich_comp, Adjusted.P.value < 0.05)
print(WvC_pos_enrich_comp$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvC_pos_enrich_comp <- split_go(WvC_pos_enrich_comp)

# Now genes showing negative DGE under Wild vs CTmax conditions
WvC_neg_genes <- dplyr::filter(combined_data, WvC.logFC < 0 & !is.na(gene_name))
WvC_neg_genes <- dplyr::distinct(WvC_neg_genes, gene_name)
WvC_neg_genes <- as.vector(WvC_neg_genes)
length(WvC_neg_genes$gene_name)

# Run enrichR, searching the biological process and molecular function databases
WvC_neg_enriched <- enrichr(WvC_neg_genes$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
WvC_neg_enrich_bio <- as_tibble(WvC_neg_enriched[["GO_Biological_Process_2018"]])
WvC_neg_enrich_bio <- dplyr::filter(WvC_neg_enrich_bio, Adjusted.P.value < 0.05)
print(WvC_neg_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvC_neg_enrich_bio <- split_go(WvC_neg_enrich_bio)

WvC_neg_enrich_mol <- as_tibble(WvC_neg_enriched[["GO_Molecular_Function_2018"]])
WvC_neg_enrich_mol <- dplyr::filter(WvC_neg_enrich_mol, Adjusted.P.value < 0.05)
print(WvC_neg_enrich_mol$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvC_neg_enrich_mol <- split_go(WvC_neg_enrich_mol)

WvC_neg_enrich_comp <- as_tibble(WvC_neg_enriched[["GO_Cellular_Component_2018"]])
WvC_neg_enrich_comp <- dplyr::filter(WvC_neg_enrich_comp, Adjusted.P.value < 0.05)
print(WvC_neg_enrich_comp$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvC_neg_enrich_comp <- split_go(WvC_neg_enrich_comp)

# Turn the scientific notation to decimal for input into Revigo
options(scipen=999)
# Use options(scipen=0) to turn scientific notation back on
# Writing out the 6 GO term spreadsheets for input into Revigo. Only including the GO ID and adjusted p value, with no column names
write_delim(as.data.frame(rbind(cbind(WvH_neg_enrich_bio$GO_ID, WvH_neg_enrich_bio$Adjusted.P.value), cbind(WvH_neg_enrich_mol$GO_ID, WvH_neg_enrich_mol$Adjusted.P.value), cbind(WvH_neg_enrich_comp$GO_ID, WvH_neg_enrich_comp$Adjusted.P.value))), "WvH_neg_GO.txt", delim = "\t", col_names = FALSE)
# WvH pos is empty
write_delim(as.data.frame(rbind(cbind(WvC_neg_enrich_bio$GO_ID, WvC_neg_enrich_bio$Adjusted.P.value), cbind(WvC_neg_enrich_bio$GO_ID, WvC_neg_enrich_bio$Adjusted.P.value), cbind(WvC_neg_enrich_bio$GO_ID, WvC_neg_enrich_bio$Adjusted.P.value))), "WvC_neg_GO.txt", delim = "\t", col_names = FALSE)
write_delim(as.data.frame(rbind(cbind(WvC_pos_enrich_bio$GO_ID, WvC_pos_enrich_bio$Adjusted.P.value), cbind(WvC_pos_enrich_bio$GO_ID, WvC_pos_enrich_bio$Adjusted.P.value), cbind(WvC_pos_enrich_bio$GO_ID, WvC_pos_enrich_bio$Adjusted.P.value))), "WvC_pos_GO.txt", delim = "\t", col_names = FALSE)
# CvH neg is empty
write_delim(as.data.frame(rbind(cbind(CvH_pos_enrich_bio$GO_ID, CvH_pos_enrich_bio$Adjusted.P.value), cbind(CvH_pos_enrich_bio$GO_ID, CvH_pos_enrich_bio$Adjusted.P.value), cbind(CvH_pos_enrich_bio$GO_ID, CvH_pos_enrich_bio$Adjusted.P.value))), "CvH_pos_GO.txt", delim = "\t", col_names = FALSE)
# Now two more text files for GO terms showing positive and negative DGE under CTmax conditions
# Negatively expressed GO terms are empty
write_delim(as.data.frame(cbind(pos_enrich_bio$GO_ID, pos_enrich_bio$Adjusted.P.value)), "CTmax_pos_GO.txt", delim = "\t", col_names = FALSE)
# Write results for CTmax-unique GO terms, positive or negative
write_delim(as.data.frame(rbind(cbind(CTmax_enrich_bio$GO_ID, CTmax_enrich_bio$Adjusted.P.value), cbind(CTmax_enrich_mol$GO_ID, CTmax_enrich_mol$Adjusted.P.value), cbind(CTmax_enrich_comp$GO_ID, CTmax_enrich_comp$Adjusted.P.value))), "CTmax_unique_GO.txt", delim = "\t", col_names = FALSE)

# Write a reference table with the full CTmax-unique GO term information, positive+negtive
write_delim(as.data.frame(rbind(CTmax_enrich_bio, CTmax_enrich_mol, CTmax_enrich_comp)), "CTmax_unique_GO_full.txt", delim = "\t", col_names = FALSE)

# Writing out the whole GO term spreadsheets into their own files
write_delim(as.data.frame(rbind(WvH_neg_enrich_bio, WvH_neg_enrich_mol, WvH_neg_enrich_comp)), "WvH_neg_GO_full.txt", delim = "\t", col_names = TRUE)
write_delim(as.data.frame(rbind(WvC_neg_enrich_bio, WvC_neg_enrich_mol, WvC_neg_enrich_comp)), "WvC_neg_GO_full.txt", delim = "\t", col_names = TRUE)
write_delim(as.data.frame(rbind(WvC_pos_enrich_bio, WvC_pos_enrich_mol, WvC_pos_enrich_comp)), "WvC_pos_GO_full.txt", delim = "\t", col_names = TRUE)
#write_delim(CvH_neg_enrich_bio, "CvH_neg_GO_full.txt", delim = "\t", col_names = TRUE) # Empty
write_delim(as.data.frame(rbind(CvH_pos_enrich_bio, CvH_pos_enrich_mol, CvH_pos_enrich_comp)), "CvH_pos_GO_full.txt", delim = "\t", col_names = TRUE)
#write_delim(neg_enrich_bio, "CTmax_neg_GO_full.txt", delim = "\t", col_names = TRUE) # Empty
write_delim(pos_enrich_bio, "CTmax_pos_GO_full.txt", delim = "\t", col_names = TRUE)

write_delim(combined_data, "clusters_w_annotations.txt", delim = "\t")

#######################################################################################
# Repeat the process of using enrichR for the three pairwise comparisons, but now search for overall lists of genes (not filtered by positive or negative)
# Now genes showing negative DGE under Wild vs CTmax conditions
WvC_overall_genes <- dplyr::filter(combined_data, abs(WvC.logFC) > 0 & !is.na(gene_name))
WvC_overall_genes <- dplyr::distinct(WvC_overall_genes, gene_name)
WvC_overall_genes <- as.vector(WvC_overall_genes)
length(WvC_overall_genes$gene_name)

# Run enrichR, searching the biological process and molecular function databases
WvC_overall_enriched <- enrichr(WvC_overall_genes$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
WvC_overall_enrich_bio <- as_tibble(WvC_overall_enriched[["GO_Biological_Process_2018"]])
WvC_overall_enrich_bio <- dplyr::filter(WvC_overall_enrich_bio, Adjusted.P.value < 0.05)
print(WvC_overall_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvC_overall_enrich_bio <- split_go(WvC_overall_enrich_bio)

WvC_overall_enrich_mol <- as_tibble(WvC_overall_enriched[["GO_Molecular_Function_2018"]])
WvC_overall_enrich_mol <- dplyr::filter(WvC_overall_enrich_mol, Adjusted.P.value < 0.05)
print(WvC_overall_enrich_mol$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvC_overall_enrich_mol <- split_go(WvC_overall_enrich_mol)

WvC_overall_enrich_comp <- as_tibble(WvC_overall_enriched[["GO_Cellular_Component_2018"]])
WvC_overall_enrich_comp <- dplyr::filter(WvC_overall_enrich_comp, Adjusted.P.value < 0.05)
WvC_overall_enrich_comp <- split_go(WvC_overall_enrich_comp)

# Count the number of GO terms in each enrichment database
nrow(WvC_overall_enrich_bio)
nrow(WvC_overall_enrich_mol)
nrow(WvC_overall_enrich_comp)

# Genes showing negative DGE under Wild vs Handling conditions
WvH_overall_genes <- dplyr::filter(combined_data, abs(WvH.logFC) > 0 & !is.na(gene_name))
WvH_overall_genes <- dplyr::distinct(WvH_overall_genes, gene_name)
WvH_overall_genes <- as.vector(WvH_overall_genes)
length(WvH_overall_genes$gene_name)

# Run enrichR, searching the biological process and molecular function databases
WvH_overall_enriched <- enrichr(WvH_overall_genes$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
WvH_overall_enrich_bio <- as_tibble(WvH_overall_enriched[["GO_Biological_Process_2018"]])
WvH_overall_enrich_bio <- dplyr::filter(WvH_overall_enrich_bio, Adjusted.P.value < 0.05)
print(WvH_overall_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvH_overall_enrich_bio <- split_go(WvH_overall_enrich_bio)

WvH_overall_enrich_mol <- as_tibble(WvH_overall_enriched[["GO_Molecular_Function_2018"]])
WvH_overall_enrich_mol <- dplyr::filter(WvH_overall_enrich_mol, Adjusted.P.value < 0.05)
print(WvH_overall_enrich_mol$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
WvH_overall_enrich_mol <- split_go(WvH_overall_enrich_mol)

WvH_overall_enrich_comp <- as_tibble(WvH_overall_enriched[["GO_Cellular_Component_2018"]])
WvH_overall_enrich_comp <- dplyr::filter(WvH_overall_enrich_comp, Adjusted.P.value < 0.05)
WvH_overall_enrich_comp <- split_go(WvH_overall_enrich_comp)

# Count the number of GO terms in each enrichment database
nrow(WvH_overall_enrich_bio)
nrow(WvH_overall_enrich_mol)
nrow(WvH_overall_enrich_comp)

# Genes showing negative DGE under CTmax vs Handling conditions
CvH_overall_genes <- dplyr::filter(combined_data, abs(CvH.logFC) > 0 & !is.na(gene_name))
CvH_overall_genes <- dplyr::distinct(CvH_overall_genes, gene_name)
CvH_overall_genes <- as.vector(CvH_overall_genes)
length(CvH_overall_genes$gene_name)

# Run enrichR, searching the biological process and molecular function databases
CvH_overall_enriched <- enrichr(CvH_overall_genes$gene_name, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
CvH_overall_enrich_bio <- as_tibble(CvH_overall_enriched[["GO_Biological_Process_2018"]])
CvH_overall_enrich_bio <- dplyr::filter(CvH_overall_enrich_bio, Adjusted.P.value < 0.05)
print(CvH_overall_enrich_bio$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
CvH_overall_enrich_bio <- split_go(CvH_overall_enrich_bio)

CvH_overall_enrich_mol <- as_tibble(CvH_overall_enriched[["GO_Molecular_Function_2018"]])
CvH_overall_enrich_mol <- dplyr::filter(CvH_overall_enrich_mol, Adjusted.P.value < 0.05)
print(CvH_overall_enrich_mol$Term)
# Use my custom function to split go terms into definitions and the formal term, useful for inputting data into revigo
CvH_overall_enrich_mol <- split_go(CvH_overall_enrich_mol)

CvH_overall_enrich_comp <- as_tibble(CvH_overall_enriched[["GO_Cellular_Component_2018"]])
CvH_overall_enrich_comp <- dplyr::filter(CvH_overall_enrich_comp, Adjusted.P.value < 0.05)
CvH_overall_enrich_comp <- split_go(CvH_overall_enrich_comp)

# Count the number of GO terms in each enrichment database
nrow(CvH_overall_enrich_bio)
nrow(CvH_overall_enrich_mol)
nrow(CvH_overall_enrich_comp)


# Write out the overall pairwise comparison results
write_delim(as.data.frame(rbind(WvC_overall_enrich_bio, WvC_overall_enrich_mol, WvC_overall_enrich_comp)), "WvC_overall_GO_full.txt", delim = "\t", col_names = TRUE)
write_delim(as.data.frame(rbind(WvH_overall_enrich_bio, WvH_overall_enrich_mol, WvH_overall_enrich_comp)), "WvH_overall_GO_full.txt", delim = "\t", col_names = TRUE)
write_delim(as.data.frame(rbind(CvH_overall_enrich_bio, CvH_overall_enrich_mol, CvH_overall_enrich_comp)), "CvH_overall_GO_full.txt", delim = "\t", col_names = TRUE)

################################################################################

# Check the numbers of the numbers of clusters in different situations, for the manuscript table
# Number of distinct Wild v CTmax clusters
nrow(filter(combined_data, abs(WvC.logFC) != 0) %>% distinct(corset_cluster))
# Distinct Wild v CTmax clusters with annotations
nrow(filter(combined_data, abs(WvC.logFC) != 0 & !is.na(gene_name)) %>% distinct(corset_cluster))
# Distinct Wild v CTmax clusters higher in Wild
nrow(filter(combined_data, WvC.logFC > 0) %>% distinct(corset_cluster))
# Distinct Wild v CTmax clusters higher in Wild with annotations
nrow(filter(combined_data, WvC.logFC > 0 & !is.na(gene_name)) %>% distinct(corset_cluster))
# Distinct Wild v CTmax clusters higher in CTmax 
nrow(filter(combined_data, WvC.logFC < 0) %>% distinct(corset_cluster))
# Distinct Wild v CTmax clusters higher in CTmax with annotations
nrow(filter(combined_data, WvC.logFC < 0 & !is.na(gene_name)) %>% distinct(corset_cluster))


# Follow the same process for Wild v Handling
nrow(filter(combined_data, abs(WvH.logFC) != 0) %>% distinct(corset_cluster))
# Distinct Wild v Handling clusters with annotations
nrow(filter(combined_data, abs(WvH.logFC) != 0 & !is.na(gene_name)) %>% distinct(corset_cluster))
# Distinct Wild v Handling clusters higher in Wild
nrow(filter(combined_data, WvH.logFC > 0) %>% distinct(corset_cluster))
# Distinct Wild v Handling clusters higher in Wild with annotations
nrow(filter(combined_data, WvH.logFC > 0 & !is.na(gene_name)) %>% distinct(corset_cluster))
# Distinct Wild v Handling clusters higher in Handling 
nrow(filter(combined_data, WvH.logFC < 0) %>% distinct(corset_cluster))
# Distinct Wild v Handling clusters higher in Handling with annotations
nrow(filter(combined_data, WvH.logFC < 0 & !is.na(gene_name)) %>% distinct(corset_cluster))


# Follow the same process for CTmax v Handling
nrow(filter(combined_data, abs(CvH.logFC) != 0) %>% distinct(corset_cluster))
# Distinct CTmax v Handling clusters with annotations
nrow(filter(combined_data, abs(CvH.logFC) != 0 & !is.na(gene_name)) %>% distinct(corset_cluster))
# Distinct CTmax v Handling clusters higher in CTmax
nrow(filter(combined_data, CvH.logFC > 0) %>% distinct(corset_cluster))
# Distinct CTmax v Handling clusters higher in CTmax with annotations
nrow(filter(combined_data, CvH.logFC > 0 & !is.na(gene_name)) %>% distinct(corset_cluster))
# Distinct CTmax v Handling clusters higher in Handling 
nrow(filter(combined_data, CvH.logFC < 0) %>% distinct(corset_cluster))
# Distinct CTmax v Handling clusters higher in Handling with annotations
nrow(filter(combined_data, CvH.logFC < 0 & !is.na(gene_name)) %>% distinct(corset_cluster))

# The same process for CTmax-unique genes
# Positive CTmax genes
nrow(filter(combined_data, pos_CTmax == 1 & is.na(neg_CTmax)) %>% distinct(corset_cluster))
# Positive CTmax genes with annotations
nrow(filter(combined_data, pos_CTmax == 1 & is.na(neg_CTmax) & !is.na(gene_name)) %>% distinct(corset_cluster))

# Negative CTmax genes
nrow(filter(combined_data, neg_CTmax == 1 & is.na(pos_CTmax)) %>% distinct(corset_cluster))
# Negative CTmax genes with annotations
nrow(filter(combined_data, neg_CTmax == 1 & is.na(pos_CTmax) & !is.na(gene_name)) %>% distinct(corset_cluster))

# Overall CTmax genes
nrow(filter(combined_data, neg_CTmax == 1 | pos_CTmax == 1) %>% distinct(corset_cluster))
# Overall CTmax genes with annotations
nrow(filter(combined_data, neg_CTmax == 1 | pos_CTmax == 1 & !is.na(gene_name)) %>% distinct(corset_cluster))
