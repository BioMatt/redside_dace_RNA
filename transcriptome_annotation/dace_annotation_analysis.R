# THis script is to take the output from Trinotate (Filtered for E<1e-3), and create a usable spreadsheet of transcript function for use with differential gene expression information

library(tidyverse)
options(stringsAsFactors = FALSE)

tidyverse_update()

trinotate <- read_delim("dace_annotation_report_filtered.txt", delim = "\t", col_names = TRUE)

# Filter out transcripts for which no blasted annotation information exists
trinotate <- dplyr::filter(trinotate, sprot_Top_BLASTX_hit != "." | sprot_Top_BLASTP_hit != ".")

# Removing the RNAmmer, transcript, and peptide columns, since they do not carry any information
trinotate <- dplyr::select(trinotate, -RNAMMER, -transcript, -peptide)

# Removing eggnog and Kegg terms, then removing gene ontology terms. While interesting, they each will not be used downstream. EnrichR can be used to find GO terms from sprot gene IDs.
trinotate <- dplyr::select(trinotate, -eggnog, -Kegg, -gene_ontology_BLASTP, -gene_ontology_Pfam)

# Checking how many genes have multiple transcripts
trinotate %>% 
  group_by('#gene_id', transcript_id) %>% 
  tally(sort = TRUE)


# Separate th sprot top blastx hit column into various other columns ID1, ID2, Q, perctID, and so on. Separated by the ^
# While some transcripts have many, many annotations associated with them, they are sorted with the lowest E-values reported first. Therefore splitting the blastx string and keeping only the first one will allow us to keep the annotation we are most confident in.
blastx <- trinotate %>% 
  separate(sprot_Top_BLASTX_hit, c("sp.BX_UniProt.ID.1", "sp.BX_UniProt.ID.2", "sp.BX_Q.H", "sp.BX_perctID", 
                                   "sp.BX_E", "sp.BX_RecName", "sp.BX_Tax.Lineage"), sep = "\\^") %>% 
  select(-sp.BX_UniProt.ID.2)  %>% # remove column since same as ...ID.1
  dplyr::rename(sp.BX_UniProt.ID = sp.BX_UniProt.ID.1) %>%  # rename UniProt_ID column
  separate(sp.BX_E, c("E", "sp.BX_E"), sep = "\\:") %>% # Separating the E's apart from E values to we can sort by quality
  select(-E)
blastx$sp.BX_E <- as.numeric(blastx$sp.BX_E)

# Cleaning up the transcript names and percent ID columns by separating the useful information from the blastx headers
blastx <- blastx %>% 
  separate(sp.BX_RecName, c("recname", "sp.BX_transcript_name"), sep = "=") %>% 
  select(-recname)
# Removing the semicolons from transcript descriptions just because they don't look good
blastx$sp.BX_transcript_name <- str_replace(blastx$sp.BX_transcript_name, ";", "")


blastx <- blastx %>% 
  separate(sp.BX_perctID, c("sp.BX_perctID", "ID"), sep = "%") %>% 
  select(-ID)
blastx$sp.BX_perctID <- as.numeric(blastx$sp.BX_perctID)


# Reading in the original Blastx output to capture bit scores, which Trinotate does not return
bit_scores <- read_delim("blastx.outfmt6", delim = "\t", col_names = FALSE)
bit_scores <- dplyr::select(bit_scores, transcript_id = X1, sp.BX_perctID = X3, sp.BX_E = X11, sp.BX_bit_score = X12)

# Add the bit scores to the annotation spreadsheet (they'll appear as the last column), joining by transcript ID, percent ID, and E value, of which all 3 must match before a bit score is attached to a transcript.
blastx <- left_join(blastx, bit_scores)
# Remove the bit score spreadsheet to keep the global environment clean
rm(bit_scores)
# Following Pearson 2013 (doi: 10.1002/0471250953.bi0301s42), percent ID is not being used as a filtering threshold. Instead, E values < 1e-6 and bit scores greater than 50 are being used.
# Out of 333,374 entries so far, filtering with these thresholds will leave 304,471 entries.
dplyr::filter(blastx, sp.BX_E < 1e-6 & sp.BX_bit_score > 50) %>% 
  tally()
blastx <- dplyr::filter(blastx, sp.BX_E < 1e-6 & sp.BX_bit_score > 50)

# Like what we did with blastx, we separate the blastp hit into different columns with informative column names
blastp <- blastx %>% 
  separate(sprot_Top_BLASTP_hit, c("sp.BP_UniProt.ID.1", "sp.BP_UniProt.ID.2", "sp.BP_Q.H", "sp.BP_perctID", 
                                   "sp.BP_E", "sp.BP_RecName", "sp.BP_Tax.Lineage"), sep = "\\^") %>% 
  select(-sp.BP_UniProt.ID.2)  %>% # remove column since same as ...ID.1
  dplyr::rename(sp.BP_UniProt.ID = sp.BP_UniProt.ID.1) %>%  # rename UniProt_ID column
  separate(sp.BP_E, c("E", "sp.BP_E"), sep = "\\:") %>% # Separating the E's apart from E values to we can sort by quality
  select(-E)

# Cleaning up the transcript names and percent ID columns by separating the useful information from the blastx headers
blastp <- blastp %>% 
  separate(sp.BP_RecName, c("recname", "sp.BP_transcript_name"), sep = "=") %>% 
  select(-recname)
blastp <- blastp %>% 
  separate(sp.BP_perctID, c("sp.BP_perctID", "ID"), sep = "%") %>% 
  select(-ID)
blastp$sp.BP_perctID <- as.numeric(blastp$sp.BP_perctID)
blastp$sp.BP_E <- as.numeric(blastp$sp.BP_E)

# Removing the semicolons from transcript descriptions just because they don't look good
blastp$sp.BP_transcript_name <- str_replace(blastp$sp.BP_transcript_name, ";", "")

# Now that the blastp annotations are cleaned up a bit, we will not filter them because they will act as supporting information for the blastx annotations, which were already filtered.

# Cleaning up the Pfam annotations
pfam <- blastp %>% 
  separate(Pfam, c("Pfam_ID", "Pfam_Protein", "Pfam_description", "Pfam_region", "Pfam_E"), sep = "\\^") %>% 
  select(-Pfam_region) %>% 
  separate(Pfam_E, c("E", "Pfam_E"), sep = "\\:") %>% # Separating the E's apart from E values to we can sort by quality
  select(-E)

# Cleaning up the signalP signal peptide predictions. Only keeping whether we believe the molecule is a signal peptide or not
signalp <- pfam %>% 
  separate(SignalP, c("signalP_ID", "signalP_score", "signalP_pct", "SignalP_Prediction"), sep = "\\^") %>% 
  select(-signalP_ID, -signalP_score, -signalP_pct)

# We are only interested in the transmembrane helices of proteins insofar as they support or refute evidence of other annotations. So tmhmm will be converted to a yes or no, like SignalP
tmhmm <- signalp %>% 
  mutate(TmHMM = if_else(TmHMM != ".", "YES", NULL))

# Writing the uniprot gene IDs to the uniprot retrieve ID/mapping tool for conversion to gene names, which are much more likely to be found in the EnrichR database down the line
write_delim(as.data.frame(unique(tmhmm$sp.BX_UniProt.ID)), "uniprot_IDs.txt", col_names = FALSE, delim = "\t")

# Read in the uniprot gene names with gene IDs from above, for attachment to the annotation spreadsheet
uniprot_gene_names <- read_delim("uniprot_ID_2_Name.txt", delim = "\t", col_names = TRUE)
# Attach the Uniprot gene names to the annotation spreadsheet
tmhmm <- left_join(tmhmm, uniprot_gene_names)
rm(uniprot_gene_names)

# Adding a supporting information column, carrying info from blastp, pfam, signalP, and tmhmm. This is intended to support the Blastx gene IDs, gene names, and transcript descriptions
tmhmm$Supporting_Info <- paste0("BlastP.TranscriptName:", tmhmm$sp.BP_transcript_name, ", Pfam.Protein:", tmhmm$Pfam_Protein, ", Pfam.Description:", tmhmm$Pfam_description, ", SignalPepdite?:", tmhmm$SignalP_Prediction, ", transmembrane?:", tmhmm$TmHMM) 

# Creating a new combined column of transcript ID and uniprot ID to look at how many duplicates exist
tmhmm <- mutate(tmhmm, combined = paste0(tmhmm$transcript_id, ":", tmhmm$sp.BX_UniProt.ID))
# 155,547 transcript and uniprot ID combinations are unique, out of 304,754 total
# It appears that different protein predictions from blastp have duplicated many transcripts.
length(unique(tmhmm$combined))

# Taking advantage of grouping within tidyverse, and grouping by transcript ID and uniprot ID combinations
# After grouping, ranking blastp protein e values from lowest to highest, and taking the lowest. If they are equal, taking the first row. This leaves us with the best protein information to go along with the transcript information
tmhmm <- tmhmm %>% 
  group_by(combined) %>% 
  filter(rank(sp.BP_E, ties.method="first")==1)

# Pulling only Gene IDs (blastx), gene names, descriptions, and supporting information, in conjunction with trinity transcript and gene IDs. This is for combining with differential gene expression analysis. In addition, renaming the columns to more descriptive names in more consistent formats.
annotations <- tibble(gene_id = tmhmm$'#gene_id', transcript_id = tmhmm$transcript_id, uniprot_id = tmhmm$sp.BX_UniProt.ID, gene_name = tmhmm$Gene_Name, transcript_description = tmhmm$sp.BX_transcript_name, supporting_info = tmhmm$Supporting_Info)

# Write out the file. 
write_delim(annotations, "dace_transcriptome_annotations.txt", delim = "\t", col_names = TRUE)
