library(tidyverse)
library(extrafont)
#font_import()
loadfonts(device = "win")
library(egg)
library(scales)
library(cowplot)
library(grid)


posCTmax_DGE <- read_delim("posCTmaxDGE_REVIGO.txt", delim = "\t") %>% 
  filter(eliminated == 0) %>% 
  arrange(desc(Enrichment_Database), Number_of_Genes)
posCTmax_DGE$term_ID <- factor(posCTmax_DGE$term_ID, levels = posCTmax_DGE$term_ID)



posCTmax_DGE_plot <- ggplot(data = posCTmax_DGE, aes(x = term_ID, y = Number_of_Genes, colour = Enrichment_Database, fill = Enrichment_Database)) +
  geom_col() +
  scale_x_discrete(labels = posCTmax_DGE$description) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  coord_flip() +
  theme_classic() +
  xlab(element_blank()) +
  ylab("Number of Clusters") +
  scale_colour_manual(name = "Enrichment Database", labels = c("Biological Process", "Molecular Function"), values = c("#4575b4", "#fdae61"), aesthetics = c("colour", "fill")) +
theme(text=element_text(size=16,  family="Arial"), axis.ticks = element_blank(), legend.position = "none")
posCTmax_DGE_plot

posCTmax_DEU <- read_delim("posCTmaxDEU_REVIGO.txt", delim = "\t") %>% 
  filter(eliminated == 0) %>% 
  arrange(factor(Enrichment_Database, levels = c("Cellular Component 2018", "Molecular Function 2018", "Biological Process 2018")), Number_of_Genes)
posCTmax_DEU$term_ID <- factor(posCTmax_DEU$term_ID, levels = posCTmax_DEU$term_ID)


posCTmax_DEU_plot <- ggplot(data = posCTmax_DEU, aes(x = term_ID, y = Number_of_Genes, colour = Enrichment_Database, fill = Enrichment_Database, group = Enrichment_Database)) +
  geom_col() +
  scale_x_discrete(labels = posCTmax_DEU$description) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme_classic() +
  xlab(element_blank()) +
  ylab(element_blank()) +
  scale_colour_manual(name = "Enrichment Database", labels = c("Biological Process 2018" = "Biological Process", "Molecular Function 2018" = "Molecular Function", "Cellular Component 2018" = "Cellular Component"), values = c("#4575b4", "#fdae61", "#d73027"), aesthetics = c("colour", "fill"), breaks = c("Biological Process 2018", "Molecular Function 2018", "Cellular Component 2018")) +
  theme(text=element_text(size=16,  family="Arial"), axis.ticks = element_blank()) +
  coord_flip()
posCTmax_DEU_plot

negCTmax_DEU <- read_delim("negCTmaxDEU_REVIGO.txt", delim = "\t") %>% 
  arrange(desc(Enrichment_Database), Number_of_Genes)
negCTmax_DEU$GO_ID <- factor(negCTmax_DEU$GO_ID, levels = negCTmax_DEU$GO_ID)

neg_CTmax_DEU_plot <- ggplot(data = negCTmax_DEU, aes(x = GO_ID, y = Number_of_Genes, colour = Enrichment_Database, fill = Enrichment_Database)) +
  geom_col() +
  scale_x_discrete(labels = negCTmax_DEU$GO_term) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  coord_flip() +
  theme_classic() +
  xlab(element_blank()) +
  ylab("Number of Clusters") +
  scale_colour_manual(name = "Enrichment Database", labels = c("Biological Process", "Molecular Function"), values = c("#4575b4", "#fdae61"), aesthetics = c("colour", "fill")) +
  theme(text=element_text(size=16,  family="Arial"), axis.ticks = element_blank(), legend.position = "none")
neg_CTmax_DEU_plot

combined_GO_terms <- posCTmax_DGE_plot | (posCTmax_DEU_plot / neg_CTmax_DEU_plot)
combined_GO_terms

posCTmax_DEU_plot / neg_CTmax_DEU_plot



# Trying the egg package to change plot heights. Using a ratio of 6/45 for the negative DEU plot because it has 6 GO terms to the positive DEU's 45 GO terms
ggarrange(posCTmax_DEU_plot, neg_CTmax_DEU_plot, ncol = 1, heights = c(1, (6/45)), widths = 0.5)

# The above plot and the one below are exported at 1500 width x 817 as SVG files for editing with Inkscape. Plot titles were added, arrows, and the legend was moved
posCTmax_DGE_plot

#######################################################################################################
# Make a plot of all DEU genes, regardless of positive or negative
# Read in the individual enrichment database results, pre-Revigo, to get gene counts per GO term
overall_bio_DEU <- read_tsv("overall_enriched_bioGO.txt") %>% 
  separate(Overlap, into = c("gene_count", "total_genes"), sep = "/")

overall_mol_DEU <- read_tsv("overall_enriched_molGO.txt") %>% 
  separate(Overlap, into = c("gene_count", "total_genes"), sep = "/")

overall_comp_DEU <- read_tsv("overall_enriched_cellGO.txt") %>% 
  separate(Overlap, into = c("gene_count", "total_genes"), sep = "/")

overall_GO_noRevigo <- rbind(overall_bio_DEU, overall_mol_DEU, overall_comp_DEU)

# Read in the overall Revigo results
DEU_overall_Revigo <- read_delim("overall_Revigo_DEU.txt", delim = "\t") %>% 
  filter(eliminated == 0) %>% 
  dplyr::rename(GO_ID = term_ID)

# Attach the original GO term and the Revigo-filtered GO terms together, to get both gene counts and non-redundant terms
DEU_overall_Revigo <- left_join(DEU_overall_Revigo, overall_GO_noRevigo) %>% 
  select(-GO_term)

# Re-format columns as factors or numeric. Re-order for a nice plot
DEU_overall_Revigo$Enrichment_Database <- factor(DEU_overall_Revigo$Enrichment_Database, levels = c("Biological Process 2018", "Molecular Function 2018", "Cellular Component 2018"))
DEU_overall_Revigo$gene_count <- as.numeric(DEU_overall_Revigo$gene_count)

DEU_overall_Revigo <- DEU_overall_Revigo %>% 
  arrange(desc(Enrichment_Database), gene_count)


DEU_overall_Revigo$GO_ID <- factor(DEU_overall_Revigo$GO_ID, levels = DEU_overall_Revigo$GO_ID)
DEU_overall_Revigo$description <- factor(DEU_overall_Revigo$description, levels = DEU_overall_Revigo$description)

# Plot the overall DEU results
DEU_overall_plot <- ggplot(data = DEU_overall_Revigo, aes(x = GO_ID, y = gene_count, colour = Enrichment_Database, fill = Enrichment_Database)) +
  geom_col() +
  scale_x_discrete(labels = DEU_overall_Revigo$description) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  coord_flip() +
  theme_classic() +
  xlab(element_blank()) +
  ylab("Number of Clusters") +
  scale_colour_manual(name = "Enrichment Database", labels = c("Biological Process", "Molecular Function", "Cellular Component"), values = c("#4575b4", "#fdae61", "#d73027"), aesthetics = c("colour", "fill"), breaks = c("Biological Process 2018", "Molecular Function 2018", "Cellular Component 2018")) +
  theme(text=element_text(size=16,  family = "Arial"), axis.ticks = element_blank(), legend.position = c(0.75, 0.075))

ggsave(filename = "DEU_overall_GO.pdf", plot = DEU_overall_plot, dpi = 2000, device = cairo_pdf, height = 10, width = 12)

###########################################################################################################
# Create a plot for CTmax-unique GO terms, positive or negative
# Read in the data, remove extraneous columns, filter for non-redundant GO-terms, add enrichment database info
CTmax_bio <- read_tsv("CTmax_bio_unique_Revigo.txt") %>% 
  select(-`...12`) %>% 
  filter(Eliminated == FALSE) %>% 
  mutate(enrichment_database = "biological_process")

CTmax_mol <- read_tsv("CTmax_mol_unique_Revigo.txt") %>% 
  filter(Eliminated == FALSE) %>% 
  mutate(enrichment_database = "molecular_function")

CTmax_cell <- read_tsv("CTmax_comp_unique_Revigo.txt") %>% 
  filter(Eliminated == FALSE) %>% 
  mutate(enrichment_database = "cellular_component")

# Read in the enrichR information for the CTmax GO terms, get a table of the number of genes for each GO term
CTmax_genes <- read_tsv("CTmax_unique_GO_full.txt", col_names = FALSE) %>% 
  select(X2, X3)
colnames(CTmax_genes) <- c("TermID", "fraction")
CTmax_counts <- stringr::str_split_fixed(CTmax_genes$fraction, "/", n = 2)
CTmax_counts <- as.data.frame(CTmax_counts)
colnames(CTmax_counts) <- c("gene_count", "total_genes")
CTmax_genes <- cbind(CTmax_genes, CTmax_counts)

# Combine the datasets
CTmax_unique_terms <- rbind(CTmax_bio, CTmax_mol, CTmax_cell)
# Match up the GO terms from Revigo with the number of genes in each term from enrichR
CTmax_unique_terms <- left_join(CTmax_unique_terms, CTmax_genes)

# Use factors and numeric variables to re-order by first enrichment database, then number of genes
CTmax_unique_terms$enrichment_database <- factor(CTmax_unique_terms$enrichment_database, levels = c("biological_process", "molecular_function", "cellular_component"))
CTmax_unique_terms$gene_count <- as.numeric(CTmax_unique_terms$gene_count)

# Re-order the terms by enrichment database and number of genes
CTmax_unique_terms <- CTmax_unique_terms %>% 
  arrange(desc(enrichment_database), gene_count)

CTmax_unique_terms$TermID <- factor(CTmax_unique_terms$TermID, levels = c(CTmax_unique_terms$TermID))
CTmax_unique_terms$Name <- factor(CTmax_unique_terms$Name, levels = c(CTmax_unique_terms$Name))
CTmax_unique_terms$enrichment_database <- factor(CTmax_unique_terms$enrichment_database, levels = c("biological_process", "molecular_function", "cellular_component"))


CTmax_unique_plot <- ggplot(data = CTmax_unique_terms, aes(x = TermID, y = gene_count, colour = enrichment_database, fill = enrichment_database)) +
  geom_col() +
  scale_x_discrete(labels = CTmax_unique_terms$Name) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  coord_flip() +
  theme_classic() +
  xlab(element_blank()) +
  ylab("Number of Clusters") +
  scale_colour_manual(name = "Enrichment Database", labels = c("Biological Process", "Molecular Function", "Cellular Component"), values = c("#4575b4", "#fdae61", "#d73027"), aesthetics = c("colour", "fill"), breaks = c("biological_process", "molecular_function", "cellular_component")) +
  theme(text=element_text(size=16,  family="Arial"), axis.ticks = element_blank(), legend.position = c(0.75, 0.090))
#CTmax_unique_plot

ggsave(filename = "CTmax_unique_GO.pdf", plot = CTmax_unique_plot, dpi = 2000, width = 12, height = 9, device = cairo_pdf)

# Get the legend
legend_plot <- ggplot(data = CTmax_unique_terms, aes(x = TermID, y = gene_count, colour = enrichment_database, fill = enrichment_database)) +
                                geom_col() +
                                scale_x_discrete(labels = CTmax_unique_terms$Name) +
                                scale_y_continuous(expand = expansion(mult = c(0, .1))) +
                                coord_flip() +
                                theme_classic() +
                                xlab(element_blank()) +
                                ylab("Number of Clusters") +
                                scale_colour_manual(name = "Enrichment Database", labels = c("Biological Process", "Molecular Function", "Cellular Component"), values = c("#4575b4", "#fdae61", "#d73027"), aesthetics = c("colour", "fill"), breaks = c("biological_process", "molecular_function", "cellular_component")) +
                                theme(text=element_text(size=16,  family="Arial"), axis.ticks = element_blank(), legend.position = "right")

ggsave(filename = "CTmax_unique_GO_legend.pdf", plot = legend_plot, dpi = 2000, device = cairo_pdf)
