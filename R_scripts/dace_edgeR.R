# This script reads in Corset-clustered reads for redside dace and calculates differential gene expression (DGE) using edgeR
# Taken from the edgeR user's manual: https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pheatmap")

library(tidyverse)
library(edgeR)
library(pheatmap)
library(ggfortify)
library(ggrepel)
library(extrafont)
#font_import()
loadfonts()
#extrafont::font_import()
library(patchwork)
library(scales)

# Reading the corset gene counts
corset_counts <- read_tsv("dace-counts.txt", col_names = TRUE)
corset_counts <- column_to_rownames(corset_counts, var = "gene")

# Read the experiment metadata
metadata <- read_tsv("dace_metadata.txt", col_names = TRUE)

data_edger <- DGEList(counts = corset_counts, group = factor(metadata$treatment, levels = c("CT_max", "Handle", "Wild")))

# Filter the count data for genes with any expression
keep <- filterByExpr(data_edger)
table(keep)
data_edger <- data_edger[keep,]

# Normalize the data for library size
data_edger <- calcNormFactors(data_edger)
data_edger$samples

# Take a quick initial look at the data
plotMDS(data_edger, col = rep(1:3, each = 10), top = nrow(data_edger))
mds <- plotMDS(data_edger, col = rep(1:3, each = 10), top = nrow(data_edger), plot = FALSE)

# Create my design matrix for running the EdgeR model
treatment <- factor(metadata$treatment, levels = c("CT_max", "Handle", "Wild"))
RIN <- metadata$RIN

design <- model.matrix(~ 0 + treatment + RIN)
no_RIN <- model.matrix(~ 0 + treatment)

# Estimate dispersion in the data
data_edger <- estimateDisp(data_edger, design, robust = TRUE)
plotBCV(data_edger)

# Taking a look at the effect of RIN score on the model fit
voom(data_edger, design = design, plot = TRUE)
voom(data_edger, design = no_RIN, plot = TRUE)


# Making some contrasts between the three treatments
my.contrasts <- makeContrasts(Wild.v.Handle = treatmentWild - treatmentHandle,
                              Wild.v.CTmax = treatmentWild - treatmentCT_max,
                              CTmax.v.Handle = treatmentCT_max - treatmentHandle,
                              levels = design)


# Use a GLM with quasi-likelihoods to test for differential expression
quasi_glm_fit <- glmQLFit(data_edger, design, robust = TRUE)
plotQLDisp(quasi_glm_fit)

# Pulling results for Wild versus Handled fish
quasi_glm_test_WvH <- glmQLFTest(quasi_glm_fit, contrast = my.contrasts[,"Wild.v.Handle"])
topTags(quasi_glm_test_WvH, p.value = 0.05)
plotMD(quasi_glm_test_WvH)
WvH_res <- as.data.frame(topTags(quasi_glm_test_WvH, n = Inf, p.value = 0.05))
WvH_res <- rownames_to_column(WvH_res, var = "corset_cluster")
nrow(filter(WvH_res, logFC > 0))
nrow(filter(WvH_res, logFC < 0))
#WvH_res <- dplyr::filter(WvH_res, abs(logFC) > 2)

# Write out the data
write_delim(WvH_res, "dace_WvH_results.txt", delim = "\t", col_names = TRUE)

# Pulling results for Wild versus CTmax-treated fish
quasi_glm_test_WvC <- glmQLFTest(quasi_glm_fit, contrast = my.contrasts[,"Wild.v.CTmax"])
topTags(quasi_glm_test_WvC, p.value = 0.05)
plotMD(quasi_glm_test_WvC)
WvC_res <- as.data.frame(topTags(quasi_glm_test_WvC, n = Inf, p.value = 0.05))
WvC_res <- rownames_to_column(WvC_res, var = "corset_cluster")
nrow(filter(WvC_res, logFC > 0))
nrow(filter(WvC_res, logFC < 0))
#WvC_res <- dplyr::filter(WvC_res, abs(logFC) > 2)

# Write out the data
write_delim(WvC_res, "dace_WvC_results.txt", delim = "\t", col_names = TRUE)

# Pulling results for CTmax-treated fish versus Wild fish
quasi_glm_test_CvH <- glmQLFTest(quasi_glm_fit, contrast = my.contrasts[,"CTmax.v.Handle"])
topTags(quasi_glm_test_CvH, p.value = 0.05)
plotMD(quasi_glm_test_CvH)
CvH_res <- as.data.frame(topTags(quasi_glm_test_CvH, n = Inf, p.value = 0.05))
CvH_res <- rownames_to_column(CvH_res, var = "corset_cluster")
nrow(filter(CvH_res, logFC > 0))
nrow(filter(CvH_res, logFC < 0))
#CvH_res <- dplyr::filter(CvH_res, abs(logFC) > 2)

# Write out the data
write_delim(CvH_res, "dace_CvH_results.txt", delim = "\t", col_names = TRUE)

# Now that the pairwise comparison results for all three treatments are done, we pull the transcripts that only show DGE in CTmax conditions, not Handling conditions. The idea being to subtract out those transcripts that reflect a general stress response and focus on a heat-specific response.
CTmax_only <- anti_join(CvH_res, WvH_res, by = "corset_cluster")
CTmax_only <- semi_join(CTmax_only, WvC_res, by = "corset_cluster")
# Write out the CTmax only data, only keeping cluster information
write_delim(dplyr::select(CTmax_only, corset_cluster), "dace_CTmax_unique_clusters.txt", delim = "\t", col_names = TRUE)

# I do not know how legitimate the method above was for isolating CTmax-specific transcripts.
# Instead, taking transcripts that show positive DGE under CTmax conditions only. So positive DGE with the CTmax versus handling comparison, and negative DGE with the wild versus CTmax comparison
pos_CTmax <- dplyr::filter(CvH_res, logFC > 0 & FDR < 0.05)
pos_CTmax <- semi_join(pos_CTmax, dplyr::filter(WvC_res, logFC < 0 & FDR < 0.05), by = "corset_cluster")
# Keep only cluster information because with 2 different DGE datasets represented by these clusters, we only know that LFC was positive but not to what magnitude
pos_CTmax <- dplyr::select(pos_CTmax, corset_cluster)
# Write out the data
write_delim(pos_CTmax, "dace_pos_CTmax_clusters.txt", delim = "\t", col_names = TRUE)

# Now looking at transcripts that show negative DGE under CTmax conditions only. So negative DGE with CTmax versus handling comparison, and positive DGE with the wild versus CTmax comparison.
neg_CTmax <- dplyr::filter(CvH_res, logFC < 0 & FDR < 0.05)
neg_CTmax <- semi_join(neg_CTmax, dplyr::filter(WvC_res, logFC > 0 & FDR < 0.05), by = "corset_cluster")
neg_CTmax <- dplyr::select(neg_CTmax, corset_cluster)
# Write out the data
write_delim(neg_CTmax, "dace_neg_CTmax_clusters.txt", delim = "\t", col_names = TRUE)


###########################################################################################
# Use pretty heatmaps (pheatmaps) to visualize the data
# Following the recommendations of the edgeR user's guide, section 2.15 (https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf), turning the data into counts per million for plotting as heat maps
# First, increasing the memory available to R
memory.limit(size=11000)
# Selecting the transcripts to be plotted in the heatmap in several ways. First, by the top 1000 in mean expression across all individuals
select <- order(rowSums(cpm_data), decreasing = TRUE)[1:1000]
# Now by the top 4000 in variance in expression across all individuals
select <- order(apply(cpm_data, MARGIN=1, FUN=var, na.rm=TRUE), decreasing = TRUE)[1:4000]
cpm_heatmap <- pheatmap(cpm_data[select[1:length(select)],], cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE, main = "Log CPM Heatmap")
ggsave("log_cpm_DGE_heatmap.pdf", plot = cpm_heatmap, device = "pdf", dpi = 2000)


# Selecting any differentially expressed gene
select <- unique(c(WvC_res$corset_cluster, WvH_res$corset_cluster, CvH_res$corset_cluster))

cpm_data <- cpm_data[select,]
cpm_data <- order(rowMeans(cpm_data), decreasing = TRUE)

# Selecting any differentially expressed gene took too much memory, so instead just making a separate heatmap for each comparison
pheatmap(cpm_data[WvC_res$corset_cluster,], cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE, main = "Wild v CTmax Heatmap")
pheatmap(cpm_data[WvH_res$corset_cluster,], cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE, main = "Wild v Handling Heatmap")
pheatmap(cpm_data[CvH_res$corset_cluster,], cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE, main = "CTmax v Handling Heatmap")

# Trying a heatmap with log2 fold changes to show differences in a more stark way
lfc_data <- rename(WvC_res, WvC_lfc = logFC)
lfc_data <- dplyr::select(lfc_data, corset_cluster, WvC_lfc)
lfc_data <- full_join(lfc_data, WvH_res)
lfc_data <- lfc_data %>% 
  rename(WvH_lfc = logFC) %>% 
  select(-logCPM, -F, -PValue, -FDR)
lfc_data <- full_join(lfc_data, CvH_res)
lfc_data <- lfc_data %>% 
  rename(CvH_lfc = logFC) %>% 
  select(-logCPM, -F, -PValue, -FDR)
lfc_data <- column_to_rownames(lfc_data, var = "corset_cluster")
lfc_data <- replace_na(lfc_data, list(WvC_lfc = 0, WvH_lfc = 0, CvH_lfc = 0))

lfc_data <- as.matrix(lfc_data)

select <- order(rowMeans(abs(lfc_data)), decreasing = TRUE)[1:200]


heatmap <- pheatmap(lfc_data[select,], cluster_cols = TRUE, cluster_rows = TRUE, show_rownames = TRUE, main = "log2-Fold Change Heatmap", fontsize = 1)
ggsave("lfc_heatmap_top200.pdf", plot = heatmap, device = "pdf", dpi = 4000)

# Run a PCA with prcomp
pca <- prcomp(t(cpm(data_edger, log = TRUE)), scale. = TRUE)
summary(pca)
autoplot(pca, colour = as.numeric(as.factor((metadata$treatment))), label = TRUE, shape = FALSE, loadings = TRUE)

###################################################################
# Making the MDS nicer
mds_plot <- tibble(ind = metadata$sample, treatment = metadata$treatment, x = mds$x, y = mds$y)

ggmds_plot <- ggplot(data = mds_plot, aes(x = x, y = y, colour = treatment, label = ind)) +
  geom_point(size = 2) +
  theme_classic() +
  geom_label_repel(aes(label = ind, family="Arial", size= (8*(5/14))), hjust = 0, vjust = 0.1) +
  scale_colour_manual(values = c("#d73027", "#fdae61", "#4575b4")) +
  xlab(expression("Log"[2]*"-Fold Change Dimension 1")) +
  ylab(expression("Log"[2]*"-Fold Change Dimension 2")) +
  theme(legend.position = "none", axis.ticks = element_blank(), axis.text = element_blank(), text=element_text(size=16,  family="Arial"))
ggmds_plot
ggsave("MDS_plot.png", plot = ggmds_plot, device = "png", dpi = 2000)

#########################################################
# Making a heatmap of clusters that showed only positive or negative DGE in CTmax (compared to both controls)
# Creating the heat map, removing the dendrogram for the clusters
# Re-order individuals to look nicer
other_ctmax_clusters <- rbind(pos_CTmax, neg_CTmax)
ctmax_heatmap2 <- pheatmap(cpm_data[other_ctmax_clusters$corset_cluster,c("Handle.1", "Handle.2", "Handle.3", "Handle.4", "Handle.5", "Handle.6", "Handle.7", "Handle.8", "Handle.9", "Handle.10", "CTmax.1", "CTmax.2", "CTmax.3", "CTmax.4", "CTmax.5", "CTmax.6", "CTmax.7", "CTmax.8", "CTmax.9", "CTmax.10", "Wild.1", "Wild.2", "Wild.3", "Wild.4", "Wild.5", "Wild.6", "Wild.7", "Wild.8", "Wild.9", "Wild.10")], cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE, treeheight_row = 0, clustering_method="complete", fontsize = 14, fontfamily = "Arial", gaps_col = c(10, 20), labels_col = c(1:10, 1:10, 1:10), angle_col = 0)
ggsave("log_cpm_CTmax-unique_DGE_heatmap.png", plot = ctmax_heatmap2, device = "png", dpi = 2000)


#################################################################
# Pulling counts per million for each cluster
cpm_data <- cpm(data_edger, log = FALSE)

# Plotting prpf38b
prpf38b <- tibble(individual = colnames(cpm_data), cpm = cpm_data["Cluster-41845.39916",], Treatment = c(rep("CTmax", 10), rep("Handle", 10), rep("Wild", 10)))

prpf38b_plot <- ggplot(data = prpf38b, aes(x = Treatment, y = cpm, group = Treatment, colour = Treatment)) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("#4575b4", "#fdae61", "#d73027")) +
  theme_classic() +
  xlab("Treatment") +
  ylab("Counts per Million") +
  theme(text=element_text(size=16,  family="Arial"))
ggsave("prpf38b_cpm.png", plot = prpf38b_plot, device = "png", dpi = 2000)

# Plotting JUN, transcription factor AP-1
jun <- tibble(individual = colnames(cpm_data), cpm = cpm_data["Cluster-71016.0",], Treatment = c(rep("CTmax", 10), rep("Handle", 10), rep("Wild", 10)))
jun_plot <- ggplot(data = jun, aes(x = Treatment, y = cpm, group = Treatment, colour = Treatment)) +
  geom_jitter(size = 4, alpha = 0.8, width = 0.1) +
  scale_colour_manual(values = c("#d73027", "#fdae61", "#4575b4")) +
  theme_classic() +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(text=element_text(size=16,  family="Arial"), legend.position = "none") +
  ggtitle(label = "JUN") +
  scale_y_continuous(labels = label_number())
jun_plot
#ggsave("jun_cpm.png", plot = jun_plot, device = "png", dpi = 2000)

junb <- tibble(individual = colnames(cpm_data), cpm = cpm_data["Cluster-63072.0",], Treatment = c(rep("CTmax", 10), rep("Handle", 10), rep("Wild", 10)))
junb_plot <- ggplot(data = junb, aes(x = Treatment, y = cpm, group = Treatment, colour = Treatment)) +
  geom_jitter(size = 4, alpha = 0.8, width = 0.1) +
  scale_colour_manual(values = c("#d73027", "#fdae61", "#4575b4")) +
  theme_classic() +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(text=element_text(size=16,  family="Arial"), legend.position = "none") +
  ggtitle("jun-B") +
  scale_y_continuous(labels = label_number())
junb_plot

jund <- tibble(individual = colnames(cpm_data), cpm = cpm_data["Cluster-41845.6375",], Treatment = c(rep("CTmax", 10), rep("Handle", 10), rep("Wild", 10)))
jund_plot <- ggplot(data = jund, aes(x = Treatment, y = cpm, group = Treatment, colour = Treatment)) +
  geom_jitter(size = 4, alpha = 0.8, width = 0.1) +
  scale_colour_manual(values = c("#d73027", "#fdae61", "#4575b4")) +
  theme_classic() +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(text=element_text(size=16,  family="Arial"), legend.position = "none") +
  ggtitle("jun-D") +
  scale_y_continuous(labels = label_number())
jund_plot


smad4 <- tibble(individual = colnames(cpm_data), cpm = cpm_data["Cluster-41845.619",], Treatment = c(rep("CTmax", 10), rep("Handle", 10), rep("Wild", 10)))
smad4_plot <- ggplot(data = smad4, aes(x = Treatment, y = cpm, group = Treatment, colour = Treatment)) +
  geom_jitter(size = 4, alpha = 0.8, width = 0.1) +
  scale_colour_manual(values = c("#d73027", "#fdae61", "#4575b4")) +
  theme_classic() +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(text=element_text(size=16,  family="Arial"), legend.position = "none") +
  ggtitle("SMAD4") +
  scale_y_continuous(labels = label_number())
smad4_plot



tgfbr <- tibble(individual = colnames(cpm_data), cpm = cpm_data["Cluster-41845.9758",], Treatment = c(rep("CTmax", 10), rep("Handle", 10), rep("Wild", 10)))
tgfbr_plot <- ggplot(data = tgfbr, aes(x = Treatment, y = cpm, group = Treatment, colour = Treatment)) +
  geom_jitter(size = 4, alpha = 0.8, width = 0.1) +
  scale_colour_manual(values = c("#d73027", "#fdae61", "#4575b4")) +
  theme_classic() +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(text=element_text(size=16,  family="Arial"), legend.position = "none") +
  ggtitle("TGFR1") +
  scale_y_continuous(labels = label_number())
tgfbr_plot
######################################################################
# Taking a look at many inducible transcription factors
# Plotting IER2
ier2 <- tibble(individual = colnames(cpm_data), cpm = cpm_data["Cluster-73275.0",], Treatment = c(rep("CTmax", 10), rep("Handle", 10), rep("Wild", 10)))
ier2_plot <- ggplot(data = ier2, aes(x = Treatment, y = cpm, group = Treatment, colour = Treatment)) +
  geom_jitter(size = 4, alpha = 0.8, width = 0.1) +
  scale_colour_manual(values = c("#d73027", "#fdae61", "#4575b4")) +
  theme_classic() +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(text=element_text(size=16,  family="Arial"), legend.position = "none") +
  ggtitle("IER2") +
  scale_y_continuous(labels = label_number())
ier2_plot

# Plotting c-FOS
cFOS <- tibble(individual = colnames(cpm_data), cpm = cpm_data["Cluster-83055.0",], Treatment = c(rep("CTmax", 10), rep("Handle", 10), rep("Wild", 10)))
cFOS_plot <- ggplot(data = cFOS, aes(x = Treatment, y = cpm, group = Treatment, colour = Treatment)) +
  geom_jitter(size = 4, alpha = 0.8, width = 0.1) +
  scale_colour_manual(values = c("#d73027", "#fdae61", "#4575b4")) +
  theme_classic() +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(text=element_text(size=16,  family="Arial"), legend.position = "none") +
  ggtitle("c-Fos") +
  scale_y_continuous(labels = label_number())
cFOS_plot

fosb <- tibble(individual = colnames(cpm_data), cpm = cpm_data["Cluster-21696.1",], Treatment = c(rep("CTmax", 10), rep("Handle", 10), rep("Wild", 10)))
fosb_plot <- ggplot(data = fosb, aes(x = Treatment, y = cpm, group = Treatment, colour = Treatment)) +
  geom_jitter(size = 4, alpha = 0.8, width = 0.1) +
  scale_colour_manual(values = c("#d73027", "#fdae61", "#4575b4")) +
  theme_classic() +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(text=element_text(size=16,  family="Arial"), legend.position = "none") +
  ggtitle("fosB") +
  scale_y_continuous(labels = label_number())
fosb_plot

# Plotting MYC
myc <- tibble(individual = colnames(cpm_data), cpm = cpm_data["Cluster-41845.21943",], Treatment = c(rep("CTmax", 10), rep("Handle", 10), rep("Wild", 10)))
myc_plot <- ggplot(data = myc, aes(x = Treatment, y = cpm, group = Treatment, colour = Treatment)) +
  geom_jitter(size = 4, alpha = 0.8, width = 0.1) +
  scale_colour_manual(values = c("#d73027", "#fdae61", "#4575b4")) +
  theme_classic() +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(text=element_text(size=16,  family="Arial"), legend.position = "none") +
  ggtitle("MYC") +
  scale_y_continuous(labels = label_number())
myc_plot

# Plotting SERPINH1
serpin <- tibble(individual = colnames(cpm_data), cpm = cpm_data["Cluster-73213.0",], Treatment = c(rep("CTmax", 10), rep("Handle", 10), rep("Wild", 10)))
serpin_plot <- ggplot(data = serpin, aes(x = Treatment, y = cpm, group = Treatment, colour = Treatment)) +
  geom_jitter(size = 4, alpha = 0.8, width = 0.1) +
  scale_colour_manual(values = c("#d73027", "#fdae61", "#4575b4")) +
  theme_classic() +
  xlab("Treatment") +
  theme(text=element_text(size=16,  family="Arial"), legend.position = "none") +
  ggtitle("SERPINH1") +
  scale_y_continuous(labels = label_number())
serpin_plot







acute_stress <- jun_plot /  (junb_plot | jund_plot) / (ier2_plot | myc_plot) / (cFOS_plot | fosb_plot)
acute_stress
ggsave("acute_stress.png", plot = acute_stress, device = "png", dpi = 2000, width = 5.5, height = 12)
