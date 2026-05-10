############################################################
# Project: Microbial and Genetic Diversity of Medaka Fish
#          Under Osmotic Stress
#
# Purpose:
# This script performs:
#   1. 16S amplicon processing using DADA2
#   2. ASV table generation
#   3. Taxonomy assignment using SILVA
#   4. Alpha diversity analysis
#   5. Microbial composition plots
#   6. Beta diversity analysis using Bray-Curtis PCoA
#   7. PERMANOVA statistical testing
#   8. COI-based phylogenetic tree construction
#
# Input files required in D:/Dixcy:
#   - *_1.fastq.gz and *_2.fastq.gz files
#   - silva_nr99_v138.1_train_set.fa.gz
#   - SraRunTable.csv
#   - COI_Oryzias_sequences.fasta
#
# Important:
# This workflow uses forward reads only for ASV generation because
# paired-end merging failed for some samples. This allowed all 10
# samples to be retained for downstream analysis.
############################################################


############################################################
# STEP 1: Install and load required R packages
#
# Why this step?
# R needs specific packages for microbiome analysis, plotting,
# statistics, and phylogenetic tree construction.
############################################################

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages
BiocManager::install(c("dada2", "phyloseq", "Biostrings", "msa"), ask = FALSE, update = FALSE)

# Install CRAN packages
install.packages(c("vegan", "ggplot2", "dplyr", "tidyr", "tibble", "readr", "ape", "seqinr"))

# Load packages
library(dada2)
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(ape)
library(seqinr)
library(Biostrings)
library(msa)


############################################################
# STEP 2: Set working directory and identify FASTQ files
#
# Why this step?
# We tell R where the project files are stored and identify
# forward and reverse sequencing files.
############################################################

path <- "D:/Dixcy"
setwd(path)

# Forward reads end with _1.fastq.gz
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))

# Reverse reads end with _2.fastq.gz
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))

# Display file lists
fnFs
fnRs

# Extract sample names from FASTQ filenames
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names


############################################################
# STEP 3: Check sequencing read quality
#
# Why this step?
# Quality plots help decide where to trim low-quality read ends.
# In this project, truncLen = c(240, 200) worked well.
############################################################

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


############################################################
# STEP 4: Filter and trim reads
#
# Why this step?
# This removes poor-quality reads, reads with ambiguous bases,
# and trims reads to selected lengths.
############################################################

# Create folder for filtered reads
filt_path <- file.path(path, "filtered")
dir.create(filt_path, showWarnings = FALSE)

# Define output filenames for filtered reads
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter and trim
out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(240, 200),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE
)

out

# Save filtering summary
filter_summary <- as.data.frame(out)
filter_summary$Percent_retained <- round(
  filter_summary$reads.out / filter_summary$reads.in * 100, 2
)

filter_summary
write.csv(filter_summary, "filtering_summary.csv")


############################################################
# STEP 5: Learn sequencing error rates
#
# Why this step?
# DADA2 learns the error pattern from the sequencing data.
# This helps distinguish true biological sequences from noise.
############################################################

errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

# Plot learned error models
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)


############################################################
# STEP 6: Denoise reads
#
# Why this step?
# Denoising converts noisy sequencing reads into exact sequence
# variants called ASVs.
############################################################

dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)

# Check first sample result
dadaFs[[1]]
dadaRs[[1]]


############################################################
# STEP 7: Check paired-end merging
#
# Why this step?
# Normally forward and reverse reads are merged. In this dataset,
# some samples failed to merge, so the final ASV table was made
# using forward reads only.
############################################################

mergers <- mergePairs(
  dadaFs, filtFs,
  dadaRs, filtRs,
  verbose = TRUE
)

head(mergers[[1]])

# NOTE:
# If some samples show 0 merged reads, continue with the
# forward-only workflow below.


############################################################
# STEP 8: Create ASV table using forward reads only
#
# Why this step?
# Forward-only analysis keeps all 10 samples and avoids losing
# samples that failed paired-end merging.
############################################################

seqtab <- makeSequenceTable(dadaFs)

dim(seqtab)
table(nchar(getSequences(seqtab)))


############################################################
# STEP 9: Remove chimeric sequences
#
# Why this step?
# Chimeras are artificial sequences formed during PCR.
# Removing them improves data quality.
############################################################

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = FALSE,
  verbose = TRUE
)

dim(seqtab.nochim)

# Proportion of reads retained after chimera removal
sum(seqtab.nochim) / sum(seqtab)

# Save ASV table
write.csv(seqtab.nochim, "ASV_table_nochim_forward_only.csv")
saveRDS(seqtab.nochim, "seqtab_nochim_forward_only.rds")


############################################################
# STEP 10: Track reads through the DADA2 pipeline
#
# Why this step?
# This table shows how many reads were retained at each stage.
# It is useful for reports and presentations.
############################################################

getN <- function(x) sum(getUniques(x))

track <- cbind(
  input = out[, 1],
  filtered = out[, 2],
  denoisedF = sapply(dadaFs, getN),
  nonchim = rowSums(seqtab.nochim)
)

rownames(track) <- sample.names

track
write.csv(track, "read_tracking_forward_only.csv")


############################################################
# STEP 11: Assign taxonomy using SILVA database
#
# Why this step?
# The ASVs are DNA sequences. Taxonomy assignment tells us
# which bacteria they likely belong to.
############################################################

# Required file:
# silva_nr99_v138.1_train_set.fa.gz
# It must be present in D:/Dixcy

taxa <- assignTaxonomy(
  seqtab.nochim,
  "silva_nr99_v138.1_train_set.fa.gz",
  multithread = FALSE
)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

# Save taxonomy table
write.csv(taxa, "taxonomy_table_forward_only.csv")
saveRDS(taxa, "taxonomy_forward_only.rds")


############################################################
# STEP 12: Summarise taxonomy assignment
#
# Why this step?
# This shows how many ASVs were classified at each taxonomic level.
############################################################

tax_summary <- data.frame(
  Level = colnames(taxa),
  Assigned_ASVs = colSums(!is.na(taxa))
)

tax_summary
write.csv(tax_summary, "taxonomy_assignment_summary.csv", row.names = FALSE)


############################################################
# STEP 13: Prepare metadata from SRA Run Table
#
# Why this step?
# FASTQ files contain sequence data only. Metadata tells us
# which sample is seawater control, intermediate stress, or
# freshwater stress.
############################################################

sra <- read.csv("SraRunTable.csv")

names(sra)

# Create clean metadata table using only samples present in ASV table
metadata_updated <- sra %>%
  filter(Run %in% rownames(seqtab.nochim)) %>%
  select(
    SampleID = Run,
    Run,
    LibraryName = Library.Name,
    SampleName = Sample.Name,
    Replicate = replicate,
    Tissue = isolation_source,
    Host = Host
  )

# Assign condition using SampleName:
# SW  = seawater control
# SFW = intermediate hypotonic stress
# FW  = freshwater hypotonic stress
metadata_updated$Condition <- ifelse(
  grepl("^SW", metadata_updated$SampleName),
  "SW_control",
  ifelse(
    grepl("^SFW", metadata_updated$SampleName),
    "SFW_intermediate",
    "FW_stress"
  )
)

metadata_updated$Stress <- ifelse(
  metadata_updated$Condition == "SW_control",
  "Seawater_control",
  ifelse(
    metadata_updated$Condition == "SFW_intermediate",
    "Intermediate_hypotonic_stress",
    "Freshwater_hypotonic_stress"
  )
)

rownames(metadata_updated) <- metadata_updated$SampleID
metadata_updated$SampleID <- rownames(metadata_updated)

metadata_updated
write.csv(metadata_updated, "metadata_updated.csv", row.names = TRUE)


############################################################
# STEP 14: Fix sample names and create phyloseq object
#
# Why this step?
# phyloseq combines ASV table, taxonomy table, and metadata
# into one object for easy analysis.
############################################################

# Clean row names in ASV table if they still contain filtered file suffix
rownames(seqtab.nochim) <- gsub(
  "_F_filt.fastq.gz",
  "",
  rownames(seqtab.nochim),
  fixed = TRUE
)

metadata <- read.csv("metadata_updated.csv", row.names = 1)

# Make sure sample names match
all(rownames(seqtab.nochim) == rownames(metadata))

# Create phyloseq object
otu <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)
tax <- tax_table(taxa)
sam <- sample_data(metadata)

ps <- phyloseq(otu, tax, sam)

ps
saveRDS(ps, "final_phyloseq_object_with_metadata.rds")


############################################################
# STEP 15: Calculate alpha diversity
#
# Why this step?
# Alpha diversity measures diversity within each sample.
# Observed = richness
# Shannon = richness + evenness
# Simpson = dominance/evenness
############################################################

alpha_div <- estimate_richness(
  ps,
  measures = c("Observed", "Shannon", "Simpson")
)

alpha_div$SampleID <- rownames(alpha_div)

metadata2 <- data.frame(sample_data(ps))
metadata2$SampleID <- rownames(metadata2)

alpha_div_meta <- dplyr::left_join(
  alpha_div,
  metadata2,
  by = "SampleID"
)

alpha_div_meta
write.csv(alpha_div_meta, "alpha_diversity_by_condition.csv", row.names = FALSE)


############################################################
# STEP 16: Plot richness by condition
#
# What this shows:
# Number of ASVs in each osmotic condition.
############################################################

ggplot(alpha_div_meta, aes(x = Condition, y = Observed, fill = Condition)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 3) +
  theme_minimal() +
  labs(
    title = "Observed ASV Richness Across Osmotic Conditions",
    x = "Condition",
    y = "Observed ASVs"
  )

ggsave("richness_by_condition.png", width = 8, height = 5, dpi = 300)


############################################################
# STEP 17: Plot Shannon diversity by condition
#
# What this shows:
# Whether microbial diversity and evenness change under stress.
############################################################

ggplot(alpha_div_meta, aes(x = Condition, y = Shannon, fill = Condition)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 3) +
  theme_minimal() +
  labs(
    title = "Shannon Diversity Across Osmotic Conditions",
    x = "Condition",
    y = "Shannon Diversity Index"
  )

ggsave("shannon_by_condition.png", width = 8, height = 5, dpi = 300)


############################################################
# STEP 18: Plot Simpson diversity by condition
#
# What this shows:
# Whether one/few microbes dominate the community.
############################################################

ggplot(alpha_div_meta, aes(x = Condition, y = Simpson, fill = Condition)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 3) +
  theme_minimal() +
  labs(
    title = "Simpson Diversity Across Osmotic Conditions",
    x = "Condition",
    y = "Simpson Diversity Index"
  )

ggsave("simpson_by_condition.png", width = 8, height = 5, dpi = 300)


############################################################
# STEP 19: Create alpha diversity summary table
#
# Why this step?
# Gives mean and standard deviation for each condition.
############################################################

alpha_summary <- alpha_div_meta %>%
  group_by(Condition) %>%
  summarise(
    Mean_Richness = mean(Observed),
    SD_Richness = sd(Observed),
    Mean_Shannon = mean(Shannon),
    SD_Shannon = sd(Shannon),
    Mean_Simpson = mean(Simpson),
    SD_Simpson = sd(Simpson),
    n = n()
  )

alpha_summary
write.csv(alpha_summary, "alpha_diversity_summary_by_condition.csv", row.names = FALSE)


############################################################
# STEP 20: Plot microbial phyla across conditions
#
# What this shows:
# Broad bacterial groups present in each osmotic condition.
############################################################

ps_phylum <- tax_glom(ps, taxrank = "Phylum")
ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

phylum_df <- psmelt(ps_phylum_rel)

ggplot(phylum_df, aes(x = Condition, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(
    title = "Microbial Phyla Across Osmotic Conditions",
    x = "Condition",
    y = "Relative Abundance"
  )

ggsave("phyla_by_condition.png", width = 9, height = 5, dpi = 300)


############################################################
# STEP 21: Plot top 10 microbial genera by condition
#
# What this shows:
# More specific bacterial genera that dominate each group.
############################################################

ps_genus <- tax_glom(ps, taxrank = "Genus")
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

genus_df <- psmelt(ps_genus_rel)

top10_genera <- genus_df %>%
  filter(!is.na(Genus)) %>%
  group_by(Genus) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(MeanAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

genus_df_top10 <- genus_df %>%
  mutate(Genus = ifelse(Genus %in% top10_genera, Genus, "Others"))

ggplot(genus_df_top10, aes(x = Condition, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(
    title = "Top 10 Microbial Genera Across Osmotic Conditions",
    x = "Condition",
    y = "Relative Abundance"
  )

ggsave("top10_genera_by_condition.png", width = 9, height = 5, dpi = 300)


############################################################
# STEP 22: Save top genera table
#
# Why this step?
# Gives the dominant genera in each condition.
############################################################

genus_condition_summary <- genus_df %>%
  filter(!is.na(Genus)) %>%
  group_by(Condition, Genus) %>%
  summarise(
    Mean_Relative_Abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Condition, desc(Mean_Relative_Abundance))

top_genera_by_condition <- genus_condition_summary %>%
  group_by(Condition) %>%
  slice_head(n = 10)

top_genera_by_condition
write.csv(top_genera_by_condition, "top_genera_by_condition.csv", row.names = FALSE)


############################################################
# STEP 23: Beta diversity using Bray-Curtis PCoA
#
# Why this step?
# Beta diversity compares microbial community differences
# between samples/groups.
############################################################

bray_dist <- phyloseq::distance(ps, method = "bray")

ord_pcoa <- ordinate(ps, method = "PCoA", distance = bray_dist)

plot_ordination(ps, ord_pcoa, color = "Condition") +
  geom_point(size = 4) +
  theme_minimal() +
  labs(
    title = "PCoA Plot of Gill Microbiome Across Osmotic Conditions",
    color = "Condition"
  )

ggsave("pcoa_by_condition.png", width = 7, height = 5, dpi = 300)


############################################################
# STEP 24: PERMANOVA statistical test
#
# Why this step?
# PERMANOVA tests whether microbial communities are
# significantly different between conditions.
############################################################

otu_mat <- as(otu_table(ps), "matrix")

if (taxa_are_rows(ps)) {
  otu_mat <- t(otu_mat)
}

meta_df <- data.frame(sample_data(ps))

adonis_condition <- adonis2(
  otu_mat ~ Condition,
  data = meta_df,
  method = "bray",
  permutations = 999
)

adonis_condition
write.csv(as.data.frame(adonis_condition), "PERMANOVA_condition_results.csv")


############################################################
# STEP 25: COI genetic diversity / phylogenetic analysis
#
# Why this step?
# COI is a common DNA barcode gene used for fish identification
# and phylogenetic relationship analysis.
############################################################

# Required file:
# COI_Oryzias_sequences.fasta

coi_seq <- readDNAStringSet("COI_Oryzias_sequences.fasta")

coi_seq
width(coi_seq)

# Align COI sequences using MUSCLE
coi_alignment <- msa(
  coi_seq,
  method = "Muscle"
)

coi_alignment

# Convert alignment into DNAbin object
coi_alignment_dnabin <- as.DNAbin(coi_alignment)

# Check aligned lengths
sapply(coi_alignment_dnabin, length)

# Calculate genetic distance using Kimura 2-parameter model
coi_dist <- dist.dna(
  coi_alignment_dnabin,
  model = "K80",
  pairwise.deletion = TRUE
)

# Build Neighbor-Joining tree
coi_tree <- nj(coi_dist)

# Plot tree
plot(
  coi_tree,
  main = "Phylogenetic Tree of Oryzias Species Based on COI Gene",
  cex = 0.8
)

add.scale.bar()

# Save tree as PDF
pdf("COI_Oryzias_phylogenetic_tree.pdf", width = 8, height = 6)

plot(
  coi_tree,
  main = "Phylogenetic Tree of Oryzias Species Based on COI Gene",
  cex = 0.8
)

add.scale.bar()

dev.off()

# Save tree as PNG
png("COI_Oryzias_phylogenetic_tree.png", width = 1200, height = 900, res = 150)

plot(
  coi_tree,
  main = "Phylogenetic Tree of Oryzias Species Based on COI Gene",
  cex = 0.8
)

add.scale.bar()

dev.off()


############################################################
# STEP 26: Check final output files
#
# Why this step?
# Confirms that the important final outputs were generated.
############################################################

list.files(pattern = "condition|COI|PERMANOVA|alpha|richness|shannon|simpson|phyla|genera|pcoa")


############################################################
# END OF SCRIPT
#
# Main final outputs:
#   - alpha_diversity_by_condition.csv
#   - alpha_diversity_summary_by_condition.csv
#   - richness_by_condition.png
#   - shannon_by_condition.png
#   - simpson_by_condition.png
#   - phyla_by_condition.png
#   - top10_genera_by_condition.png
#   - pcoa_by_condition.png
#   - PERMANOVA_condition_results.csv
#   - COI_Oryzias_phylogenetic_tree.png
#   - COI_Oryzias_phylogenetic_tree.pdf
############################################################
