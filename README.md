# Microbial and Genetic Diversity of Medaka Fish Under Osmotic Stress

## Project Overview

This project analyses how osmotic stress affects the gill-associated microbiome of medaka fish and adds a simple genetic diversity component using COI barcode sequences.

The microbial diversity analysis uses publicly available 16S rRNA amplicon sequencing data from **NCBI SRA BioProject PRJNA702883**. The dataset contains gill microbiome samples from *Oryzias melastigma* under different osmotic conditions. Because exact public microbiome datasets for *Oryzias dancena* are limited, *Oryzias melastigma*, a closely related medaka species, was used as a comparative model for osmotic stress.

The genetic diversity section uses COI barcode sequences from selected *Oryzias* species and *Danio rerio* as an outgroup to construct a phylogenetic tree.

---

## Project Title

**Microbial and Genetic Diversity of Medaka Fish Under Osmotic Stress**

---

## Objectives

1. Download and process 16S rRNA sequencing reads from NCBI SRA.
2. Generate an ASV table using the DADA2 pipeline in R.
3. Assign taxonomy using the SILVA reference database.
4. Calculate alpha diversity: Observed ASV richness, Shannon diversity, and Simpson diversity.
5. Study microbial distribution at phylum and genus levels.
6. Perform beta diversity analysis using Bray-Curtis dissimilarity and PCoA.
7. Test whether osmotic condition significantly affects microbial composition using PERMANOVA.
8. Add a genetic diversity component using COI barcode-based phylogenetic analysis.

---

## Dataset Information

### Microbiome Dataset

- **Database:** NCBI SRA
- **BioProject:** PRJNA702883
- **Host:** *Oryzias melastigma*
- **Sample type:** Gill microbiome
- **Assay type:** 16S rRNA amplicon sequencing
- **Platform:** Illumina MiSeq
- **Analysis type:** Microbial diversity under osmotic stress

### Samples Used

| Condition | Description | Samples |
|---|---|---|
| SW_control | Seawater control | SRR13758749, SRR13758750, SRR13758751, SRR13758752 |
| SFW_intermediate | Intermediate hypotonic stress | SRR13758972, SRR13758973, SRR13758974 |
| FW_stress | Freshwater hypotonic stress | SRR13758978, SRR13758979, SRR13758981 |

### Genetic Dataset

- **Database:** NCBI Nucleotide
- **Gene used:** COI / COX1 barcode gene
- **Species included:** Selected *Oryzias* species and *Danio rerio* as an outgroup
- **Output:** Neighbor-joining phylogenetic tree

---

## Repository Structure

```text
.
├── data/
│   ├── SraAccList.txt
│   ├── SraRunTable.csv
│   ├── metadata.csv
│   ├── metadata_updated.csv
│   └── COI_Oryzias_sequences.fasta
│
├── results/
│   ├── ASV_table_nochim_forward_only.csv
│   ├── taxonomy_table_forward_only.csv
│   ├── taxonomy_assignment_summary.csv
│   ├── alpha_diversity_by_condition.csv
│   ├── alpha_diversity_summary_by_condition.csv
│   ├── top_15_genera_table.csv
│   ├── top_genera_by_condition.csv
│   └── PERMANOVA_condition_results.csv
│
├── plots/
│   ├── richness_by_condition.png
│   ├── shannon_by_condition.png
│   ├── simpson_by_condition.png
│   ├── pcoa_by_condition.png
│   ├── phyla_by_condition.png
│   ├── phylum_relative_abundance.png
│   ├── top10_genera_by_condition.png
│   ├── top10_genera_relative_abundance.png
│   └── COI_Oryzias_phylogenetic_tree.png
│
├── scripts/
│   └── medaka_microbiome_analysis.R
│
├── presentation/
│   └── Medaka_Oryzias_Microbial_Genetic_Diversity_10min_PPT.pptx
│
├── docs/
│   └── Oryzias_Medaka_Project_R_Code_Step_by_Step_Document.docx
│
└── README.md
```

The exact folder structure can be modified depending on how the files are arranged locally.

---

## Software and R Packages Used

### Required Software

- R / RStudio
- SRA Toolkit
- Optional: WSL/Linux terminal for easier FASTQ handling

### R Packages

```r
library(dada2)
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ape)
library(seqinr)
library(msa)
library(Biostrings)
```

---

## Workflow Summary

```text
NCBI SRA FASTQ files
        ↓
Quality profile check
        ↓
Filtering and trimming
        ↓
DADA2 denoising
        ↓
Forward-read ASV table
        ↓
Chimera removal
        ↓
Taxonomy assignment using SILVA
        ↓
Phyloseq object creation
        ↓
Alpha diversity analysis
        ↓
Microbial distribution analysis
        ↓
Bray-Curtis PCoA
        ↓
PERMANOVA
        ↓
COI barcode phylogenetic tree
```

---

## Main Analysis Steps

### 1. Download SRA Data

The run accessions were saved in `SraAccList.txt`.

Example:

```bash
prefetch --option-file SraAccList.txt
```

Convert SRA files to FASTQ:

```bash
for run in SRR13758749 SRR13758750 SRR13758751 SRR13758752 SRR13758972 SRR13758973 SRR13758974 SRR13758978 SRR13758979 SRR13758981
 do
   fasterq-dump $run --split-files
 done
```

Compress FASTQ files:

```bash
gzip *.fastq
```

---

### 2. DADA2 Processing

The reads were processed using DADA2. Because paired-end merging failed for some samples, a **forward-read-only workflow** was used to retain all 10 samples.

Important filtering settings:

```r
truncLen = c(240, 200)
maxN = 0
maxEE = c(2, 2)
truncQ = 2
```

Forward reads were used to generate the final ASV table.

Final ASV result:

```text
10 samples × 915 ASVs
```

Chimera removal retained approximately:

```text
90.7% of reads
```

---

### 3. Taxonomy Assignment

Taxonomy was assigned using the SILVA v138.1 reference database.

Taxonomy assignment summary:

| Taxonomic level | Assigned ASVs |
|---|---:|
| Kingdom | 881 |
| Phylum | 645 |
| Class | 618 |
| Order | 587 |
| Family | 536 |
| Genus | 373 |

---

### 4. Alpha Diversity

Alpha diversity was calculated using:

- Observed ASV richness
- Shannon diversity index
- Simpson diversity index

The freshwater stress group showed reduced richness, Shannon diversity, and Simpson diversity compared with seawater control and intermediate conditions.

Important output plots:

- `richness_by_condition.png`
- `shannon_by_condition.png`
- `simpson_by_condition.png`

---

### 5. Microbial Distribution

Microbial composition was analysed at phylum and genus levels.

Dominant phyla included:

1. Proteobacteria
2. Fusobacteriota
3. Firmicutes

Dominant genera included:

1. *Pseudomonas*
2. *Vibrio*
3. *Cetobacterium*
4. *Shinella*
5. *Catenococcus*

Important output plots:

- `phyla_by_condition.png`
- `top10_genera_by_condition.png`

---

### 6. Beta Diversity and PCoA

Beta diversity was calculated using Bray-Curtis dissimilarity.

PCoA was used to visualize microbial community differences among osmotic conditions.

Important output plot:

- `pcoa_by_condition.png`

The PCoA plot showed clear separation between seawater control, intermediate stress, and freshwater stress samples.

---

### 7. PERMANOVA

PERMANOVA was performed to test whether osmotic condition significantly affected microbial community composition.

Result:

```text
R² = 0.6567
F = 6.6951
p = 0.002
```

Interpretation:

Osmotic condition significantly influenced the gill microbial community. Around 65.7% of microbial community variation was explained by condition.

---

### 8. Genetic Diversity Using COI Barcode

COI barcode sequences were downloaded from NCBI for selected *Oryzias* species and *Danio rerio*.

Steps:

1. Read COI FASTA sequences.
2. Align sequences using MUSCLE.
3. Calculate Kimura 2-parameter genetic distance.
4. Build a neighbor-joining tree.
5. Save the tree as PNG and PDF.

Important output files:

- `COI_Oryzias_phylogenetic_tree.png`
- `COI_Oryzias_phylogenetic_tree.pdf`

Interpretation:

The COI tree shows genetic relationships among *Oryzias* species. *Danio rerio* was used as an outgroup. This section supports the genetic diversity component of the project.

---

## Key Results

1. Freshwater hypotonic stress reduced microbial richness and diversity.
2. Proteobacteria dominated the gill microbiome across conditions.
3. Freshwater stress samples showed strong microbial community shifts.
4. Bray-Curtis PCoA showed clear condition-wise separation.
5. PERMANOVA confirmed a statistically significant condition effect.
6. COI barcode analysis showed phylogenetic relationships among related medaka species.

---

## Final Conclusion

The gill microbiome of medaka fish changed clearly under osmotic stress. Freshwater hypotonic stress was associated with lower microbial richness, Shannon diversity, and Simpson diversity. Microbial composition also shifted across conditions, with Proteobacteria being the dominant phylum. PERMANOVA confirmed that osmotic condition significantly affected gill microbial community structure. The COI barcode tree added a genetic diversity component by showing evolutionary relationships among *Oryzias* species.

---

## Limitations

- Exact public microbiome data for *Oryzias dancena* were limited, so *Oryzias melastigma* was used as a related medaka model.
- A forward-read-only DADA2 workflow was used because paired-end merging failed for some samples.
- One COI sequence was partial; therefore, pairwise deletion was used during distance calculation.
- The project focuses on bioinformatics analysis and does not include wet-lab validation.

---

## Suggested Future Work

1. Include more samples and replicate groups.
2. Use full paired-end data after adapter trimming optimization.
3. Perform differential abundance analysis.
4. Add stress-response genes such as NKCC/SLC12A2, ATP1A1, aquaporins, and HSP70.
5. Compare microbiome changes with host gene expression data.

---

## Presentation Files

A 10-minute undergraduate-friendly presentation was prepared with the important plots and key results.

Suggested slide order:

1. Title
2. Background
3. Dataset and study design
4. Bioinformatics workflow
5. Alpha diversity
6. Microbial distribution
7. Beta diversity and PERMANOVA
8. Genetic diversity using COI tree
9. Integrated interpretation
10. Conclusion and thank you

---

## Author

Prepared as an undergraduate-friendly bioinformatics project on microbial and genetic diversity of medaka fish under osmotic stress.
