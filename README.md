# Transcriptomics
Assignment 2
# Introduction
The yeast Saccharomyces cerevisiae is primarily recognized for its role in fermentative metabolism; however, certain strains have evolved the ability to transition into a "flor" or "velum" growth phase (Legras et al., 2016). This process is central to the biological aging of specialized wines, such as Sherry, where the yeast must transition from a fermentative state to a highly organized, buoyant biofilm at the air-liquid interface which is unlike typical fermentative yeast that settles at the bottom of the vessel (Mardanov et al., 2020). This biofilm serves as a protective barrier, preventing the excessive oxidation of the wine while the yeast switches its metabolism to the oxidation of ethanol and glycerol. The transition to a velum state is triggered by a combination of nutrient depletion, specifically the exhaustion of fermentable sugars—and the presence of high ethanol concentrations (Mardanov et al., 2020).
Biofilm development in S. cerevisiae proceeds through distinct stages that include initial surface adhesion and structural maturation. These transitions are accompanied by coordinated shifts in gene expression that reprogram cellular metabolism, stress tolerance, and cell wall architecture (Bojsen et al., 2012; Bouyx et al., 2021). As biofilms mature, transcriptional programs increasingly favor stress adaptation, nutrient scavenging, and structural reinforcement, reflecting a transition from growth-oriented physiology to community-based resilience (Sauer et al., 2022).
Early biofilm states are typically enriched for metabolic adjustment and initial stress response pathways, whereas mature biofilms exhibit stronger enrichment of cell wall organization, extracellular matrix assembly, and environmental resistance programs (Sauer et al., 2022).
For the alignment and quantification, pseudo-aligners like Salmon and Kallisto offer significant speed advantages, the STAR + RSEM pipeline was prioritized for this study because it maps every read to specific coordinates on the chromosomes, accounting for introns and splice junctions unlike pseudo-aligners. Differential expression analysis alone does not fully capture the biological implications of transcriptional reprogramming. Functional enrichment strategies, including over-representation analysis (ORA), allow differentially expressed genes to be categorized into defined biological processes and pathways (Carbon et al., 2021; "The Gene Ontology resource: enriching a GOld mine," 2021). In S. cerevisiae, curated annotations from the Saccharomyces Genome Database provide high-confidence gene-to-function mappings that enable biologically meaningful interpretation of transcriptomic datasets (Cherry et al., 2012).
High-throughput RNA sequencing (RNA-seq) provides a powerful approach for quantifying genome-wide gene expression changes during biofilm development. So, in this assignment, RNAseq transcriptomic profiling was used to identify and characterize gene expression changes across three major stages of Saccharomyces cerevisiae biofilm development (Early, thin and mature). Also, ORA was carried out to categorize top differentially expressed genes according to biological process. This with the statistical analysis assisted in identifying the transcriptional programs involved in biofilm progression. 
# Methodology
## Data acquisition
The dataset of the Saccharomyces cerevisiae flor yeast strain (PRJNA592304) was gotten from the NCBI SRA and the nine samples divided across three stages were downloaded. The raw sequence data was then retrieved into a FASTQ format using the fasterq-dump utility from the SRA Toolkit. The reference genome and annotation were also downloaded from the data base in fasta and gtf format.
``` bash
fasterq-dump --split-files – “each sample name”
```
## Quality Control
Initial sequence quality was assessed using FastQC (v0.12.1) to evaluate per-base quality scores, GC content, adapter contamination and other parameters. This is further viewed using multiplot to view the nine samples. 
```bash
fastqc *.fastq -t 4 -o
```
## Genomic Alignment and Quantification
Reads were mapped to the Saccharomyces cerevisiae reference genome (R64-1-1) using STAR (Spliced Transcripts Alignment to a Reference, v2.7.11b). STAR was selected for its splice-aware alignment capabilities, which are essential for accurately mapping eukaryotic transcripts that may contain introns. Following alignment, transcript abundance was quantified using RSEM (RNA-Seq by Expectation-Maximization, v1.3.3). RSEM provides gene-level counts, producing a raw count matrix for downstream analysis. So the star index was built with this;
``` bash
STAR \
  --runThreadN 4 \
  --genomeDir ../star_index \
  --readFilesIn SRR10551657.fastq \
  --outFileNamePrefix ../star_alignments/SRR10551657_ \
  --outSAMtype BAM SortedByCoordinate
```
After which the RSEM reference was prepared; 
```bash
rsem-prepare-reference \
  --gtf Saccharomyces_cerevisiae.R64-1-1.113.gtf \
  Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
  rsem_reference/yeast_R64
```
which was followed my rsem calculate expression;
```bash
rsem-calculate-expression \
  --star \
  --estimate-rspd \
  -p 4 \
  ~/yeast_biofilm/raw_data/SRR10551657.fastq \
  rsem_reference \
  SRR10551657
```
## Differential Expression Analysis
Statistical analysis was performed in the R studio. Transcript-level estimates from RSEM were imported into the DESeq2 using tximport. DESeq2 was chosen because it employs a negative binomial distribution model and empirical Bayes shrinkage to estimate dispersions and fold changes (Love et al., 2015).
```bash
dir <- "~/ubunututo windows/rsem_output"
samples <- read.csv("~/ubunututo windows/metadata.csv")
files <- file.path(dir, paste0(samples$sample, ".genes.results"))
names(files) <- samples$sample
Import data using trixmport
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
```
## Functional Annotation and Over-Representation Analysis (ORA)
Biological implications of the DEGs were characterized using the clusterProfiler package on R. Afterwards, Over-Representation Analysis (ORA) based on the Gene Ontology (GO) database was performed focusing on Biological Processes (BP). The org.Sc.sgd.db package provided the yeast-specific genome annotations (Wu et al., 2021).
# Results and Discussion 
## Data Quality and Transcriptomic Structure
After the data has been downloaded in a fastq format, fastQC was employed on each of the nine samples to ensure the quality of the raw reads, confirming the mean sequence quality (Phred scores) and the absence of significant adapter contamination before proceeding to alignment. It is of great importance to carry out this phase to ensure that other downstream analysis encounter little or no errors.
Following alignment with STAR and quantification via RSEM, Principal Component Analysis (PCA), a dimensionality reduction was performed to assess the overall data structure.
![image_alt](https://github.com/Sodjoh/Transcriptomics/blob/main/Screenshot%20from%202026-03-01%2001-10-42.png)
![image_alt](https://github.com/Sodjoh/Transcriptomics/blob/main/Screenshot%20from%202026-03-01%2012-29-52.png)
In figure 3, it could be seen that the three replicates for each group (the dots of the same color) are tightly packed together which shows high reproducibility. The PCA plot clearly demonstrates a progressive transcriptomic shift along PC1, accounting for 71% of the total variance. This separation strongly correlates with the chemical metadata; as the environment moves from an Early stage with 12.4% ethanol to a Mature stage of 9.6% characterized by extreme aldehyde toxicity (668.8 mg/l) and total glucose depletion.
![image_alt](https://github.com/Sodjoh/Transcriptomics/blob/main/Rplot01.png)
