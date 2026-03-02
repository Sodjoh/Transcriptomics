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
![image_alt](https://github.com/Sodjoh/Transcriptomics/blob/main/Rplot.png)
## Differential Gene Expression Analysis
Two pairwise comparisons were conducted to capture the evolution of the biofilm: Thin vs. Early and Mature vs. Early. Using a threshold of padj< 0.05 and absolute log2FC > 1, this resulted into several thousand differentially expressed genes (DEGs) of about 6807 genes.
|6807 genes|Mature Vs Early|Thin Vs Early|
|:---|---|---:|
|Padj|<0.05|<0.05|
|LFC > 0 (up)|1638, 24%|1097, 16%|
|LFC < 0 (down)|1507, 22%|1072, 16%|
|Outliers|3, 0.044%|3, 0.044%|
|Low counts|129, 1.9%|395, 5.8%|
### MA Plots
MA plots demonstrated symmetric distributions of log2 fold changes around zero, with a substantial proportion of genes exhibiting large effect sizes (|log2FC| > 1). In figure 4, it shows a significant increase in red points (significant genes). This shows that as the biofilm matures, more genes are recruited into the survival response even as the condition became unfavorable. Similar to the findings by Mardanov et al. (2020), the surge in significant genes in the Mature stage directly correlates with the "hostile" high-aldehyde environment described in the study.
![image_alt](https://github.com/Sodjoh/Transcriptomics/blob/main/Rplot01.png)
![image_alt](https://github.com/Sodjoh/Transcriptomics/blob/main/Rplot02.png)
### Volcano Plot
The volcano plot visualizes the statistical significance against the magnitude of change for the Mature vs. Early comparison. Red points represent significantly Upregulated genes, while blue points represent Downregulated genes padj < 0.05, log2FC > 1. This result illustrates a huge transcriptomic shift. A significant number of genes were upregulated, reflecting a proactive response to the increasingly hostile environment. Specifically, the magnitude of change was most pronounced in genes related to stress response, coinciding with the peak Aldehyde concentration of 668.8 mg/l and total Glucose depletion (< 0.1 g/l) observed in the chemical metadata.
![image_alt](https://github.com/Sodjoh/Transcriptomics/blob/main/Rplot03.png)
### Heatmap
These top 10 genes represent the strongest biological responses to the 109-day maturation period. It could be seen in figure 6 that there is a difference between the Early (green) and Mature (red) columns. The top genes like IRT1 and ETS2-1 are highly upregulated (red) in the mature phase but downregulated (blue) in the early phase.
![image_alt](https://github.com/Sodjoh/Transcriptomics/blob/main/Rplot04.png)
### Venn Diagram
The Venn diagram illustrates the distribution of DEGs (padj < 0.05, |log2FC| > 1) for the Thin vs. Early and Mature vs. Early comparisons. While 728 genes represent a core biofilm response, the 1,347 genes unique to the Mature stage highlight the extensive metabolic recruitment required to survive peak aldehyde toxicity and nutrient starvation at maturation.
![image_alt](https://github.com/Sodjoh/Transcriptomics/blob/main/Rplot05.png)
### Functional Annotation
The Over-Representation Analysis (ORA) provides the biological reason behind the shift in the stages of development. In the Thin stage, the high significance of "carbohydrate catabolic processes" confirms that the yeast tried to source for energy as the glucose level dropped. By the time it got to the matured stage, it could be seen in the parameter (Energy derivation by oxidation of organic compounds) that the yeast has already gotten a new source of energy. As noted by Mardanov et al. (2020), flor yeast must transition to the aerobic oxidation of secondary carbon sources to survive. The additional metadata provide the molecular evidence for this, as the yeast maintains energy production through the consumption of Ethanol (decreasing from 12.4% to 9.6%) and Glycerol (8.5 to 6.8 g/l) once primary sugars are exhausted. There is another parameter that is worthy of note, Transmembrane transport; this remained constant through the development, which could indicate that the yeast is constantly modifying its cell wall.
![image_alt](https://github.com/Sodjoh/Transcriptomics/blob/main/Rplot06.png)
# Conlusion
This study successfully implemented a STAR-RSEM-DESeq2 pipeline to investigate the transcriptomic shift of Saccharomyces cerevisiae during biofilm (velum) maturation. The analysis revealed that the transition from Early to Mature stages is not merely a cessation of growth, but an innate metabolic adaptation. The differential expression of key genes, visualized through PCA, heatmaps and Volcano plots, highlights the yeast's response to the harsh environment of wine aging. Ultimately, these results demonstrate how S. cerevisiae manages the delicate balance of protecting wine from oxidation while surviving extreme chemical conditions.
# References
Bojsen, R. K., Andersen, K. S., & Regenberg, B. (2012). Saccharomyces cerevisiae—a model to uncover molecular mechanisms for yeast biofilm biology. FEMS Immunology & Medical Microbiology, 65(2), 169-182. 

Bouyx, C., Schiavone, M., & François, J. M. (2021). FLO11, a developmental gene conferring impressive adaptive plasticity to the yeast Saccharomyces cerevisiae. Pathogens, 10(11), 1509. 

Carbon, S., Douglass, E., Good, B. M., Unni, D. R., Harris, N. L., & Mungall, C. J. (2021). et alThe gene ontology resource: enriching a GOld mine. Nucleic Acids Res, 49(D1), D325-D323. 

Cherry, J. M., Hong, E. L., Amundsen, C., Balakrishnan, R., Binkley, G., Chan, E. T., Christie, K. R., Costanzo, M. C., Dwight, S. S., & Engel, S. R. (2012). Saccharomyces Genome Database: the genomics resource of budding yeast. Nucleic acids research, 40(D1), D700-D705. 

The Gene Ontology resource: enriching a GOld mine. (2021). Nucleic acids research, 49(D1), D325-D334. 

Legras, J.-L., Moreno-Garcia, J., Zara, S., Zara, G., Garcia-Martinez, T., Mauricio, J. C., Mannazzu, I., Coi, A. L., Bou Zeidan, M., & Dequin, S. (2016). Flor yeast: new perspectives beyond wine aging. Frontiers in Microbiology, 7, 503. 

Love, M. I., Anders, S., Kim, V., & Huber, W. (2015). RNA-Seq workflow: gene-level exploratory analysis and differential expression. F1000Res, 4, 1070. https://doi.org/10.12688/f1000research.7035.1 

Mardanov, A. V., Eldarov, M. A., Beletsky, A. V., Tanashchuk, T. N., Kishkovskaya, S. A., & Ravin, N. V. (2020). Transcriptome profile of yeast strain used for biological wine aging revealed dynamic changes of gene expression in course of flor development. Frontiers in Microbiology, 11, 538. 

Sauer, K., Stoodley, P., Goeres, D. M., Hall-Stoodley, L., Burmølle, M., Stewart, P. S., & Bjarnsholt, T. (2022). The biofilm life cycle: expanding the conceptual model of biofilm formation. Nature Reviews Microbiology, 20(10), 608-620. 

Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., Feng, T., Zhou, L., Tang, W., Zhan, L., Fu, X., Liu, S., Bo, X., & Yu, G. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. Innovation (Camb), 2(3), 100141. https://doi.org/10.1016/j.xinn.2021.100141
