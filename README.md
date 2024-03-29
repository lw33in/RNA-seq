# RNA-seq
## RNA-seq Experimental Design Best Practice
- **Sample Preparation**:
  * Highly recommend snap-frozen cell pellets
- **Library Construction**:
  * mRNA library construction for coding mRNA studies
  * Total RNA library construction for long noncoding RNA studies
- **Data Type**:
  * Paired-end 150bp
  * Short-read sequencing is most commonly used. Long-read sequencing is recommended for isoform-emphasized analysis
- **Replicates**:
    * Highly recommend >= 3 biological replicates for statistical analysis
    * Technical replicates are typically unnecessary for biology discovery projects
- **Sequencing Depth**:
    * At least 20M read pairs recommended for most cell line samples
    * At least 30M read pairs recommended for patient samples and long noncoding RNA studies
- **Batch Effects**:
    * Ideally, samples in the same studies are sequenced together. If samples need to sequenced in batches, ensure that samples are randomized and replicates are eveny distributed across batches
- **Platform**:
    * Illumina NovaSeq 6000

## RNA-seq Experimental and Analysis Workflow

> Note: This is a standard workflow designed for bulk mRNA sequencing. Based on the specific biological questions, modifications of this workflow might be needed.

- **Pre-processing**:
    * In this workflow, [nf-core/rnaseq/3.12.0](https://nf-co.re/rnaseq/3.12.0) is used to process FASTQ files, performing quality control, trimming, alignment, and produce gene expression matrix. v3.12.0 is an older version of the pipeline. If you wish to use the latest release, please visit https://nf-co.re/rnaseq/
- **Downstream analysis**:
  * [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) are two R tools used for differential gene expression analysis.


```mermaid
flowchart LR

subgraph A["Bench Experiments"]
direction TB
    A1(Experimental Biology) --> A2(Sample Preparation) --> A3(Library Construction)
    --> A4(Sequencing)
end

subgraph B["Pipeline Processing"]
direction TB
    B1(Raw Reads)
    B2(Sequence Quality Control)
    B2.1("Read Quality Metrics
        GC-content Distribution
        Ribosomal RNA Removal")
    B3(Read Alignment)
    B4(Count Quantification)

    B1 --> |FastQC & MultiQC|B2--> |FastP & SortMeRNA|B2.1 --> |STAR|B3 --> |Salmon & MarkDuplicates|B4
end

subgraph C["Downstream Analysis"]
direction TB
    C1(Count Matrix) 
    C2(Quality Control) 
    C2.1("Gene-level Filtering
        Transcript-level Filtering") 
    C3(Normalization) 
    C4("Differential Expression
        Analysis")
    C5(Functional Analysis) 
    C5.1("GO Enrichment
        KEGG Enrichment
        Drug Target Anlysis
        Protein-protein Interaction
        WGCNA
        GSEA")
    C6(Data Visualization)
    C7(Analysis Report)

    C1 --> C2 --> C2.1 --> C3--> |DESeq2 & edgeR|C4 --> C5 --> C5.1 --> C6 --> C7
end

D["Multi-omics
    Integrated Analysis"]

A ==> B ==> C ==> D
```
- **Reference Genome**:
    * Gencode v32 reference genome
    * Command line download: <code>wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz</code>
    * Direct download: [Human reference (GRCh38) dataset required for Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest/)
- **Differential Expression Analysis package**:
    * Significance Cutoff: abs|FC| >= 1.5, fdr < 0.05

## RNA-seq Computational Analysis
- `01_nextflow_RNAseq.sh` General command lines used to process raw reads using nf-core pipelines.
   * example_rnaseq_config.csv:
  ```
  sample,fastq_1,fastq_2,strandedness
  ctrl_rep1,ctrl_rep1_S1_L002_R1_001.fastq.gz,ctrl_rep1_S1_L002_R2_001.fastq.gz,auto
  ctrl_rep2,ctrl_rep2_S1_L003_R1_001.fastq.gz,ctrl_rep2_S1_L003_R2_001.fastq.gz,auto
  ctrl_rep3,ctrl_rep3_S1_L004_R1_001.fastq.gz,ctrl_rep3_S1_L004_R2_001.fastq.gz,auto
  ```
  ```
  nextflow run nf-core/rnaseq --input example_rnaseq_config.csv -profile docker --aligner star_salmon \
  -r 3.12.0 --outdir output_dir --igenomes_ignore --fasta refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  --gtf refdata-gex-GRCh38-2020-A/genes/genes.gtf --star_index refdata-gex-GRCh38-2020-A/star/ \
  --trimmer fastp --gencode --remove_ribo_rna
  ```
  > Note If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with -profile test before running the workflow on actual data.
  
- `02_gene_filtering.R` Additional gene-level filtering criteria (stringent).
- `03_DESeq2.R` Downstream differential expression analysis (DEA) using DESeq2.
- `04_edgeR.R` Downstream DEA using edgeR.
- `05_enrichment.R` Enrichment analysis on significant differentially expressed genes (DEGs).
- `06_visualization.R` Data visualization of RNA-seq analysis results.
- `07_cell_type_specific.R` Identification and exploration of cell-type specific signature genes.

## Tools 
- R 3.12.0
- nextflow 3.12.0
- FastQC 0.11.9
- MultiQC 1.14
- FastP 0.23.4
- SortMeRNA 4.3.4
- STAR 2.6.1
- Salmon 1.10.1
- MarkDuplicates 3.0.0
- DESeq2 v1.42.0














