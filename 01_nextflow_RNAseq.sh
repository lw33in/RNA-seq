# Nextflow command to pre-process RNA-seq raw data

# RNA-seq raw data: PE150, fastq.gz files
# Pipeline: nf-core/rnaseq (https://nf-co.re/rnaseq), version: 3.12.0
# Input csv file format: sample[sample_name], fastq_1[fastq.gz_file_name], fastq_2[fastq.gz_file_name],strandedness (default=auto)
# Reference genome: GRCh38 Gencode v32, download command (provided by Cellranger): wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
# Mannual download of reference genome see README.md
# FASTQ preprocessor: fastp, remove ribosomal RNA

nextflow run nf-core/rnaseq --input example_rnaseq_config.csv -profile docker --aligner star_salmon \
  -r 3.12.0 --outdir output_dir --igenomes_ignore --fasta refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  --gtf refdata-gex-GRCh38-2020-A/genes/genes.gtf --star_index refdata-gex-GRCh38-2020-A/star/ \
  --trimmer fastp --gencode --remove_ribo_rna
