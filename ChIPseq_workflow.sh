#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define directories
BASE_DIR="/home/abhi39/Projects/CHIPseq"
DATA_DIR="$BASE_DIR/data"
READS_DIR="$DATA_DIR/reads"
CHR_BED_DIR="$DATA_DIR/chr_bed"

# Reference genome and annotation
REFERENCE_DIR="$BASE_DIR/genome"
GENOME_DIR="$REFERENCE_DIR/genome"
ANNOTATION_DIR="$REFERENCE_DIR/annotation"
INDEX_DIR="$GENOME_DIR/bwa_index"

# Output directories
OUTPUT_DIR="$BASE_DIR/output"
ALIGNMENT_DIR="$OUTPUT_DIR/3_alignment"
MACS_DIR="$OUTPUT_DIR/4_peak_calling"
QC_DIR="$OUTPUT_DIR/1_qc"
TRIMMED_DIR="$OUTPUT_DIR/2_trimmed"

# Create necessary directories
mkdir -p "$DATA_DIR" "$READS_DIR" "$CHR_BED_DIR" "$REFERENCE_DIR" "$GENOME_DIR" "$ANNOTATION_DIR" "$INDEX_DIR"
mkdir -p "$OUTPUT_DIR" "$ALIGNMENT_DIR" "$MACS_DIR" "$QC_DIR" "$TRIMMED_DIR"

# Quality Control using FastQC
echo "Running FastQC..."
fastqc -o "$QC_DIR" "$READS_DIR"/*.fastq.gz

# Adapter Trimming using Cutadapt
echo "Trimming reads with Cutadapt..."
for file in "$READS_DIR"/*.fastq.gz; do
    cutadapt -q 20 -m 30 -a AGATCGGAAGAGC -o "$TRIMMED_DIR/$(basename $file .fastq.gz)_trimmed.fastq.gz" "$file"
done

# Read Alignment using Bowtie2
echo "Aligning reads to reference genome with Bowtie2..."
for file in "$TRIMMED_DIR"/*_trimmed.fastq.gz; do
    sample_name=$(basename "$file" _trimmed.fastq.gz)
    bowtie2 -x "$INDEX_DIR/genome" -U "$file" -S "$ALIGNMENT_DIR/$sample_name.sam"
done

# Convert SAM to BAM, sort, and index
echo "Processing BAM files..."
for file in "$ALIGNMENT_DIR"/*.sam; do
    sample_name=$(basename "$file" .sam)
    samtools view -Sb "$file" | samtools sort -o "$ALIGNMENT_DIR/$sample_name.sorted.bam"
    samtools index "$ALIGNMENT_DIR/$sample_name.sorted.bam"
    rm "$file"  # Remove SAM files to save space
done

# Remove duplicates using Picard
echo "Removing duplicates using Picard..."
for file in "$ALIGNMENT_DIR"/*.sorted.bam; do
    sample_name=$(basename "$file" .sorted.bam)
    picard MarkDuplicates I="$file" O="$ALIGNMENT_DIR/$sample_name.dedup.bam" M="$ALIGNMENT_DIR/$sample_name.metrics.txt" REMOVE_DUPLICATES=true
    samtools index "$ALIGNMENT_DIR/$sample_name.dedup.bam"
done

# Peak Calling with MACS2
echo "Performing peak calling with MACS2..."
macs2 callpeak -t "$ALIGNMENT_DIR/SRR502327.dedup.bam" -c "$ALIGNMENT_DIR/SRR502225.dedup.bam" -f BAM -g hs -n STAT1_IFN6h -B -q 0.01 --outdir "$MACS_DIR"
macs2 callpeak -t "$ALIGNMENT_DIR/SRR502329.dedup.bam" -c "$ALIGNMENT_DIR/SRR502228.dedup.bam" -f BAM -g hs -n STAT1_IFN30h -B -q 0.01 --outdir "$MACS_DIR"

# Convert peak signals to BigWig format for visualization
echo "Converting peaks to BigWig format..."
bedtools genomecov -bg -ibam "$ALIGNMENT_DIR/SRR502327.dedup.bam" > "$MACS_DIR/STAT1_IFN6h.bedgraph"
bedtools genomecov -bg -ibam "$ALIGNMENT_DIR/SRR502329.dedup.bam" > "$MACS_DIR/STAT1_IFN30h.bedgraph"

bedGraphToBigWig "$MACS_DIR/STAT1_IFN6h.bedgraph" "$GENOME_DIR/hg38.chrom.sizes" "$MACS_DIR/STAT1_IFN6h.bw"
bedGraphToBigWig "$MACS_DIR/STAT1_IFN30h.bedgraph" "$GENOME_DIR/hg38.chrom.sizes" "$MACS_DIR/STAT1_IFN30h.bw"

# Find differentially bound peaks
echo "Identifying differential peaks..."
bedtools intersect -v -a "$MACS_DIR/STAT1_IFN30h_peaks.narrowPeak" -b "$MACS_DIR/STAT1_IFN6h_peaks.narrowPeak" > "$MACS_DIR/STAT1_IFN30h_unique.bed"
bedtools intersect -v -a "$MACS_DIR/STAT1_IFN6h_peaks.narrowPeak" -b "$MACS_DIR/STAT1_IFN30h_peaks.narrowPeak" > "$MACS_DIR/STAT1_IFN6h_unique.bed"

# Motif analysis using HOMER
echo "Running HOMER for motif analysis..."
findMotifsGenome.pl "$MACS_DIR/STAT1_IFN6h_unique.bed" hg38 "$MACS_DIR/HOMER_IFN6h" -size 200 -bg "$MACS_DIR/STAT1_IFN30h_unique.bed"
findMotifsGenome.pl "$MACS_DIR/STAT1_IFN30h_unique.bed" hg38 "$MACS_DIR/HOMER_IFN30h" -size 200 -bg "$MACS_DIR/STAT1_IFN6h_unique.bed"

echo "ChIP-seq workflow completed successfully."
