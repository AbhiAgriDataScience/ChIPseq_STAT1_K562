# **ChIP-seq Analysis: STAT1 Binding in K562 Cells**

## **Overview**

This repository contains a bioinformatics pipeline for analyzing ChIP-seq data to identify STAT1 binding sites in the K562 human leukemia cell line treated with Interferon-alpha (IFNα). The study is based on ENCODE ChIP-seq data from GEO dataset **GSE31477**. The workflow involves data preprocessing, peak calling, differential analysis, motif discovery, and functional annotation.

## **Dataset Information**

The dataset consists of four ChIP-seq samples retrieved from **NCBI SRA (GSE31477)**:

| Sample ID | Description | Condition | Treatment | Purpose |
| ----- | ----- | ----- | ----- | ----- |
| SRR502225 | K562 Input (control) \- IFNα 6h | Control | No | Input DNA control for normalization |
| SRR502228 | K562 Input (control) \- IFNα 30h | Control | No | Input DNA control for normalization |
| SRR502327 | K562 STAT1 ChIP \- IFNα 6h | Treatment | Yes | Captures STAT1 binding sites after 6h IFNα treatment |
| SRR502329 | K562 STAT1 ChIP \- IFNα 30h | Treatment | Yes | Captures STAT1 binding sites after 30h IFNα treatment |

## **Workflow Overview**

This project follows a standard ChIP-seq analysis workflow:

### **1\. Quality Control (QC)**

* Raw FASTQ files are assessed using **FastQC**.  
* Adapter sequences and low-quality bases are removed using **Cutadapt**.

### **2\. Read Alignment**

* Trimmed reads are aligned to the **GRCh38 human genome** using **Bowtie2**.  
* BAM files are sorted and indexed using **Samtools**.

### **3\. Remove Duplicates (if necessary)**

* Duplicate reads are checked using **Picard MarkDuplicates**.  
* Duplication metrics are generated for assessment.

### **4\. Peak Calling**

* **MACS2** is used to identify STAT1 binding sites.  
* Input control samples are used for normalization.

### **5\. Differential Peak Analysis**

* Goal: Compare STAT1 binding across IFN-α time points.  
* **DiffBind (R)** was considered but not used due to single replicate per condition.  
* **macs2 bdgdiff** was attempted but required control bedgraphs.  
* **bedtools** was used to identify peaks unique to IFNα30h and IFNα6h.

### **6\. Motif Analysis**

* **HOMER** is used to find STAT1-associated motifs unique to each time point.  
* Background files for each time point were created for accurate enrichment analysis.

### **7\. Visualization**

* Peak signals converted to **BigWig format** for **IGV/UCSC genome browser**.  
* NarrowPeak and summit files are used for motif analysis.  
* HOMER results are visualized to interpret motif enrichment.

### **8\. Functional Annotation of Differentially Bound STAT1 Peaks**

* Annotation of differentially bound STAT1 peaks.  
* Genes without Entrez IDs were filtered before GO enrichment.  
* **HOMER (findGO)** was used for Gene Ontology (GO) enrichment.  
* GO enrichment results were visualized.

### **9\. Transcription Factor Enrichment Analysis**

* Classification of peaks:  
  * Distance-based classification (enhancer vs. promoter-associated peaks).  
  * Annotation-based classification (enhancer vs. promoter-associated peaks).  
* **HOMER** is used for TF enrichment on enhancer peaks.  
* Genes near enhancer-associated peaks are extracted.  
* **ChIP-X Enrichment Analysis 3 (ChEA3)** is used for TF enrichment.  
* Visualization of TF enrichment results.  
* Comparison of enriched TFs to known STAT1 target genes.

### **10\. Validation and Refinement of TF Enrichment Analysis**

* Overlay enriched TFs with motif analysis (**HOMER results**).  
* Compare TF motifs from **HOMER** with **ChEA3** enriched TFs.  
* Assess overlap of transcription factors identified using different approaches.  
* Build a **TF-TF regulatory network** based on known interactions.

### **11\. Pathway Analysis (GO/KEGG Enrichment)**

* Genes near enhancer-associated peaks are extracted.  
* **KEGG and GO enrichment analysis** is performed to identify enriched pathways.

## **Repository Structure**

ChIPseq\_STAT1\_K562/  
├── data/                  \# Raw and processed data files  
├── genome/  
├── scripts/               \# Detailed ChIPseq workflow with explanation in jupyter notebook  
├── output/               \# Outputs from each analysis step  
│   ├── 1\_qc\_results/                \# Quality control reports  
│   ├── 2\_trimmed\_reads/   \# Trimmed reads  
│   ├── 3\_alignment/         \# BAM files and indexes  
│   ├── 4\_peak calling/             \# MACS2 peak calling results  
│   ├── 5\_visualization/     \# Plots and genome browser tracks  
│   ├── 6\_peak\_statistics/       
│   ├── 7\_motif\_analysis/    \# HOMER motif enrichment results  
│   ├── 8\_peak\_annotation/     \# Plots and genome browser tracks  
│   ├── 9\_GO\_KEGG/        \# GO/KEGG enrichment results  
├── README.md              \# Project documentation  
├── requirements.txt       \# Software dependencies  
└── ChIPseq\_workflow.sh    \# Main workflow script

## **Software and Dependencies**

To reproduce this analysis, install the following dependencies:

### **Python:**

pip install numpy pandas matplotlib seaborn macs2 bedtools

### **R:**

install.packages(c("DiffBind", "ggplot2", "ChIPseeker"))

### **Other Tools:**

* **FastQC**  
* **Cutadapt**  
* **Bowtie2**  
* **Samtools**  
* **MACS2**  
* **bedtools**  
* **HOMER**  
* **ChEA3**

## **How to Run the Analysis**

Clone the repository:  
git clone https://github.com/yourusername/ChIPseq\_STAT1\_K562.git

1. cd ChIPseq\_STAT1\_K562  
2. Install dependencies using `requirements.txt`:  
   pip install \-r requirements.txt  
3. Run the main workflow script:  
   bash ChIPseq\_workflow.sh

## **Authors and Acknowledgments**

* **Abhishek Shrestha**   
* **NCBI SRA project (GSE31477)**

## **License**

This project is licensed under the MIT License \- see the LICENSE file for details.

