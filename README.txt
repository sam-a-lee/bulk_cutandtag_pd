### OVERVIEW ###

These scripts analysis bulk cell type CUT&Tag data that was generated through FANS of striatum from indiviudals with (N=5) and without (N=5) Parkinson's disease.

Three different cell types were isolated using FANS: neurons (NeurN+), microglia (IRF8+), and oligodendrocytes (Sox10+).

A total of 100000 nuclei of each cell type were used to perform nanobody CUT&Tag for H3K27ac. 

CUT&Tag libraries were sequenced to a target depth of 20 million reads per sample at 50bp paired end sequencing on an Illumina Nextseq2000 sequencer by the Imperial College London BRC Genomics facility. 

In total, 28 libraries were analysed. All three samples were available for all samples, except for PD833, for which only microglia were available.

Libraries were run in technical duplicate (lanes one and two).


### DATA STUCTURE ###

Reads are dividing into lane one (/1) and lane two (/2). 

Reads are already demultiplexed.

i5 and i7 adapters are provided in files following the format SAMPLE_S1_L001_I1_001.fastq.gz or SAMPLE_S1_L001_I2_001.fastq.gz. These can be ignored during processing. 

Actual reads are stored in files following the format SAMPLE_S1_L001_R1_001.fastq.gz and SAMPLE_S1_L001_R2_001.fastq.gz



### PREPROCESSING ###

Quality of technical duplicates was examined with fastqc prior to merging fastq files. 