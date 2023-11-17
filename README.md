# Ashwood-MIA-RNAseq-2023
Osman et al MIA placenta and fetal brain RNAseq analysis

Analysis of 3'Quant RNA Seq of placenta and fetal brain. 
Raw Fastq QC, UMI collapsing, alignment, and gene counts was performed using the Master_QuantSeqAnalyisV2.sh script. 
Reads were then read into R as a counts matrix, filtered, and differential expression analysis between tissues, treatments and sexes was performed using the LimmaVoom_DEG_collapsUMI_9_24_19.R script.


