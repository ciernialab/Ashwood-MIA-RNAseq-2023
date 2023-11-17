#!/bin/bash


#######################################################################################
#RNA QuantSeq Analysis: https://www.lexogen.com/quantseq-data-analysis/
#9/23/2019
#######################################################################################
#######################################################################################

#make poly A tail file to filter out 
#printf ">polyA\nAAAAAAAAAAAAA\n>polyT\nTTTTTTTTTTTTT\n" | gzip - >  polyA.fa.gz

#make sample files.txt for each sample
#from within fastq file folder: ls -R *fastq.gz > fastq.txt
#trim with:
#cat fastq.txt | awk -F '_' '{print $1}' > samples.txt
#remove undetermined manually

#######################################################################################
#Trim and UMI Index
#######################################################################################
#intall packages if don't have them:
#install fastqc with brew or conda
#install multiqc with: pip install multiqc
#From Lexogen: The programs umi2index and collapse_UMI_bam are distributed as binaries compiled on Ubuntu 16.04, CentOS 6.6, CentOS 7, and MACos. The HTS library from samtools (http://www.htslib.org/download/) is required by collapse_UMI_bam. After installation, the HTS library should be added to the standard library search path. Alternatively, environment variable LD_LIBRARY_PATH can be set to point to the location of the HTS library.
#from bbmap:https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/ download bbduk.sh
#######################################################################################
#make folders
#path to raw fastq files: /Users/aciernia/Desktop/MIA_RNAseq/raw_sequences/Project_PAAC_L8_H1878P_Ciernia/
#path to analysis: /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/
#path to output file: /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/
#scripts: /Users/aciernia/Desktop/MIA_RNAseq/scripts/
#run code from in /Users/aciernia/Desktop/MIA_RNAseq/

mkdir -p /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/fastqc_pretrim
mkdir -p /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/fastqc_posttrim
mkdir -p /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/trimmed_sequences
#run trim.sh srcipt:
#./trim.sh

#it contains this:
for sample in `cat /Users/aciernia/Desktop/MIA_RNAseq/samples.txt`
do
R1=${sample}
echo ${R1} "pre qc"

#fastqc on pre-trimmed file
fastqc /Users/aciernia/Desktop/MIA_RNAseq/raw_sequences/Project_PAAC_L8_H1878P_Ciernia/${R1}*.fastq.gz --outdir /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/fastqc_pretrim

#index UMIs > add UMI sequence for each read to the read identifier
#The program umi2index pre-processes fastq files by extracting the UMI from a read and storing it as a read index. This step is performed right after demultiplexing to ensure that artificial UMI sequences do not influence subsequent analysis of endogenous read sequences. The program is invoked as follows:
#umi2index <fastq_gz_file_in> <fastq_gz_file_out>
echo ${R1} "umi2index"
./umi2index /Users/aciernia/Desktop/MIA_RNAseq/raw_sequences/Project_PAAC_L8_H1878P_Ciernia/${R1}*.fastq.gz /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/${R1}.UMI.fastq.gz

### remove the adapter contamination, polyA read through, and low quality tails
##https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
#In ktrim=r mode, once a reference kmer is matched in a read, that kmer and all the bases to the right will be trimmed, leaving only the bases to the left; this is the normal mode for adapter trimming
## quality-trim to Q10 using the Phred algorithm,
#set Java's memory usage to 24g with -Xmx24g
echo ${R1} "trimming"

/Users/aciernia/Desktop/programs/bbmap/bbduk.sh -Xmx24g in=/Users/aciernia/Desktop/MIA_RNAseq/analysis_files/${R1}.UMI.fastq.gz out=/Users/aciernia/Desktop/MIA_RNAseq/analysis_files/${R1}.trim.fastq.gz ref=/Users/aciernia/Desktop/scripts/polyA.fa.gz,/Users/aciernia/Desktop/scripts/truseq_rna.fa.gz k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 threads=12


#fastqc on post-trimmed file
echo ${R1} "post qc"
fastqc /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/${R1}.trim.fastq.gz --outdir /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/fastqc_posttrim


done

#######################################################################################
#collect all qc together with:

multiqc /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/fastqc_pretrim/ --filename /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/PreTrim_multiqc_report.html --ignore-samples Undetermined* --interactive

multiqc /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/fastqc_posttrim/ --filename /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/PostTrim_multiqc_report.html --ignore-samples Undetermined* --interactive


#move to trimmed sequences folder
mv /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/*.trim.fastq.gz /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/trimmed_sequences
#######################################################################################
##align to mm10 with STAR
######################################################################################
#install STAR
#https://github.com/alexdobin/STAR
#git clone https://github.com/alexdobin/STAR.git
#brew install gcc
#Build STAR:
# note that the path to c++ executable has to be adjusted to its current version > 9
#cd source
#make STARforMacStatic CXX=/usr/local/Cellar/gcc/9.1.0/bin/g++-9
#run with ./STAR -h
#or add to path: export PATH=$PATH:/the_dir_with_STAR
#export PATH=$PATH:/Users/aciernia/Desktop/programs/STAR/source
#echo $PATH
#echo 'export PATH=/Users/aciernia/Desktop/programs/STAR/source:$PATH' >>~/.bash_profile
#source .bash_profile
# build indicies with mm10STARbuild.sh
#######################################################################################

#run trim.sh srcipt:
#./STAR_alignmm10.sh
mkdir -p /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out

#it contains this:
for sample in `cat /Users/aciernia/Desktop/MIA_RNAseq/samples.txt`
do
R1=${sample}

echo ${R1} "unzipping"

gunzip /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/trimmed_sequences/${R1}.trim.fastq.gz

echo ${R1} "mapping started"

#allows 10000 files to be open at once:
ulimit -n 10000

STAR --runThreadN 12 --genomeDir /Users/aciernia/Desktop/programs/STAR_libs/mm10/star_indices/ --readFilesIn  /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/trimmed_sequences/${R1}.trim.fastq --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}

#rezip fastq
gzip /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/trimmed_sequences/${R1}.trim.fastq

echo ${R1} "mapping completed"
done


#######################################################################################
#mapping qc
#takes Log.final.out files and compiles into HTML

multiqc /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/*Log.final.out --filename /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/STARAlignmentLogs.html --interactive


#######################################################################################
#UMI Filter after alignment
#######################################################################################
#From Lexogen: The programs umi2index and collapse_UMI_bam are distributed as binaries compiled on Ubuntu 16.04, CentOS 6.6, CentOS 7, and MACos. The HTS library from samtools (http://www.htslib.org/download/) is required by collapse_UMI_bam. After installation, the HTS library should be added to the standard library search path. Alternatively, environment variable LD_LIBRARY_PATH can be set to point to the location of the HTS library.

#######################################################################################
#make folders

#it contains this:
for sample in `cat /Users/aciernia/Desktop/MIA_RNAseq/samples.txt`
do
R1=${sample}
echo ${R1} "index bam"

#Indexed bam files are necessary for many visualization and downstream analysis tools

samtools index /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}Aligned.sortedByCoord.out.bam /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}Aligned.sortedByCoord.out.bai

echo ${R1} "flagstat on unfiltered bam"
#flagstat on reads
samtools flagstat /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}Aligned.sortedByCoord.out.bam > /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}.preUMIcollapse.flagstat.txt

#collapse UMIs
#Program collapse_UMI_bam takes a bam file as input which is generated by aligning a fastq file to a reference genome. This fastq file is derived from <fastq_gz_file_out>, usually by adapter and quality trimming. The output bam file is used for downstream analysis. This version of collapse_UMI_bam has been tested with bam files generated by the star aligner
#collapse_UMI_bam <bam_file_in> <bam_file_out>

echo ${R1} "umi2index"
./collapse_UMI_bam /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}Aligned.sortedByCoord.out.bam /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}.collapseUMI.bam


echo ${R1} "coordinate sort and index"
samtools sort  /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}.collapseUMI.bam -o /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}.collapseUMI.sort.bam

samtools index /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}.collapseUMI.sort.bam /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}.collapseUMI.sort.bai

#remove unmapped reads
#samtools view -h -F 4 -b ${R1}.collapseUMI.sort.bam > ${R1}.collapseUMI.mapped.bam

echo ${R1} "flagstat on filtered bam"
#flagstat on reads
samtools flagstat /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}.collapseUMI.sort.bam > /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}.collapseUMI.flagstat.txt


done

#######################################################################################
#collect all qc flagstats together with:

multiqc /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/*.flagstat.txt --filename /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/STARAlignment.html --interactive

#######################################################################################
#additional alignment QC with Qualimap
#######################################################################################
#intall packages if don't have them:
#install qualimap
#add the location of the Qualimap tool to our PATH variable
#echo 'export PATH=/Users/aciernia/Desktop/programs/qualimap_v2.2.1/:$PATH' >>~/.bash_profile
#source ~/.bash_profile
#qualimap rnaseq --help
#######################################################################################
#make folders
mkdir -p /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/qualimapQC_preUMI
mkdir -p /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/qualimapQC_postUMI

#By default, Qualimap will try to open a GUI to run Qualimap, so we need to run the unset DISPLAY command
unset DISPLAY



#it contains this:
for sample in `cat /Users/aciernia/Desktop/MIA_RNAseq/samples.txt`
do
R1=${sample}
echo ${R1} "qualimapQC RNAseq"

##-outdir: output directory for html report
#-a: Counting algorithm - uniquely-mapped-reads(default) or proportional (each multi-mapped read is weighted according to the number of mapped locations)
#-bam: path/to/bam/file(s)
#-p: Sequencing library protocol - strand-specific-forward, strand-specific-reverse or non-strand-specific (default)
#-gtf: path/to/gtf/file - needs to match the genome build and GTF used in alignment
#--java-mem-size=: set Java memory

#Note that Qualimap must be run with the -outdir option as well as -outformat HTML (which is on by default). MultiQC uses files found within the raw_data_qualimapReport folder (as well as genome_results.txt).

qualimap rnaseq -bam /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}Aligned.sortedByCoord.out.bam -gtf /Users/aciernia/Desktop/programs/STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf -outdir /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/qualimapQC_preUMI/${R1}.preUMI --java-mem-size=8G -p strand-specific-forward

qualimap rnaseq -bam /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}.collapseUMI.sort.bam -gtf /Users/aciernia/Desktop/programs/STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf -outdir /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/qualimapQC_postUMI/${R1}.postUMI --java-mem-size=8G -p strand-specific-forward


done

#######################################################################################
#collect all qc together with:

multiqc /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/qualimapQC_preUMI/*.preUMI/ --filename /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/PreUMI_multiQualimapRNAseq.html --interactive
multiqc /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/qualimapQC_postUMI/*.postUMI/ --filename /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/PostUMI_multiQualimapRNAseq.html --interactive

#######################################################################################
#count with subread feature counts
#######################################################################################
#install if don't have:
#http://bioinf.wehi.edu.au/featureCounts/
#download binary and add to path
#echo 'export PATH=/Users/aciernia/Desktop/programs/subread-1.6.5-MacOSX-x86_64/bin:$PATH' >>~/.bash_profile
#source ~/.bash_profile
#######################################################################################

#counting

mkdir -p /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/counts


#it contains this:
for sample in `cat /Users/aciernia/Desktop/MIA_RNAseq/samples.txt`
do
R1=${sample}


#quant seq kit is FWD stranded

echo ${R1} "count started"
#0 (unstranded), 1 (stranded) and 2 (reversely stranded)
# Number of CPU threads
#set for gene_id > ensembl id mm10
#exon level counts

featureCounts -T 12 -s 1 -t exon -g gene_id -a /Users/aciernia/Desktop/programs/STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf -o /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/counts/${R1}.preUMI.counts.txt /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}Aligned.sortedByCoord.out.bam

featureCounts -T 12 -s 1 -t exon -g gene_id -a /Users/aciernia/Desktop/programs/STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf -o /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/counts/${R1}.collapsedUMI.counts.txt /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/star_out/${R1}.collapseUMI.sort.bam

echo ${R1} "count completed"
done

#######################################################################################
#collect all qc together with:
multiqc /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/counts/*.preUMI.counts.txt.summary --filename /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/Counts.preUMI --interactive
multiqc /Users/aciernia/Desktop/MIA_RNAseq/analysis_files/counts/*.collapsedUMI.counts.txt.summary --filename /Users/aciernia/Desktop/MIA_RNAseq/analysis_output/Counts.collapsedUMI --interactive
#######################################################################################
#######################################################################################
