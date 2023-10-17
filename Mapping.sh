####--------- ALIGNMENT AND MAPPING OF LONG READ RNA AND DNA ----------######
#Author: Kilian Kleemann
#Date: August 2023
#Resources: 
	# Epi2Me labs
	# Minimap2 - Pomoxis
	# Stringtie
	# Gencode reference


#Installation of anaconda3
#kilian@kilian-aurora-R15:~$ export PATH="$PATH:/home/kilian/anaconda3/bin/"
#kilian@kilian-aurora-R15:~$ source ~/.profile 
#kilian@kilian-aurora-R15:~$ echo $PATH

#Issue with x64 M1 chip !! do not copy the notes into commad line
CONDA_SUBDIR=osx-64 conda create -n ONT_pipeline   # create a new environment
conda activate ONT_pipeline
conda env config vars set CONDA_SUBDIR=osx-64
conda deactivate
conda activate ONT_pipeline
conda install pomoxis
conda install samtools
conda install stringtie
conda install gffcompare


#Setting environment in linux
CONDA_SUBDIR=linux-64 conda create -n ONT_pipeline   # create a new environment
conda activate ONT_pipeline
conda env config vars set CONDA_SUBDIR=linux-64
conda deactivate
conda activate ONT_pipeline
conda install pomoxis
conda install samtools
conda install stringtie
conda install gffcompare

#Switiching conda environments
conda activate ONT_pipelines

#Download references #download reference from https://www.gencodegenes.org/human/ Genome sequence (GRCh38.p14)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz 
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
gzip -d *.gz


#single file to try
#Standard
minimap2 -a ref.mmi 230724_S1_RNA_flongle/fastq_pass/APT278_pass_d1f70bd3_d535a438_1.fastq.gz > APT278_pass_d1f70bd3_d535a438_1_alignment.sam 
#Oxford nanopore reads
minimap2 -ax splice ref.fa nanopore-cdna.fa > aln.sam        # Nanopore 2D cDNA-seq
minimap2 -ax splice -uf -k14 ref.fa direct-rna.fq > aln.sam  # Nanopore Direct RNA-seq



#Make samplelist of all fastq files
echo 'Detected fastq samples'
find ./fastq_pass -name "*.fastq.gz" -maxdepth 1 -type f -exec basename "{}" \; |  cut -d '.' -f1 | sort -u > sample_list.txt
cat sample_list.txt


mkdir mapping
#Process all files in the directory
cat sample_list.txt | while read sample; do
	minimap2 -ax map-ont /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.p13.genome.fa fastq_pass/${sample}.fastq.gz | samtools sort -o mapping/${sample}.sorted.bam
done


############ Alignment and RAW file processing RNA #############
#Oxford nanopore reads cDNA reads
mkdir cDNA_ANALYSIS/S1_cDNA/mapping_cDNA_specific
minimap2 -ax splice -t 8 reference/GRCh38.p13.genome.fa cDNA_ANALYSIS/S1_cDNA/fastq_pass/concat_output.fastq.gz  > cDNA_ANALYSIS/S1_cDNA/mapping_cDNA_specific/concat_output_ax_splice_cDNA.sam  # Nanopore 2D cDNA-seq

mkdir cDNA_ANALYSIS/WT_cDNA/mapping_cDNA_specific
minimap2 -ax splice -t 8 reference/GRCh38.p13.genome.fa cDNA_ANALYSIS/WT_cDNA/fastq_pass/concat_output.fastq.gz  > cDNA_ANALYSIS/WT_cDNA/mapping_cDNA_specific/concat_output_ax_splice_cDNA.sam  # Nanopore 2D cDNA-seq


#Oxford nanopore reads cDNA reads MULTIMAPPING option
#mkdir mapping_cDNA_specific
#minimap2 -ax splice -N 100 -t 8 reference/GRCh38.p13.genome.fa cDNA_ANALYSIS/S1_cDNA/fastq_pass/concat_output.fastq.gz  > cDNA_ANALYSIS/S1_cDNA/mapping_cDNA_specific/concat_output_ax_splice_cDNA.sam  # Nanopore 2D cDNA-seq

#mkdir mapping_cDNA_specific
#minimap2 -ax splice -N 100 -t 8 reference/GRCh38.p13.genome.fa cDNA_ANALYSIS/S1_cDNA/fastq_pass/concat_output.fastq.gz  > cDNA_ANALYSIS/S1_cDNA/mapping_cDNA_specific/concat_output_ax_splice_cDNA.sam  # Nanopore 2D cDNA-seq

#Sorting sam S1 files
samtools sort cDNA_ANALYSIS/WT_cDNA/mapping_cDNA_specific/concat_output_ax_splice_cDNA.sam -o cDNA_ANALYSIS/WT_cDNA/mapping_cDNA_specific/concat_output_ax_splice_cDNA.sorted.bam

#Sorting sam WT files
samtools sort cDNA_ANALYSIS/S1_cDNA/mapping_cDNA_specific/concat_output_ax_splice_cDNA.sam -o cDNA_ANALYSIS/S1_cDNA/mapping_cDNA_specific/concat_output_ax_splice_cDNA.sorted.bam


#Index all bamfiles 
samtools index --threads 4 -M RNA_ANALYSIS/WT_RNA/mapping_RNA_specific/*.bam
samtools index --threads 4 -M RNA_ANALYSIS/S1_RNA/mapping_RNA_specific/*.bam

#Show quality of alignments
samtools flagstat RNA_ANALYSIS/WT_RNA/mapping_RNA_specific/*.bam > RNA_ANALYSIS/S1_RNA/mapping_RNA_specific/minimap2_alignment_summary.txt 
samtools flagstat RNA_ANALYSIS/S1_RNA/mapping_RNA_specific/*.bam > RNA_ANALYSIS/S1_RNA/mapping_RNA_specific/minimap2_alignment_summary.txt 






############ Alignment and RAW file processing RNA #############
#Oxford nanopore reads direct RNA reads
mkdir RNA_ANALYSIS/S1_RNA/mapping_RNA_specific
minimap2 -ax splice -uf -k14 -t 8 reference/GRCh38.p13.genome.fa RNA_ANALYSIS/S1_RNA/fastq_pass/concat_output.fastq.gz > RNA_ANALYSIS/S1_RNA/mapping_RNA_specific/concat_output_ax_splice_directRNA.sam  # Nanopore Direct RNA-seq

mkdir RNA_ANALYSIS/WT_RNA/mapping_RNA_specific
minimap2 -ax splice -uf -k14 -t 8 reference/GRCh38.p13.genome.fa RNA_ANALYSIS/WT_RNA/fastq_pass/concat_output.fastq.gz > RNA_ANALYSIS/WT_RNA/mapping_RNA_specific/concat_output_ax_splice_directRNA.sam  # Nanopore Direct RNA-seq


#Sorting sam S1 files
samtools sort RNA_ANALYSIS/S1_RNA/mapping_RNA_specific/concat_output_ax_splice_directRNA.sam -o RNA_ANALYSIS/S1_RNA/mapping_RNA_specific/concat_output_ax_splice_directRNA.sorted.bam

#Sorting sam WT files
samtools sort RNA_ANALYSIS/WT_RNA/mapping_RNA_specific/concat_output_ax_splice_directRNA.sam -o RNA_ANALYSIS/WT_RNA/mapping_RNA_specific/concat_output_ax_splice_directRNA.sorted.bam


#Index all bamfiles 
samtools index --threads 4 -M RNA_ANALYSIS/WT_RNA/mapping_RNA_specific/*.bam
samtools index --threads 4 -M RNA_ANALYSIS/S1_RNA/mapping_RNA_specific/*.bam

#Show quality of alignments
samtools flagstat RNA_ANALYSIS/WT_RNA/mapping_RNA_specific/*.bam > RNA_ANALYSIS/S1_RNA/mapping_RNA_specific/minimap2_alignment_summary.txt 
samtools flagstat RNA_ANALYSIS/S1_RNA/mapping_RNA_specific/*.bam > RNA_ANALYSIS/S1_RNA/mapping_RNA_specific/minimap2_alignment_summary.txt 


#Quantification with stringtie / normal GTF
mkdir RNA_ANALYSIS/stringtie_output_RNA_specific
stringtie -e -B -p 8 \
	-G reference/GRCh38.refGene.gtf \
	-o  RNA_ANALYSIS/stringtie_output_RNA_specific/WT_RNA/concat_output.gtf  \
	RNA_ANALYSIS/WT_RNA/mapping_RNA_specific/concat_output_ax_splice_directRNA.sorted.bam 

stringtie -e -B -p 8 \
	-G reference/GRCh38.refGene.gtf \
	-o  RNA_ANALYSIS/stringtie_output_RNA_specific/S1_RNA/concat_output.gtf  \
	RNA_ANALYSIS/S1_RNA/mapping_RNA_specific/concat_output_ax_splice_directRNA.sorted.bam 

#Running Python script to extract raw counts
cd RNA_ANALYSIS/stringtie_output_RNA_specific
python /media/kilian/OhneTitel/Melani_THP1_August_2023/scripts_analysis/prepDE.py3


############ Visualization #############
#making BW files for vizualisation
echo "Making BIGWIG_files"
mkdir BIGWIG_files
cat sample_list.txt | while read sample; do
	bamCoverage -b mapping/${sample}.sorted.bam -o BIGWIG_files/${sample}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads 2000 --numberOfProcessors 4 --skipNonCoveredRegions --samFlagExclude 16
done


mkdir BIGWIG_files
bamCoverage -b mapping/concat_output.sorted.bam -o BIGWIG_files/concat_output.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads 2000 --numberOfProcessors 4 --skipNonCoveredRegions --samFlagExclude 16



#Processing stringtie output!
mkdir stringtie_merged
mkdir stringtie_abundance
cat sample_list.txt | while read sample; do
	stringtie --merge -p 8 -G /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.refGene.gtf -o stringtie_merged/stringtie_merged.gtf stringtie_output/${sample}.gtf -A stringtie_abundance/${sample}.tab
done

#Merging gtf files per experiment
mkdir stringtie_merged
stringtie --merge -p 8 -G /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.refGene.gtf -o stringtie_merged/concat_output.gtf S1_cDNA/stringtie_output/concat_output/concat_output.gtf WT_cDNA/stringtie_output/concat_output/concat_output.gtf 

mkdir stringtie_merged
stringtie --merge -p 8 -G /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.refGene.gtf -o stringtie_merged/concat_output.gtf S1_RNA/stringtie_output/concat_output/concat_output.gtf WT_RNA/stringtie_output/concat_output/concat_output.gtf 


#Comparing the GTF files
mkdir stringtie_merged/gff_compare_check
gffcompare -G -r /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.refGene.gtf -o merged stringtie_merged/stringtie_merged.gtf

#Make counts for DEG analysis
mkdir ballgown
cat sample_list.txt | while read sample; do
	mkdir ballgown/${sample}
	stringtie -e -B -p 8 -G stringtie_merged/stringtie_merged.gtf \
	-o ballgown/${sample}/${sample}.gtf \
	mapping/${sample}.sorted.bam
done

#Quantification with StringTie #DONT ADD THE -L parameter -transcripts not counted properly 
mkdir stringtie_output
cat sample_list.txt | while read sample; do
	stringtie -p 8 \
	-G /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.refGene.gtf \
	-o stringtie_output/${sample}.gtf \
	mapping/${sample}.sorted.bam 
done


#WT-cDNA
mkdir stringtie_output
stringtie -e -B -p 8 \
	-G /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.refGene.gtf \
	-o  stringtie_output/WT_cDNA_concat_output/WT_cDNA_concat_output.gtf  \
	WT_cDNA/mapping/concat_output.sorted.bam 

#S1-cDNA
mkdir stringtie_output
stringtie -e -B -p 8 \
	-G /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.refGene.gtf \
	-o  stringtie_output_S1_cDNA/S1_cDNA_concat_output/S1_cDNA_concat_output.gtf  \
	S1_cDNA/mapping/concat_output.sorted.bam 

#S1-RNA
mkdir stringtie_output
stringtie -e -B -p 8 \
	-G /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.refGene.gtf \
	-o stringtie_output/S1_RNA_concat_output/S1_RNA_concat_output.gtf \
	S1_RNA/mapping/concat_output.sorted.bam 

#WT-RNA
mkdir stringtie_output
stringtie -e -B -p 8 \
	-G /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.refGene.gtf \
	-o stringtie_output/WT_RNA_concat_output/WT_RNA_concat_output.gtf \
	WT_RNA/mapping/concat_output.sorted.bam 



#Processing aligned reads with repeat gtf file
mkdir stringtie_output_rmsk_TE
cat sample_list.txt | while read sample; do
	stringtie -p 8 \
	-G /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38_GENCODE_rmsk_TE.gtf \
	-o stringtie_output_rmsk_TE/${sample}.gtf \
	-L mapping/${sample}.sorted.bam 
done




#Extract counts with featureCounts
mkdir featurecounts
featureCounts -O -T 8 -a stringtie_output/concat_output/concat_output.gtf -o featurecounts/concat_output.txt mapping/concat_output.sorted.bam



#Extracting counts from BAM files and Stringtie GTF (convert to BED) then annotate with HOMER
#Convert Strigntie GTF to BED
convert2bed -i gtf < stringtie_output/concat_output/concat_output.gtf > stringtie_output/concat_output/concat_output.bed
convert2bed -i gtf < reference/GRCh38.refGene.gtf > reference/GRCh38.refGene.bed

#Extracting Data
bedtools multicov -bams S1_cDNA/mapping/*.bam WT_cDNA/mapping/*.bam -bed stringtie_output/concat_output_edit.bed > Extracted_counts_transcripts.txt
#cDNA
bedtools multicov -bams S1_cDNA/mapping/*.bam WT_cDNA/mapping/*.bam -bed /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.refGene.bed > Extracted_counts_transcripts_hg38_ref.txt
#RNA
bedtools multicov -bams S1_RNA/mapping/*.bam WT_RNA/mapping/*.bam -bed /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.refGene.bed > Extracted_counts_transcripts_hg38_ref.txt

#Merging bed files
bedtools merge -i Master_ATAC_peak_merged.bed -d 300 -c 7 -o count > Master_ATAC_peak_merged_d300_c7_all.bed

#Peakannotation with Homer
annotatePeaks.pl Master_ATAC_peak_merged_d300_c7_all.bed mm10 > Annotated_peaks_merged_cDNA.txt


#BigwigMerge for vizualization UCSCtools
bigWigMerge BIGWIG_files/H3K9ac3.bw BIGWIG_files/H3K9ac4.bw BIGWIG_files/H3K9ac5.bw  Merged_bigwig_H3K9ac.bedGraph




#Setting up Epi2Me nextflow analysis
nextflow run epi2me-labs/wf-transcriptomes --help

curl get.nextflow.io | bash
sudo mv nextflow /usr/bin


#FOR DIRECT RNA ANALYSIS 
OUTPUT=~/output;
cat sample_list.txt | while read sample; do
	nextflow run epi2me-labs/wf-transcriptomes \
  --fastq  RNA_ANALYSISfastq_pass/${sample}/${sample}.fastq \
  --ref_genome reference/GRCh38.p13.genome.fa \
  --ref_annotation reference/GRCh38.refGene.gtf \
  --out_dir outdir -w workspace_dir \
  --direct_rna false
done





#Investigation of transposable elements from long-read sequencing
#Required TE concensus sequence in fasta format:
#Converting GTF to bed
gtf2bed reference/GRCh38_GENCODE_rmsk_TE.gtf > reference/GRCh38_GENCODE_rmsk_TE.bed
bedtools getfasta -fi reference/GRCh38.p13.genome.fa -bed reference/hg38_rmsk_ucsd.bed -fo reference/GRCh38.p13.TE.genome.fa

############## TELR ##############
conda install telr 
telr -i mapping/concat_output.sorted.bam  -r reference/GRCh38.p13.TE.genome.fa -l reference/GRCh38.p13.TE.genome.fa 
                

#TALON ANALYSIS
#RNA
mkdir RNA_ANALYSIS/S1_RNA/mapping_TALON
minimap2 -ax splice -uf -k14 -t 4 --MD reference/GRCh38.p13.genome.fa RNA_ANALYSIS/S1_RNA/fastq_pass/concat_output.fastq.gz > RNA_ANALYSIS/S1_RNA/mapping_TALON/concat_output_ax_splice_MD_flag_directRNA.sam  # Nanopore Direct RNA-seq with MD flag

mkdir RNA_ANALYSIS/WT_RNA/mapping_TALON
minimap2 -ax splice -uf -k14 -t 4 --MD reference/GRCh38.p13.genome.fa RNA_ANALYSIS/WT_RNA/fastq_pass/concat_output.fastq.gz > RNA_ANALYSIS/WT_RNA/mapping_TALON/concat_output_ax_splice_MD_flag_directRNA.sam  # Nanopore Direct RNA-seq with MD flag


#labelling reads to identify potential internal priming contructs
mkdir RNA_ANALYSIS/S1_RNA/TALON_labeled
talon_label_reads --f RNA_ANALYSIS/S1_RNA/mapping_TALON/concat_output.sam \
    --g reference/GRCh38.p13.genome.fa  \
    --t 4 \
    --ar 20 \
    --deleteTmp \
    --o RNA_ANALYSIS/S1_RNA/TALON_labeled/concat_output.sam 

mkdir RNA_ANALYSIS/WT_RNA/TALON_labeled
talon_label_reads --f RNA_ANALYSIS/WT_RNA/mapping_TALON/concat_output.sam \
        --g reference/GRCh38.p13.genome.fa  \
        --t 4 \
        --ar 20 \
        --deleteTmp \
        --o RNA_ANALYSIS/WT_RNA/TALON_labeled/concat_output.sam     

talon_reformat_gtf -g /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38_GENCODE_rmsk_TE.gtf > /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38_GENCODE_TALON_rmsk_TE.gtf


talon_initialize_database --f  /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38_GENCODE_rmsk_TE_reformatted.gtf \
  --g hg38_rmsk_ucsd \
  --a hg38 \
  --o hg38 

#--f CONFIG_FILE: Dataset config file: dataset name, sample description, platform, sam file (comma-delimited)  
mkdir TALON_output
talon  --f config.csv \
       --db /Volumes/OhneTitel/Melani_THP1_August_2023/RNA_ANALYSIS/WT_RNA/hg38.db \
       --build hg38_rmsk_ucsd \
       --threads 4 \
       --o TALON_output/





#Extracting counts from full length L1 bed track - did not find anything 
#Extracting Data
#cDNA
bedtools multicov -bams S1_cDNA/mapping/*.bam WT_cDNA/mapping/*.bam -bed /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38_Human_L1_full.bed > Extracted_counts_L1_full_length_transcripts.txt
#RNA
bedtools multicov -bams S1_RNA/mapping/*.bam WT_RNA/mapping/*.bam -bed /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38_Human_L1_full.bed > Extracted_counts_L1_full_length_transcripts_RNA.txt











