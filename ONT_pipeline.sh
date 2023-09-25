####--------- ALIGNMENT AND MAPPING OF LONG READ RNA AND DNA ----------######
#Author: Kilian Kleemann
#Date: August 2023
#Resources: 
	# Epi2Me labs
	# Minimap2 - Pomoxis
	# Stringtie
	# Gencode reference

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

#Switiching conda environments
conda activate ONT_pipelines

#Download references #download reference from https://www.gencodegenes.org/human/ Genome sequence (GRCh38.p14)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz 
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
gzip -d *.gz


#pipeline for multiple fastq files 
#Make samplelist of all fastq files
echo 'Detected fastq samples'
find ./fastq_pass -name "*.fastq.gz" -maxdepth 1 -type f -exec basename "{}" \; |  cut -d '.' -f1 | sort -u > sample_list.txt
cat sample_list.txt


mkdir mapping
#Process all files in the directory
cat sample_list.txt | while read sample; do
	minimap2 -ax map-ont <directory>GRCh38.p13.genome.fa <directory>/${sample}.fastq.gz | samtools sort -o mapping/${sample}.sorted.bam
done

#Index all bamfiles 
samtools index -M mapping/*

#making BW files for vizualisation
echo "Making BIGWIG_files"
mkdir BIGWIG_files
cat sample_list.txt | while read sample; do
	bamCoverage -b mapping/${sample}.sorted.bam -o BIGWIG_files/${sample}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads 2000 --numberOfProcessors 4 --skipNonCoveredRegions --samFlagExclude 16
done

#Quantification with StringTie
mkdir stringtie_output
cat sample_list.txt | while read sample; do
	stringtie -e -p 8 \
	-G /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.refGene.gtf \
	-o stringtie_output/${sample}.gtf \
	-L mapping/${sample}.sorted.bam 
done



##### Merging all fastq files - ONT outputs fastq after 4000 reads.
#Checking that all fastq files per biological replicate are in separate folders
#If NOT then merge based on sample list
cat <directory>/*.fastq.gz > concat_output.fastq.gz

#Alignment of merged fastq file
mkdir mapping 
minimap2 -ax map-ont <directory>GRCh38.p13.genome.fa <directory>/concat_output.fastq.gz | samtools sort -o mapping/concat_output.sorted.bam


mkdir BIGWIG_files
bamCoverage -b mapping/concat_output.sorted.bam -o BIGWIG_files/concat_output.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads 2000 --numberOfProcessors 4 --skipNonCoveredRegions --samFlagExclude 16


mkdir stringtie_output
stringtie -e -p 8 \
	-G /Volumes/OhneTitel/Melani_THP1_August_2023/reference/GRCh38.refGene.gtf \
	-o stringtie_output/concat_output.gtf \
	-L mapping/concat_output.sorted.bam 


#Extracting counts using python script