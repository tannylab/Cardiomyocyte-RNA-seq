
cd /path/to/experiment/data

##Trim adapter sequences, remove low quality reads and set length cutoff
for sample in `ls *.fastq.gz`
do
	base=$(basename $sample ".fastq.gz")
	/path/to/TrimGalore-0.6.0/trim_galore --phred33 --length 36 -q 5 --stringency 1 -e 0.1 \
	--path_to_cutadapt /path/to/cutadapt $sample
done	

##Move trimmed fastq files to a new directory
mkdir trimmed
mv *R1_trimmed.fq.gz trimmed
mv *trimming_report.txt trimmed

##Run FastQC on all samples to assess quality of reads and if further processing is required
##Make directory for FastQC files
mkdir trimmed/FastQC

for sample in `ls trimmed/*R1_trimmed.fq.gz`
do
	/path/to/fastqc $sample -o trimmed/FastQC
done

##Generate genome using Ensembl Rattus_norvegicus.Rnor_6.0.94 genome
/path/to/STAR-2.6.0a/bin/Linux_x86_64/STAR \
--runThreadN 12 --runMode genomeGenerate --genomeDir /path/to/rat_star_genome \
--genomeFastaFiles 	/path/to/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa

#Make direectory for aligned BAM files
mkdir aligned_bam

#Align reads to Ensembl rat genome with STAR
for sample in `ls trimmed/*R1_trimmed.fq.gz`
do
	base=$(basename $sample "_R1_trimmed.fq.gz")
	/path/to/STAR.2.7.1a \
	--runThreadN 12 --genomeDir /path/to/rat_star_genome --sjdbGTFfile /path/to/Rattus_norvegicus.Rnor_6.0.94.gtf --sjdbOverhang 49 \
	--readFilesIn $sample --readFilesCommand zcat --outFileNamePrefix path/to/aligned_bam/${base} \
	--outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 35000000000
done

##Assemble transcripts with StringTie and Ensembl Rattus_norvegicus.Rnor_6.0.94 genome
##Make directory for StringTie results
mkdir StringTie_Ensembl

for sample in `ls aligned_bam/*.bam`
do
	base=$(basename $sample "Aligned.sortedByCoord.out.bam")
	/path/to/stringtie/1.3.4d -e -B -p 4 --rf -G /path/to/Rattus_norvegicus.Rnor_6.0.94.gtf \
	-o path/to/StringTie_Ensembl/${base}/${base}.gtf $sample
done

# Rest of analysis is completed with R script


