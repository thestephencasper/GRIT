#! /bin/sh

#$ -v ROSENLAB_BASE,SOFTWARE,PROJECTS,DATA
#$ -N wgs_total_variants_mutect2
#$ -t 1:8
#$ -j yes
#$ -o wgs_total_variants_mutect2.log
#$ -l h_vmem=32G
#$ -l h_rt=120:00:00

###qsub -cwd your/script/path/somatic_variants_hg38_mutect2.sh


work_dir="your/file/path"
cd "$work_dir" || {
	ls -l "$work_dir"
	exit
}

#mkdir mutect2

infile="$(find $work_dir/colony_*.bam | sort | head -n ${SGE_TASK_ID} | tail -n 1)"
echo "Processing sample ${infile}"

gatk_folder=your/software/path/gatk-4.2.0.0

##Call the somatic variants
java -jar $gatk_folder/gatk-package-4.2.0.0-local.jar Mutect2 -I ${infile} -I $work_dir/colony_1.bam \
-tumor $(basename ${infile} .bam) -normal colony_1 -O mutect2/$(basename ${infile} .bam)_unfiltered_Mutect2.vcf \
-R Homo_sapiens_assembly38.fasta

##Annotate the somatic variants
java -jar $gatk_folder/gatk-package-4.2.0.0-local.jar FilterMutectCalls \
-R Homo_sapiens_assembly38.fasta \
-V mutect2/$(basename ${infile} .bam)_unfiltered_Mutect2.vcf \
-O mutect2/$(basename ${infile} .bam)_Mutect2.vcf

##VCF for 1)SNVs 2)passed the filter only and 3)in chromosome 1-22,X,Y 
cat mutect2/$(basename ${infile} .bam)_Mutect2.vcf | grep -w '^#\|^#CHROM\|chr[1-9]\|chr[1-2][0-9]\|chr[X]\|chr[Y]' |\
awk '{if($0 ~ /\#/) print; else if($7 == "PASS" && length($4)==1 && length($5)==1) print}' > mutect2/$(basename ${infile} .bam)_mutect2_snvs_filtered.vcf

##VCF for 1)Indels 2)passed the filter only and 3)in chromosome 1-22,X,Y 
cat mutect2/$(basename ${infile} .bam)_Mutect2.vcf | grep -w '^#\|^#CHROM\|chr[1-9]\|chr[1-2][0-9]\|chr[X]\|chr[Y]' |\
awk '{if($0 ~ /\#/) print; else if($7 == "PASS" && (length($4)!=1 || length($5)!=1)) print}' > mutect2/$(basename ${infile} .bam)_mutect2_indels_filtered.vcf

