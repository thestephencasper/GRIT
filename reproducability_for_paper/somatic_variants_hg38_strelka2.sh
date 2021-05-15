#! /bin/sh

#$ -v ROSENLAB_BASE,SOFTWARE,PROJECTS,DATA
#$ -N wgs_off_target_strelka2
#$ -t 1:8
#$ -j yes
#$ -o wgs_off_target_strelka2.log
#$ -l h_vmem=128G
#$ -l h_rt=120:00:00

###qsub -cwd your/script/path/somatic_variants_hg38_strelka2.sh


work_dir="your/file/path"
cd "$work_dir" || {
	ls -l "$work_dir"
	exit
}

#mkdir strelka2

infile="$(find $work_dir/colony_*.bam | sort | head -n ${SGE_TASK_ID} | tail -n 1)"
echo "Processing sample ${infile}"

##Call the somatic variants
strelka2_folder=your/software/path/strelka-2.9.10.centos6_x86_64
python2.7 $strelka2_folder/bin/configureStrelkaSomaticWorkflow.py  \
	--normalBam colony_1.bam \
    --tumorBam ${infile} \
    --referenceFasta hg38.fa \
    --runDir strelka2/$(basename ${infile} .bam)

##execution on a single local machine with 30 parallel jobs
python2.7 strelka2/$(basename ${infile} .bam)/runWorkflow.py -m local -j 30


#VCF for 1)SNPs 2)passed the filter only and 3)in chromosome 1-22,X,Y 
zcat strelka2/$(basename ${infile} .bam)/results/variants/somatic.snvs.vcf.gz | grep -w '^#\|^#CHROM\|chr[1-9]\|chr[1-2][0-9]\|chr[X]\|chr[Y]' |\
awk '{if($0 ~ /\#/) print; else if($7 == "PASS" && length($4)==1 && length($5)==1) print}' > strelka2/$(basename ${infile} .bam)/results/variants/$(basename ${infile} .bam)_strelka2_snvs_filtered.vcf

#VCF for 1)Indels 2)passed the filter only and 3)in chromosome 1-22,X,Y 
zcat strelka2/$(basename ${infile} .bam)/results/variants/somatic.indels.vcf.gz | grep -w '^#\|^#CHROM\|chr[1-9]\|chr[1-2][0-9]\|chr[X]\|chr[Y]' |\
awk '{if($0 ~ /\#/) print; else if($7 == "PASS" && (length($4)!=1 || length($5)!=1)) print}' > strelka2/$(basename ${infile} .bam)/results/variants/$(basename ${infile} .bam)_strelka2_indels_filtered.vcf


#colony_4 test
#zcat strelka2/colony_4/results/variants/somatic.snvs.vcf.gz | grep -w '^#\|^#CHROM\|chr[1-9]\|chr[1-2][0-9]\|chr[X]\|chr[Y]' |\
#awk '{if($0 ~ /\#/) print; else if($7 == "PASS" && length($4)==1 && length($5)==1) print}' > strelka2/colony_4/results/variants/colony_4_strelka2_snvs_filtered.vcf


#/broad/rosenlab_archive/anlu/CRISPResso2/data/wgs/strelka2/colony_4/results/variants/.snapshot