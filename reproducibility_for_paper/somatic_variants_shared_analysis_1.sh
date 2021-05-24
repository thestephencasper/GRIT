#! /bin/sh

#$ -v ROSENLAB_BASE,SOFTWARE,PROJECTS,DATA
#$ -N wgs_total_variants_shared
#$ -j yes
#$ -o wgs_total_variants_shared.log
#$ -l h_vmem=8G
#$ -l h_rt=2:00:00

# (for hu in {5,7,8,9,10}; do sample="colony_${hu}.bam"; [ -d "your/file/path" ] && qsub -cwd your/script/path/somatic_variants_shared_analysis_1.sh $sample; done)

infile="$(basename $1 .bam)"
echo "Processing sample ${infile}"

work_dir="your/file/path"
cd "$work_dir" || {
	ls -l "$work_dir"
	exit
}

#mkdir shared


################################################################################
#####STEP1: IDENTIFY THE OVERLAPPED SNVS AND INDELS BETWEEN MUTECT2 AND STRELKA2
export PATH="your/software/path/bedtools2-2.29.2/bin${PATH:+:${PATH}}"
bedtools intersect   -wa -wb\
                     -a mutect2/${infile}_mutect2_snvs_filtered.vcf \
                     -b strelka2/${infile}/results/variants/${infile}_strelka2_snvs_filtered.vcf \
                     >shared/${infile}_shared_snvs.txt
 echo ${infile}
 wc -l shared/${infile}_shared_snvs.txt

bedtools intersect   -wa -wb\
                     -a mutect2/${infile}_mutect2_indels_filtered.vcf \
                     -b strelka2/${infile}/results/variants/${infile}_strelka2_indels_filtered.vcf \
                     >shared/${infile}_shared_indels.txt

echo ${infile}
wc -l shared/${infile}_shared_indels.txt


################################################################################
#####STEP2: IDENTIFY ALTERED ALLELE FREQUENCY FROM MUTECT2 AND STRELKA2
##Extract Mutect and Strelka variants frequency
#Mutect has AF, which is slightly different from alt_cnts/total_cnts reference1: https://github.com/broadinstitute/gatk/issues/6067 reference2: https://github.com/broadinstitute/gatk/issues/7016
#Strelka2 doesn't directly provide ref counts and alt counts, needs to calculate based on tier 1 A C G T counts
cat shared/${infile}_shared_snvs.txt|awk 'BEGIN { OFS = "\t" }; {split($22, arr, "[:,]"); split($11, sp, "[:]"); print $1,$2,$4,$5,sp[3],arr[5],arr[7],arr[9],arr[11]}' |\
awk 'BEGIN{print "chr\tpos\tref\talt\tmutect2_af\tA\tC\tG\tT"}1'>shared/${infile}_shared_snvs_final.txt

cat shared/${infile}_shared_indels.txt|awk 'BEGIN { OFS = "\t" }; {split($22, arr, "[|:,-]"); split($11, sp, "[:]"); print $1,$2,$4,$5,sp[3],arr[3],arr[5]}' |\
awk 'BEGIN{print "chr\tpos\tref\talt\tmutect2_af\tref_cnts\talt_cnts"}1'>shared/${infile}_shared_indels_final.txt


################################################################################
#####STEP3: ANNOTATE VARIANTS BASED ON FUNCTIONAL ANNOTATION
#Download required database with correct genome build 
#All available database can be found here https://annovar.openbioinformatics.org/en/latest/user-guide/download/#-for-gene-based-annotation
annovar="your/software/path/annovar"
$annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene $annovar/humandb/
cat $annovar/humandb/wgEncodeOpenChromDnaseHek293tPk.narrowPeak|awk 'BEGIN { OFS = "\t" }; {print $0}' |cut -f 1,2,3,5 >$annovar/humandb/hg19_wgEncodeOpenChromDnaseHek293tPk.txt
##Liftover to $annovar/humandb/hg38_liftover_wgEncodeOpenChromDnaseHek293tPk.txt


#Generate a table output of refGene (gene fitler) and 293T DNase Hypersensitive region (region filter)
$annovar/table_annovar.pl wgs_somatic_total_snvs_mutect2_strelka2.avinput $annovar/humandb/ -buildver hg38 -out wgs_somatic_total_snvs_annovar -remove -protocol refGene,liftover_wgEncodeOpenChromDnaseHek293tPk -operation g,r -nastring . -csvout -polish -xref $annovar/example/gene_xref.txt
$annovar/table_annovar.pl wgs_somatic_total_indels_mutect2_strelka2.avinput $annovar/humandb/ -buildver hg38 -out wgs_somatic_total_indels_annovar -remove -protocol refGene,liftover_wgEncodeOpenChromDnaseHek293tPk -operation g,r -nastring . -csvout -polish -xref $annovar/example/gene_xref.txt
