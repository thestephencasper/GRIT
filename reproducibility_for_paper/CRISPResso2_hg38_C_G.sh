#! /bin/sh

#$ -v ROSENLAB_BASE,SOFTWARE,PROJECTS,DATA
#$ -N CRISPResso2_baseC_G_center-15
#$ -t 1:8
#$ -j yes
#$ -o CRISPResso2_baseC_G_center-15.log
#$ -l h_vmem=24G
#$ -l h_rt=3:00:00

##qsub -cwd your/script/path/CRISPResso2_hg38_C_G.sh

work_dir="your/file/path"
cd "$work_dir" || {
	ls -l "$work_dir"
	exit
}

infile="$(find $work_dir/*.bam | sort | head -n ${SGE_TASK_ID} | tail -n 1)"


your/software/path/CRISPRessoWGS -b "${infile}" \
			-f CRISPResso2_genelist_input.txt \
			-r hg38.fa \
			-p 5 \
			--exclude_bp_from_left 0 --exclude_bp_from_right 0 \
			--min_reads_to_use_region 2 --plot_window_size 15 \
			--name $(basename ${infile})_len41_hg38_baseC \
			--base_edit \
			--conversion_nuc_from C\
			-wc -15 -w 10 \
			--base_editor_output

your/software/path/CRISPRessoWGS -b "${infile}" \
			-f CRISPResso2_genelist_input.txt \
			-r hg38.fa \
			-p 5 \
			--exclude_bp_from_left 0 --exclude_bp_from_right 0 \
			--min_reads_to_use_region 2 --plot_window_size 15 \
			--name $(basename ${infile})_len41_hg38_baseG \
			--base_edit \
			--conversion_nuc_from G\
			-wc -15 -w 10 \
			--base_editor_output