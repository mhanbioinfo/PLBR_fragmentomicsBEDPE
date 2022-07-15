#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -t 1-00:00:00
#SBATCH -D ./ 
#SBATCH --mem=32G
#SBATCH -J bam_picardInsertSize
#SBATCH -p himem
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -o %j-%x.out
#SBATCH -e %j-%x.err

echo "Job started at "$(date) 
time1=$(date +%s)

module load samtools/1.14
module load bedtools/2.27.1
module load java/8
module load picard/2.10.9
module load perl
module load R/4.1.0

## --------------------------------------------------------- ##

INPUT_DIR="/cluster/projects/pughlab/projects/fragmentomics_bedpe_input/original_files"
#INPUT_F="TGL48_0011_Ct_T_PE_325_WG.merged.sorted.bam"
INPUT_F="TGL48_0011_Ct_T_PE_325_WG.merged.sorted.sub1pct.coordSortd.bam"
OUT_DIR="/cluster/projects/pughlab/projects/fragmentomics_bedpe_input/insert_sizes/try14_picardInsertSize_sub1pctBAM_alltogether"
#OUT_DIR="/cluster/projects/pughlab/projects/fragmentomics_bedpe_input/insert_sizes/try16_picardInsertSize_fullBAM_alltogether"

## picard CollectInsertSizeMetrics
echo "Performing picard CollectInsertSizeMetrics... "
java -jar /cluster/tools/software/picard/2.10.9/picard.jar \
    CollectInsertSizeMetrics \
    INCLUDE_DUPLICATES=true \
    I=${INPUT_DIR}/${INPUT_F} \
    O=${OUT_DIR}/${INPUT_F%.*}_picardInsertSize_wDups.txt \
    H=${OUT_DIR}/${INPUT_F%.*}_picardInserSizeHisto_wDups.pdf \
    M=0 W=600

#echo "Performing transpose on picard output... "
#cat ${OUT_DIR}/${INPUT_F%.*}_picardInsertSize_wDups.txt | \
#    grep -v "^#" | head -5 | \
#    perl ./transpose.pl - > ${OUT_DIR}/${INPUT_F%.*}_picardInsertSize_wDups_transposed.txt


time2=$(date +%s)
echo "Job ended at "$(date) 
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"

## EOF
