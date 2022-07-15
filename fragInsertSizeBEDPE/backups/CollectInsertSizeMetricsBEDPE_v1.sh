#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -t 1-00:00:00
#SBATCH -D ./ 
#SBATCH --mem=32G
#SBATCH -J CollectInsertSizeMetricsBEDPE
#SBATCH -p himem
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -o %j-%x.out
#SBATCH -e %j-%x.err

echo "Job started at "$(date) 
time1=$(date +%s)

INPUT_DIR="/cluster/projects/pughlab/projects/fragmentomics_bedpe_input/original_files"
INPUT_BEDPE="TGL48_0011_Ct_T_PE_325_WG.merged.sorted.sub1pct.coordSortd_coordSortd.bedpe.gz"

OUT_DIR="/cluster/projects/pughlab/projects/fragmentomics_bedpe_input/insert_sizes/try18_CollectInsertSizeMetricsBEDPE_v1"

BEDPE_F2316_TLENgt0="${INPUT_BEDPE%.bedpe.gz}.F2316_TLENgt0.bedpe"
BEDPE_F2316_TLENgt0_TANDEMS="${BEDPE_F2316_TLENgt0%.*}.TANDEMS.bedpe"
BEDPE_F2316_TLENgt0_NOTANDEMS="${BEDPE_F2316_TLENgt0%.*}.NOTANDEMS.bedpe"
BEDPE_F2316_TLENgt0_TANDEMS_COL16="${BEDPE_F2316_TLENgt0_TANDEMS%.*}.col16"

BEDPE_F2316_TLENgt0_NOTANDEMS_RFs="${BEDPE_F2316_TLENgt0_NOTANDEMS%.*}.RFs.bedpe"
BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_COL16="${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs%.*}.col16"

BEDPE_F2316_TLENgt0_NOTANDEMS_FRs="${BEDPE_F2316_TLENgt0_NOTANDEMS%.*}.FRs.bedpe"
BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16="${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs%.*}.col16"



BEDPE_F2316_TANDEMS_INSERTSIZE=""

## remove unmapped/singletons (4, 8), secondary (256), supplementary (2048) reads, or if TLEN == 0
zcat ${INPUT_DIR}/${INPUT_BEDPE} \
    | awk 'BEGIN {OFS="\t"} ($16 != 0 && $17 != 0) {print}' \
    | awk 'BEGIN {OFS="\t"} ((!and($14, 0x4)) && (!and($15, 0x4))) {print}' \
    | awk 'BEGIN {OFS="\t"} ((!and($14, 0x8)) && (!and($15, 0x8))) {print}' \
    | awk 'BEGIN {OFS="\t"} ((!and($14, 0x100)) && (!and($15, 0x100))) {print}' \
    | awk 'BEGIN {OFS="\t"} ((!and($14, 0x800)) && (!and($15, 0x800))) {print}' \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0}

## extract tandem reads (16, 32 both set or unset)
cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0} \
    | awk 'BEGIN {OFS="\t"} ((!and($14, 0x10)) && (!and($14, 0x20))) || ((!and($15, 0x10)) && (!and($15, 0x20))) || ((and($14, 0x10)) && (and($14, 0x20))) || ((and($15, 0x10)) && (and($15, 0x20))) {print}' \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS}

## remove tandem reads (16, 32 both set or unset)
cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0} \
    | awk 'BEGIN {OFS="\t"} ((!and($14, 0x10)) && (!and($14, 0x20))) {next} {print}' \
    | awk 'BEGIN {OFS="\t"} ((!and($15, 0x10)) && (!and($15, 0x20))) {next} {print}' \
    | awk 'BEGIN {OFS="\t"} ((and($14, 0x10)) && (and($14, 0x20))) {next} {print}' \
    | awk 'BEGIN {OFS="\t"} ((and($15, 0x10)) && (and($15, 0x20))) {next} {print}' \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS}


## TANDEM READS ####################################################
PAIR_ORIENTATION="TANDEM"

## get InsertSizes for tandem reads
echo "Getting InsertSizes metrics for TANDEM reads"
cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS} \
    | awk 'BEGIN {OFS="\t"} {print sqrt($16*$16)}' \
    | sort -n \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS_COL16}
#    | awk 'function SUB(F) {sub("^-","",$F)} SUB(16) SUB(17); {print $16, $17}' \
#    | awk '$1 != $2 {print}' \
## TLEN columns 16 and 17 always same number, just different signs

## READ_PAIRS
TANDEM_READ_PAIRS=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS_COL16} | wc -l)
echo $TANDEM_READ_PAIRS

## MEDIAN_INSERT_SIZE (exclude 0's)
TANDEM_MEDIAN=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS_COL16} \
    | awk '{ a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) printf "%.0f\n", (a[x-1]+a[x])/2; else printf "%.0f\n", a[x-1]; }') 
echo $TANDEM_MEDIAN

## MEDIAN_ABSOLUTE_DEVIATION ( abs value of deviations from median, median )
TANDEM_MAD=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS_COL16} \
    | awk -v TANDEM_MEDIAN=$TANDEM_MEDIAN '{$2=$1-TANDEM_MEDIAN; sub("^-","",$2); print $2}' \
    | sort -n \
    | awk '{ a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) printf "%.0f\n", (a[x-1]+a[x])/2; else printf "%.0f\n", a[x-1]; }')
echo $TANDEM_MAD

## MIN
TANDEM_MIN=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS_COL16} \
    | sed -n '1p')
echo $TANDEM_MIN

## MAX
TANDEM_MAX=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS_COL16} \
    | sed -n '$p')
echo $TANDEM_MAX

## MEAN_INSERT_SIZE ( core inserts only : +/-10 MADs around median insert size )
## 10 is default (https://github.com/broadinstitute/picard/blob/949d7f954f84e02e833be66e8f65553a0918de65/src/main/java/picard/analysis/CollectInsertSizeMetrics.java)
TANDEM_CORE_MIN=$(expr $TANDEM_MEDIAN - 10 \* $TANDEM_MAD) 
echo $TANDEM_CORE_MIN
TANDEM_CORE_MAX=$(expr $TANDEM_MEDIAN + 10 \* $TANDEM_MAD)
echo $TANDEM_CORE_MAX

## WIDTH_OF_X_PERCENT
## https://github.com/broadinstitute/picard/blob/959411f3a97e979c406cc068ca9200ff4e2bf6bf/src/main/java/picard/analysis/directed/InsertSizeMetricsCollector.java
## how many reads is 10% around median
expr $TANDEM_READ_PAIRS \/ 10


## RF ##########################################################
PAIR_ORIENTATION="RF"

## { .bedpe.gz should be READ1 READ2 order, 
## so FLAGS column can only be 81,83,97or99 in column 14, 
## then 145,147,161or163 in column 15 }
cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS} \
    | awk '($14 == 99 && $16 < 0) || ($14 == 97 && $16 < 0) || ($14 == 81 && $16 > 0) || ($14 == 83 && $16 > 0) || ($15 == 145 && $17 > 0) || ($15 == 147 && $17 > 0) || ($15 == 161 && $17 < 0) || ($15 == 163 && $17 < 0) {print}' \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs}
## { should only have FLAGS 81,83,97,99,145,147,161,163 left, unless there's markdup... then might have problems... }

## get InsertSizes for RF reads
echo "Getting InsertSizes metrics for RF reads"
cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs} \
    | awk 'BEGIN {OFS="\t"} {print sqrt($16*$16)}' \
    | sort -n \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_COL16}

## READ_PAIRS
RF_READ_PAIRS=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_COL16} | wc -l)
echo $RF_READ_PAIRS

## MEDIAN_INSERT_SIZE (exclude 0's)
RF_MEDIAN=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_COL16} \
    | awk '{ a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) printf "%.0f\n", (a[x-1]+a[x])/2; else printf "%.0f\n", a[x-1]; }')
echo $RF_MEDIAN

## MEDIAN_ABSOLUTE_DEVIATION ( abs value of deviations from median, median )
RF_MAD=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_COL16} \
    | awk -v RF_MEDIAN=$RF_MEDIAN '{$2=$1-RF_MEDIAN; sub("^-","",$2); print $2}' \
    | sort -n \
    | awk '{ a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) printf "%.0f\n", (a[x-1]+a[x])/2; else printf "%.0f\n", a[x-1]; }')
echo $RF_MAD

## MIN
RF_MIN=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_COL16} \
    | sed -n '1p')
echo $RF_MIN

## MAX
RF_MAX=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_COL16} \
    | sed -n '$p')
echo $RF_MAX

## MEAN_INSERT_SIZE ( core inserts only : +/-10 MADs around median insert size )
## 10 is default (https://github.com/broadinstitute/picard/blob/949d7f954f84e02e833be66e8f65553a0918de65/src/main/java/picard/analysis/CollectInsertSizeMetrics.java)
RF_CORE_MIN=$(expr $RF_MEDIAN - 10 \* $RF_MAD)
echo $RF_CORE_MIN
RF_CORE_MAX=$(expr $RF_MEDIAN + 10 \* $RF_MAD)
echo $RF_CORE_MAX

## WIDTH_OF_X_PERCENT
## https://github.com/broadinstitute/picard/blob/959411f3a97e979c406cc068ca9200ff4e2bf6bf/src/main/java/picard/analysis/directed/InsertSizeMetricsCollector.java
## how many reads is 10% around median
expr $RF_READ_PAIRS \/ 10


## FR #####################################################
PAIR_ORIENTATION="FR"

cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS} \
    | awk '($14 == 99 && $16 < 0) || ($14 == 97 && $16 < 0) || ($14 == 81 && $16 > 0) || ($14 == 83 && $16 > 0) || ($15 == 145 && $17 > 0) || ($15 == 147 && $17 > 0) || ($15 == 161 && $17 < 0) || ($15 == 163 && $17 < 0) {next} {print}' \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs}

## get InsertSizes for FR reads
echo "Getting InsertSizes metrics for FR reads"
cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs} \
    | awk 'BEGIN {OFS="\t"} {print sqrt($16*$16)}' \
    | sort -n \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16}

## READ_PAIRS
FR_READ_PAIRS=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16} | wc -l)
echo $FR_READ_PAIRS

## MEDIAN_INSERT_SIZE (exclude 0's)
FR_MEDIAN=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16} \
    | awk '{ a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) printf "%.0f\n", (a[x-1]+a[x])/2; else printf "%.0f\n", a[x-1]; }')
echo $FR_MEDIAN

## MEDIAN_ABSOLUTE_DEVIATION ( abs value of deviations from median, median )
FR_MAD=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16} \
    | awk -v FR_MEDIAN=$FR_MEDIAN '{$2=$1-FR_MEDIAN; sub("^-","",$2); print $2}' \
    | sort -n \
    | awk '{ a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) printf "%.0f\n", (a[x-1]+a[x])/2; else printf "%.0f\n", a[x-1]; }')
echo $FR_MAD

## MIN
FR_MIN=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16} \
    | sed -n '1p')
echo $FR_MIN

## MAX
FR_MAX=$(cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16} \
    | sed -n '$p')
echo $FR_MAX

## MEAN_INSERT_SIZE ( core inserts only : +/-10 MADs around median insert size )
## 10 is default (https://github.com/broadinstitute/picard/blob/949d7f954f84e02e833be66e8f65553a0918de65/src/main/java/picard/analysis/CollectInsertSizeMetrics.java)
FR_CORE_MIN=$(expr $FR_MEDIAN - 35 \* $FR_MAD)
echo $FR_CORE_MIN
FR_CORE_MAX=$(expr $FR_MEDIAN + 35 \* $FR_MAD)
echo $FR_CORE_MAX

cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16} \
    | awk -v FR_CORE_MIN=$FR_CORE_MIN -v FR_CORE_MAX=$FR_CORE_MAX '($1 >= FR_CORE_MIN && $1 <= FR_CORE_MAX) {print}' \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16}.core

#cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16}.core \
cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16} \
    | awk '{ sum += $1 } END { if (NR > 0) printf "%.6f\n", sum/NR }'


## MEAN_INSERT_SIZE and WIDTH_OF_X_PERCENT are just not matching Picard's...



## WIDTH_OF_X_PERCENT
## https://github.com/broadinstitute/picard/blob/959411f3a97e979c406cc068ca9200ff4e2bf6bf/src/main/java/picard/analysis/directed/InsertSizeMetricsCollector.java
## how many reads is 10% around median
expr $FR_READ_PAIRS \/ 10








time2=$(date +%s)
echo "Job ended at "$(date) 
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"

## EOF
