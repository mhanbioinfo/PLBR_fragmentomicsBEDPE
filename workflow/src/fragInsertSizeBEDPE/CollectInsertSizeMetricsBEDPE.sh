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

module load R/4.1.0

echo "Job started at "$(date) 
time1=$(date +%s)

SRC_DIR=$(pwd)
RSCRIPT_INSERTSIZEMETRICS="${SRC_DIR}/column16_picardInsertSizeMetrics.R"
RSCRIPT_INSERTSIZE_HISTOGRAM="${SRC_DIR}/insertSizeHistogramBEDPE.R"

INPUT_DIR="/cluster/projects/pughlab/projects/fragmentomics_bedpe_input/original_files"
INPUT_BEDPE="TGL48_0011_Ct_T_PE_325_WG.merged.sorted_coordSortd.bedpe.gz"

OUT_DIR="/cluster/projects/pughlab/projects/fragmentomics_bedpe_input/insert_sizes/try23_CollectInsertSizeMetricsBEDPE_v2_fullBAM_histogram"

HISTOGRAM_WIDTH=600

OUT_PATH_INSERTSIZE_METRIC="${OUT_DIR}/${INPUT_BEDPE%.bedpe.gz}.insertSizeMetric.parsed"
OUT_PATH_INSERTSIZE_DF="${OUT_DIR}/${INPUT_BEDPE%.bedpe.gz}.insertSizeMetric_DF"
OUT_PATH_INSERTSIZE_HISTOGRAM="${OUT_DIR}/${INPUT_BEDPE%.bedpe.gz}.insertSizeMetric_histogram.pdf"

## intermediate filenames ###############################################

BEDPE_F2316_TLENgt0="${INPUT_BEDPE%.bedpe.gz}.F2316_TLENgt0.bedpe"
BEDPE_F2316_TLENgt0_TANDEMS="${BEDPE_F2316_TLENgt0%.*}.TANDEMS.bedpe"
BEDPE_F2316_TLENgt0_NOTANDEMS="${BEDPE_F2316_TLENgt0%.*}.NOTANDEMS.bedpe"
BEDPE_F2316_TLENgt0_TANDEMS_COL16="${BEDPE_F2316_TLENgt0_TANDEMS%.*}.col16"
BEDPE_F2316_TLENgt0_TANDEMS_INSERTSIZE_DF="${BEDPE_F2316_TLENgt0_TANDEMS%.*}.insertSize_DF"

BEDPE_F2316_TLENgt0_NOTANDEMS_RFs="${BEDPE_F2316_TLENgt0_NOTANDEMS%.*}.RFs.bedpe"
BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_COL16="${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs%.*}.col16"
BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_INSERTSIZE_DF="${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs%.*}.insertSize_DF"

BEDPE_F2316_TLENgt0_NOTANDEMS_FRs="${BEDPE_F2316_TLENgt0_NOTANDEMS%.*}.FRs.bedpe"
BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16="${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs%.*}.col16"
BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_INSERTSIZE_DF="${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs%.*}.insertSize_DF"

## script ###############################################################

echo "Removing unmapped/singletons (4, 8), secondary (256), supplementary (2048) reads, or if TLEN == 0"
zcat ${INPUT_DIR}/${INPUT_BEDPE} \
    | awk 'BEGIN {OFS="\t"} ($16 != 0 && $17 != 0) {print}' \
    | awk 'BEGIN {OFS="\t"} ((!and($14, 0x4)) && (!and($15, 0x4))) {print}' \
    | awk 'BEGIN {OFS="\t"} ((!and($14, 0x8)) && (!and($15, 0x8))) {print}' \
    | awk 'BEGIN {OFS="\t"} ((!and($14, 0x100)) && (!and($15, 0x100))) {print}' \
    | awk 'BEGIN {OFS="\t"} ((!and($14, 0x800)) && (!and($15, 0x800))) {print}' \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0}

echo "Extracting tandem reads (16, 32 both set or unset)"
cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0} \
    | awk 'BEGIN {OFS="\t"} ((!and($14, 0x10)) && (!and($14, 0x20))) || ((!and($15, 0x10)) && (!and($15, 0x20))) || ((and($14, 0x10)) && (and($14, 0x20))) || ((and($15, 0x10)) && (and($15, 0x20))) {print}' \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS}

echo "Removing tandem reads (16, 32 both set or unset)"
cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0} \
    | awk 'BEGIN {OFS="\t"} ((!and($14, 0x10)) && (!and($14, 0x20))) {next} {print}' \
    | awk 'BEGIN {OFS="\t"} ((!and($15, 0x10)) && (!and($15, 0x20))) {next} {print}' \
    | awk 'BEGIN {OFS="\t"} ((and($14, 0x10)) && (and($14, 0x20))) {next} {print}' \
    | awk 'BEGIN {OFS="\t"} ((and($15, 0x10)) && (and($15, 0x20))) {next} {print}' \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS}


## TANDEM READS ####################################################
PAIR_ORIENTATION="TANDEM"

echo "Getting InsertSize column for TANDEM reads"
cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS} \
    | awk 'BEGIN {OFS="\t"} {print sqrt($16*$16)}' \
    | sort -n \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS_COL16}
#    | awk 'function SUB(F) {sub("^-","",$F)} SUB(16) SUB(17); {print $16, $17}' \
#    | awk '$1 != $2 {print}' \
## TLEN columns 16 and 17 always same number, just different signs

echo "Getting InsertSize METRICS for TANDEM reads"
Rscript $RSCRIPT_INSERTSIZEMETRICS \
    --input "${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS_COL16}" \
    --orientation "TANDEM" \
    --output "${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS_COL16}.insertSizeMetrics" \
    --insertsizedf "${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS_INSERTSIZE_DF}" \
    --histogram_width ${HISTOGRAM_WIDTH}

## RF ##########################################################
PAIR_ORIENTATION="RF"

## { .bedpe.gz should be READ1 READ2 order, 
## so FLAGS column can only be 81,83,97or99 in column 14, 
## then 145,147,161or163 in column 15 }
cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS} \
    | awk '($14 == 99 && $16 < 0) || ($14 == 97 && $16 < 0) || ($14 == 81 && $16 > 0) || ($14 == 83 && $16 > 0) || ($15 == 145 && $17 > 0) || ($15 == 147 && $17 > 0) || ($15 == 161 && $17 < 0) || ($15 == 163 && $17 < 0) {print}' \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs}
## { should only have FLAGS 81,83,97,99,145,147,161,163 left, unless there's markdup... then might have problems... might need to have "includeDuplicates" option }

echo "Getting InsertSize column for RF reads"
cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs} \
    | awk 'BEGIN {OFS="\t"} {print sqrt($16*$16)}' \
    | sort -n \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_COL16}

echo "Getting InsertSize METRICS for RF reads"
Rscript $RSCRIPT_INSERTSIZEMETRICS \
    --input "${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_COL16}" \
    --orientation "RF" \
    --output "${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_COL16}.insertSizeMetrics" \
    --insertsizedf "${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_INSERTSIZE_DF}" \
    --histogram_width ${HISTOGRAM_WIDTH}

## FR #####################################################
PAIR_ORIENTATION="FR"

cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS} \
    | awk '($14 == 99 && $16 < 0) || ($14 == 97 && $16 < 0) || ($14 == 81 && $16 > 0) || ($14 == 83 && $16 > 0) || ($15 == 145 && $17 > 0) || ($15 == 147 && $17 > 0) || ($15 == 161 && $17 < 0) || ($15 == 163 && $17 < 0) {next} {print}' \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs}

echo "Getting InsertSize column for FR reads"
cat ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs} \
    | awk 'BEGIN {OFS="\t"} {print sqrt($16*$16)}' \
    | sort -n \
    > ${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16}

echo "Getting InsertSize METRICS for FR reads"
Rscript $RSCRIPT_INSERTSIZEMETRICS \
    --input "${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16}" \
    --orientation "FR" \
    --output "${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16}.insertSizeMetrics" \
    --insertsizedf "${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_INSERTSIZE_DF}" \
    --histogram_width ${HISTOGRAM_WIDTH}

echo "Collating InsertSize METRICS for FR, RF, TANDEM"
paste "${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_COL16}.insertSizeMetrics" \
      <(cut -f2 "${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_COL16}.insertSizeMetrics") \
      <(cut -f2 "${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS_COL16}.insertSizeMetrics") \
      > ${OUT_PATH_INSERTSIZE_METRIC}

echo "Join InsertSize DF for FR, RF, TANDEM"
join -t $'\t' -a 1 -a 2 -o0,1.2,2.2 -e '0' \
    <(sort -k1,1 "${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_FRs_INSERTSIZE_DF}") \
    <(sort -k1,1 "${OUT_DIR}/${BEDPE_F2316_TLENgt0_NOTANDEMS_RFs_INSERTSIZE_DF}") \
    | join --nocheck-order -t $'\t' -a 1 -a 2 -o0,1.2,1.3,2.2 -e '0' \
    <(sort -k1,1 -) <(sort -k1,1 "${OUT_DIR}/${BEDPE_F2316_TLENgt0_TANDEMS_INSERTSIZE_DF}") \
    | sort -n \
    | awk -v HISTOGRAM_WIDTH=${HISTOGRAM_WIDTH} '$1 ~/insert_size/ || $1 <= HISTOGRAM_WIDTH {print}' \
    > ${OUT_PATH_INSERTSIZE_DF}
## { picard's output has a blank line at the end... }
echo "" >> ${OUT_PATH_INSERTSIZE_DF}

echo "Generating InsertSize Histogram"
Rscript $RSCRIPT_INSERTSIZE_HISTOGRAM \
    ${OUT_PATH_INSERTSIZE_DF} \
    ${OUT_PATH_INSERTSIZE_HISTOGRAM} \
    ${INPUT_BEDPE%.bedpe.gz}
 

time2=$(date +%s)
echo "Job ended at "$(date) 
echo "Job took $(((time2-time1)/3${HISTOGRAM_WIDTH})) hours $((((time2-time1)%3${HISTOGRAM_WIDTH})/60)) minutes $(((time2-time1)%60)) seconds"


## resources
## CollectInsertSizeMetrics (Picard) usage:
## https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard-
## CollectInsertSizeMetrics (Picard) output definitions:
## https://broadinstitute.github.io/picard/picard-metric-definitions.html#InsertSizeMetrics
## CollectInsertSizeMetrics (Picard) github code:
## https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/analysis/directed/InsertSizeMetricsCollector.java
## https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/analysis/CollectInsertSizeMetrics.java


## EOF
