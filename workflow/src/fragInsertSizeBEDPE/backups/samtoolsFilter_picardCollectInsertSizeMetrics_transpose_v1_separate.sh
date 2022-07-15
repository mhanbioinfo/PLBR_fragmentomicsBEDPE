#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -t 1-00:00:00
#SBATCH -D ./ 
#SBATCH --mem=32G
#SBATCH -J bam_picardInsertSize_sep
#SBATCH -p himem
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -o %j-%x.out
#SBATCH -e %j-%x.err

echo "Job started."

module load samtools/1.14
module load bedtools/2.27.1
module load java/8
module load picard/2.10.9
module load perl
module load R/4.1.0

## --------------------------------------------------------- ##

INPUT_DIR="/cluster/projects/pughlab/projects/fragmentomics_bedpe_input/original_files"
#INPUT_DIR="/cluster/projects/pughlab/projects/fragmentomics_bedpe_input/insert_sizes/try11_FR_RF_TANDEM_separate"

#INPUT_F="TGL48_0011_Ct_T_PE_325_WG.merged.sorted.bam"
INPUT_F="TGL48_0011_Ct_T_PE_325_WG.merged.sorted.sub1pct.coordSortd.bam"
#OUT_DIR="/cluster/projects/pughlab/projects/fragmentomics_bedpe_input/insert_sizes/try17_picardInsertSize_fullBAM_separate"
OUT_DIR="/cluster/projects/pughlab/projects/fragmentomics_bedpe_input/insert_sizes/try14_picardInsertSize_sub1pctBAM_alltogether"

## For filtering unmapped/singletons (4, 8), secondary (256), supplementary (2048) reads
FLAG1="F2316"
BAM_F2316_F="${INPUT_F%.*}_${FLAG1}.bam"

## For filtering TANDEM read pairs
FLAG2="f48"
FLAG3="F48"
FLAG16n32_READS="${BAM_F2316_F%.*}_16n32reads.txt"
FLAGno16or32_READS="${BAM_F2316_F%.*}_no16or32reads.txt"
TANDEM_READS="${BAM_F2316_F%.*}_tandemReads.txt"

BAM_FILTD_NO_TANDEMS="${BAM_F2316_F%.*}_noTandems.bam"
BAM_FILTD_ONLY_TANDEMS="${BAM_F2316_F%.*}_onlyTandems.bam"
BAM_FILTD_NO_TANDEMS_NO_TLENeq0="${BAM_FILTD_NO_TANDEMS%.*}_noTLENeq0.bam"

## For filtering RF read pairs, and TLEN=0 read pairs
RF_READS="${BAM_FILTD_NO_TANDEMS_NO_TLENeq0%.*}_RFreads.txt"
BAM_FILTD_NO_TANDEMS_NO_TLENeq0_noRF="${BAM_FILTD_NO_TANDEMS_NO_TLENeq0%.*}_noRF.bam"
BAM_FILTD_NO_TANDEMS_NO_TLENeq0_onlyRF="${BAM_FILTD_NO_TANDEMS_NO_TLENeq0%.*}_onlyRF.bam"

## --------------------------------------------------------- ##
## filtering

## samtools filter1
echo "Peforming samtools filter1... "
samtools view -${FLAG1} \
    -b ${INPUT_DIR}/${INPUT_F} \
    -o ${OUT_DIR}/${BAM_F2316_F} 

## samtools filter2 - both 16 and 32 FLAGS set
echo "Performing samtools filter2... "
samtools view -${FLAG2} \
    ${OUT_DIR}/${BAM_F2316_F} | \
    cut -f1 | sort | uniq \
    > ${OUT_DIR}/${FLAG16n32_READS}

## samtools filter2 - both 16 and 32 FLAGS UNset
echo "Performing samtools filter3... "
samtools view -${FLAG3} \
    ${OUT_DIR}/${BAM_F2316_F} | \
    cut -f1 | sort | uniq \
    > ${OUT_DIR}/${FLAGno16or32_READS}

## cat 2 read name lists
echo "Combine 2 filter lists... "
cat ${OUT_DIR}/${FLAG16n32_READS} ${OUT_DIR}/${FLAGno16or32_READS} \
    > ${OUT_DIR}/${TANDEM_READS}

## picard exclude TANDEM reads
echo "Picard FilterSamReads exclude TANDEM reads... "
java -jar /cluster/tools/software/picard/2.10.9/picard.jar \
    FilterSamReads \
    I=${OUT_DIR}/${BAM_F2316_F} \
    O=${OUT_DIR}/${BAM_FILTD_NO_TANDEMS} \
    READ_LIST_FILE=${OUT_DIR}/${TANDEM_READS} \
    FILTER=excludeReadList \
    WRITE_READS_FILES=false \
    USE_JDK_DEFLATER=true \
    USE_JDK_INFLATER=true

## picard extract TANDEM reads
echo "Picard FilterSamReads extract TANDEM reads... "
java -jar /cluster/tools/software/picard/2.10.9/picard.jar \
    FilterSamReads \
    I=${OUT_DIR}/${BAM_F2316_F} \
    O=${OUT_DIR}/${BAM_FILTD_ONLY_TANDEMS} \
    READ_LIST_FILE=${OUT_DIR}/${TANDEM_READS} \
    FILTER=includeReadList \
    WRITE_READS_FILES=false \
    USE_JDK_DEFLATER=true \
    USE_JDK_INFLATER=true

# remove TLEN eq 0 reads
echo "removing TLEN eq 0 reads... "
samtools view -H ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS} \
    > ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS%.*}_header.sam
samtools view ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS} | \
    awk '$9 != 0' | \
    cat ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS%.*}_header.sam - | \
    samtools view -b - > ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0}

## remove RF outieReads
echo "removing RF outie reads... "
samtools view ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0} | \
    awk '($2 == 97 && $9 < 0) || ($2 == 99 && $9 < 0) || ($2 == 145 && $9 > 0) || ($2 == 147 && $9 > 0) || ($2 == 161 && $9 < 0) || ($2 == 163 && $9 < 0) || ($2 == 81 && $9 > 0) || ($2 == 83 && $9 > 0)' | \
    cut -f1 \
    > ${OUT_DIR}/${RF_READS}

## picard exclude RF reads
echo "Picard FilterSamReads exclude RF reads... "
java -jar /cluster/tools/software/picard/2.10.9/picard.jar \
    FilterSamReads \
    I=${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0} \
    O=${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_noRF} \
    READ_LIST_FILE=${OUT_DIR}/${RF_READS} \
    FILTER=excludeReadList \
    WRITE_READS_FILES=false \
    USE_JDK_DEFLATER=true \
    USE_JDK_INFLATER=true

## picard extract RF reads
echo "Picard FilterSamReads extract RF reads... "
java -jar /cluster/tools/software/picard/2.10.9/picard.jar \
    FilterSamReads \
    I=${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0} \
    O=${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_onlyRF} \
    READ_LIST_FILE=${OUT_DIR}/${RF_READS} \
    FILTER=includeReadList \
    WRITE_READS_FILES=false \
    USE_JDK_DEFLATER=true \
    USE_JDK_INFLATER=true

## --------------------------------------------------------- ##
## get stats

## samtools flagstat
echo "Performing samtools flagstat... "
samtools flagstat \
    ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_noRF} \
    > ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_noRF%.*}.flagstat
samtools flagstat \
    ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_onlyRF} \
    > ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_onlyRF%.*}.flagstat
samtools flagstat \
    ${OUT_DIR}/${BAM_FILTD_ONLY_TANDEMS} \
    > ${OUT_DIR}/${BAM_FILTD_ONLY_TANDEMS%.*}.flagstat

## count flags
echo "counting flags... "
samtools view ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_noRF} | \
    cut -f2 | sort | uniq -c \
    > ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_noRF}_flags_count.txt
samtools view ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_onlyRF} | \
    cut -f2 | sort | uniq -c \
    > ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_onlyRF}_flags_count.txt
samtools view ${OUT_DIR}/${BAM_FILTD_ONLY_TANDEMS} | \
    cut -f2 | sort | uniq -c \
    > ${OUT_DIR}/${BAM_FILTD_ONLY_TANDEMS}_flags_count.txt

## count TLEN
echo "counting TLEN... "
samtools view ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_noRF} | \
    cut -f9 | grep -v '-' | sort | uniq -c | sort -k1V \
    > ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_noRF}_TLENcount.txt
samtools view ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_onlyRF} | \
    cut -f9 | grep -v '-' | sort | uniq -c | sort -k1V \
    > ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_onlyRF}_TLENcount.txt
samtools view ${OUT_DIR}/${BAM_FILTD_ONLY_TANDEMS} | \
    cut -f9 | grep -v '-' | sort | uniq -c | sort -k1V \
    > ${OUT_DIR}/${BAM_FILTD_ONLY_TANDEMS}_TLENcount.txt

## picard CollectInsertSizeMetrics
echo "Performing picard CollectInsertSizeMetrics... "
java -jar /cluster/tools/software/picard/2.10.9/picard.jar \
    CollectInsertSizeMetrics \
    INCLUDE_DUPLICATES=true \
    I=${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_noRF} \
    O=${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_noRF%.*}_picardInsertSize_wDups.txt \
    H=${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_noRF%.*}_picardInserSizeHisto_wDups.pdf \
    M=0 W=600 
java -jar /cluster/tools/software/picard/2.10.9/picard.jar \
    CollectInsertSizeMetrics \
    INCLUDE_DUPLICATES=true \
    I=${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_onlyRF} \
    O=${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_onlyRF%.*}_picardInsertSize_wDups.txt \
    H=${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_onlyRF%.*}_picardInserSizeHisto_wDups.pdf \
    M=0 W=600
java -jar /cluster/tools/software/picard/2.10.9/picard.jar \
    CollectInsertSizeMetrics \
    INCLUDE_DUPLICATES=true \
    I=${OUT_DIR}/${BAM_FILTD_ONLY_TANDEMS} \
    O=${OUT_DIR}/${BAM_FILTD_ONLY_TANDEMS%.*}_picardInsertSize_wDups.txt \
    H=${OUT_DIR}/${BAM_FILTD_ONLY_TANDEMS%.*}_picardInserSizeHisto_wDups.pdf \
    M=0 W=600

## transpose picard output
echo "Performing transpose on picard output... "
cat ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_noRF%.*}_picardInsertSize_wDups.txt | \
    grep -v "^#" | head -5 | \
    perl ./transpose.pl - > ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_noRF%.*}_picardInsertSize_wDups_transposed.txt
cat ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_onlyRF%.*}_picardInsertSize_wDups.txt | \
    grep -v "^#" | head -5 | \
    perl ./transpose.pl - > ${OUT_DIR}/${BAM_FILTD_NO_TANDEMS_NO_TLENeq0_onlyRF%.*}_picardInsertSize_wDups_transposed.txt
cat ${OUT_DIR}/${BAM_FILTD_ONLY_TANDEMS%.*}_picardInsertSize_wDups.txt | \
    grep -v "^#" | head -5 | \
    perl ./transpose.pl - > ${OUT_DIR}/${BAM_FILTD_ONLY_TANDEMS%.*}_picardInsertSize_wDups_transposed.txt

echo "Job ended."
## EOF
