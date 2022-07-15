#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -t 1-00:00:00
#SBATCH -D ./
#SBATCH --mem=32G
#SBATCH -J get_clean_bedpe4fragmentomics
#SBATCH -p himem
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -o %j-%x.out
#SBATCH -e %j-%x.err

#BEDPE_PATH="/Users/minghan/bioinfoproj/fragmentomics/pipeline/frag_ratio/data/TGL48_0001_Ct_T_PE_312_WG.merged.sorted_coordSortd.bedpe.gz"
#OUT_PATH="/Users/minghan/bioinfoproj/fragmentomics/pipeline/frag_ratio/fragmentomicsBEDPE_testing/TGL48_0001_Ct_T_PE_312_WG.merged.sorted_coordSortd.df4fragm"


# getopts ###################################################
usage(){
    echo 
    echo "Usage: bash get_clean_bedpe4fragmentomics.sh -i BEDPE_PATH -o OUT_PATH"
    echo 
}
no_args="true"

## Help 
Help()
{
    # Display Help
    echo 
    echo "Process .bedpe.gz into a lean .bedpe for input into Fragmentomics pipeline."
    echo
    echo "Usage: bash get_clean_bedpe4fragmentomics.sh -b BEDPE_PATH -o OUT_PATH"
    echo "options:"
    echo "-h   [HELP]      print help"
    echo "-i   [REQUIRED]  input .bedpe.gz (full path)"
    echo "-o   [REQUIRED]  output file (full path)"
    echo
}

## Get the options
while getopts ":hi:o:" option; do
    case "${option}" in
        h) Help
           exit;;
        i) BEDPE_PATH=${OPTARG};;
        o) OUT_PATH=${OPTARG};;
       \?) echo "Error: Invalid option"
           exit;;
    esac
    no_args="false"
done

[[ "$no_args" == "true" ]] && { usage; exit 1; }

echo "input .bedpe.gz path:    $BEDPE_PATH"
echo "output .df4fragm path:   $OUT_PATH"


# Main program ##############################################

echo "Job started at "$(date) 
time1=$(date +%s)

gzcat "${BEDPE_PATH}" \
    | awk '(!and($14,0x4)) {print}'                               `## isUnmappedQuery = FALSE` \
    | awk '(!and($14,0x8)) {print}' \
    | awk '(!and($15,0x4)) {print}' \
    | awk '(!and($15,0x8)) {print}' \
    | awk '($8>=30 && $9>=30) {print}'                            `## mapQ >= 30` \
    | awk '(and($14,0x2)) {print}'                                `## isProperPair = TRUE` \
    | awk 'BEGIN{OFS="\t"} (!and($14,0x400) && 
                            !and($15,0x400)) {print}'             `## isDuplicate = FALSE` \
    | awk '(!and($14,0x100)) {print}'                             `## isSecondaryAlignment = FALSE` \
    | awk '(!and($15,0x100)) {print}' \
    | awk '$1==$4 {print}'                                        `## remove ambiguous and supplementary reads` \
    | awk '((and($14,0x800) || and($15,0x800)) && $1==$4) {next} {print}' \
    | awk '(($14 == 113 && $15 == 177) || ($14 == 177 && $15 == 113)) {next} {print}' \
    | awk '(and($14,0x40) && and($15,0x80)) ||
           (and($14,0x80) && and($15,0x40)) {print}' \
    | awk '(!and($14,0x100) && !and($15, 0x100)) ||
           (and($14,0x100) && and($15, 0x100)) {print}' \
    | awk '(and($14,0x10) && and($15,0x20)) ||
           (and($14,0x20) && and($15,0x10)) {print}' \
    | awk 'BEGIN{OFS="\t"} {print $1, $10, $2+=1, $3, $5+=1, $6}' `## get only columns needed for GAlignmentPairs R object` \
    | awk 'BEGIN{OFS="\t"} ($2 == "+") {print $1,$3,$6,$2}
                           ($2 == "-") {print $1,$5,$4,$2}'       `## swap ranges depending on strand, get frag start-end` \
    > "${OUT_PATH}"

cat "${OUT_PATH}" | wc -l
# 24369546


time2=$(date +%s)
echo "Job ended at "$(date) 
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"

## EOF