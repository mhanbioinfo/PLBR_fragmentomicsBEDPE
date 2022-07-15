#!/bin/bash

PROJ_DIR="/Users/minghan/bioinfoproj/fragmentomics/pipeline/frag_ratio"
INPUT_ID="TGL48_0001_Ct_T_PE_312_WG.merged.sorted"
#INPUT_TYPE="bam"
#INPUT_PATH="/Users/minghan/bioinfoproj/fragmentomics/pipeline/frag_ratio/data/TGL48_0001_Ct_T_PE_312_WG.merged.sorted.bam"
INPUT_TYPE="bedpe"
INPUT_PATH="/Users/minghan/bioinfoproj/fragmentomics/pipeline/frag_ratio/data/TGL48_0001_Ct_T_PE_312_WG.merged.sorted_coordSortd.bedpe.gz"

FRAG_PL="/Users/minghan/GDrive_minghanpughlab/PughLabPMH/_projects/fragmentomics/analysis003_pipeline/fragmentomicsBEDPE/"

if [ $INPUT_TYPE == "bedpe" ]; then
    BEDPE4FRAGM="${INPUT_PATH%.bedpe.gz}.df4fragm"
    bash $FRAG_PL/get_clean_bedpe4fragmentomics.sh \
        -i $INPUT_PATH \
        -o $BEDPE4FRAGM
    INPUT_PATH="${BEDPE4FRAGM}"
    echo $INPUT_PATH
fi


Rscript $FRAG_PL/runFrag.R \
    --id $INPUT_ID \
    --input_type $INPUT_TYPE \
    --input_path $INPUT_PATH \
    --filters $FRAG_PL/extdata/filters.hg38.rda \
    --gaps $FRAG_PL/extdata/gaps.hg38.rda \
    --VNTRs $FRAG_PL/extdata/VNTRs.hg38.rda \
    --tiles $FRAG_PL/extdata/hg38_tiles.bed \
    --healthy $FRAG_PL/extdata/healthy.median.hg38.rda \
    --outdir $PROJ_DIR \
    --libdir $FRAG_PL


## EOF