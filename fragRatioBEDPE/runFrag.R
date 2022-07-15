# file: runFrag.R
# author: Derek Wong, Ph.D
# date: October 5th, 2021
# edited by: Ming Han
# date: July 15th, 2022

library(optparse)

## Set script variables
option_list <- list(
  make_option(c("--id"), type = "character", help = "sample id. Required"),
  make_option(c("--input_type"), type = "character", help = "bam or bedpe. Required"),
  make_option(c("--input_path"), type = "character", help = "Path to input file. Required."),
  make_option(c("--filters"), type = "character", help = "Path to genomic blacklist regions. Required."),
  make_option(c("--gaps"), type = "character", help = "Path to genome gaps. Required."),
  make_option(c("--tiles"), type = "character", help = "Path to 100kb tiled genome. Required."),
  make_option(c("--VNTRs"), type = "character", help = "Path to VNTRs. Required."),
  make_option(c("--healthy"), type = "character", help = "Path to panel of healthy controls. Required."),
  make_option(c("--outdir"), type = "character", help = "Path to output directory. Required."),
  make_option(c("--libdir"), type = "character", help = "Path to scripts. Required.")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

## Load required packages
library(tidyverse)
# library(multidplyr)
library(Rsamtools)
library(rtracklayer)
library(GenomicAlignments)
library(GenomicRanges)
library(devtools)
# library(Homo.sapiens)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg38)
library(readxl)
# class(Homo.sapiens)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')

## Get variables from input script
id <- opt$id

if (opt$input_type == "bam"){
  bam <- opt$input_path
} else if (opt$input_type == "bedpe"){
  bedpe <- opt$input_path
}

filters <- opt$filters
gaps <- opt$gaps
tiles <- opt$tiles
VNTRs <- opt$VNTRs
healthy <- opt$healthy
libdir <- opt$libdir
outdir <- file.path(opt$outdir, id)


# ## manual inputs for testing
# frag_main = "/Users/minghan/GDrive_minghanpughlab/PughLabPMH/_projects/fragmentomics/analysis003_pipeline/fragmentomicsBEDPE/"
# id <- "TGL48_0001_Ct_T_PE_312_WG.merged.sorted"
# bam <- "/Users/minghan/bioinfoproj/fragmentomics/pipeline/frag_ratio/data/TGL48_0001_Ct_T_PE_312_WG.merged.sorted.bam"
# bedpe = "/Users/minghan/bioinfoproj/fragmentomics/pipeline/frag_ratio/fragmentomicsBEDPE_testing/TGL48_0001_Ct_T_PE_312_WG.merged.sorted_coordSortd.df4fragm"
# filters <- paste0(frag_main, "/extdata/filters.hg38.rda")
# gaps <- paste0(frag_main, "/extdata/gaps.hg38.rda")
# tiles <- paste0(frag_main, "/extdata/hg38_tiles.bed")
# VNTRs <- paste0(frag_main, "/extdata/VNTRs.hg38.rda")
# healthy <- paste0(frag_main, "/extdata/healthy.median.hg38.rda")
# libdir <- frag_main
# outdir <- "/Users/minghan/bioinfoproj/fragmentomics/pipeline/frag_ratio/fragmentomicsBEDPE_output/"


## Create output directory
dir.create(outdir)

## Run scripts
if (opt$input_type == "bam"){
  source(paste0(libdir,"/R/git_01-read_fragmentsBAM.R"))
} else if (opt$input_type == "bedpe"){
  source(paste0(libdir,"/R/git_01-read_fragmentsBEDPE.R"))
}
source(paste0(libdir,"/R/git_02-mito_frag.R"))
source(paste0(libdir,"/R/git_03-100kb_bins.R"))
source(paste0(libdir,"/R/git_03-plot_GC_correction.R"))
source(paste0(libdir,"/R/git_04-5Mb_bins.R"))
source(paste0(libdir,"/R/git_05-summary.R"))
source(paste0(libdir,"/R/git_06-plotting.R"))
#source(paste0(libdir,"/R/git_07-gbm_prediction.R"))

q('no')

## EOF