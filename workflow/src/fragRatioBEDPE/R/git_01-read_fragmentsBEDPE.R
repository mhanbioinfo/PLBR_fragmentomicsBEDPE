# file: git_01-read_bam.R
# author: Derek Wong, Ph.D
# date: October 5th, 2021
# edited by: Ming Han
# date: July 15th, 2022

## Read in lean bedpe for fragmentomics
bedpe4fragm =
  read.table(file = bedpe, header = F, sep = '\t')
bedpe4fragm %>% dim() # [1] 24369546        4
bedpe4fragm %>% head()
names(bedpe4fragm) = c("chr", "start", "end", "strand")

bedpe4fragm.chr1to22 =
  bedpe4fragm[bedpe4fragm$chr %in% paste0("chr", 1:22),]
bedpe4fragm.chr1to22 %>% dim() # [1] 23584457        4

bedpe4fragm.mito =
  bedpe4fragm[bedpe4fragm$chr == "chrM",]
bedpe4fragm.mito %>% dim() # [1] 150   4
bedpe4fragm.mito %>% head()


## Filter reads: 90-220bp on Autosomes and mitochondrial reads
frags = makeGRangesFromDataFrame(df = bedpe4fragm.chr1to22, keep.extra.columns = F)
frags
mito = makeGRangesFromDataFrame(df = bedpe4fragm.mito, keep.extra.columns = F)
mito

length <- nrow(bedpe4fragm)
rm(bedpe4fragm)

w.all <- width(frags)
frags <- frags[which(w.all >= 90 & w.all <= 220)]
rm(w.all)

gcs <- GCcontent(Hsapiens, unstrand(frags))
frags$gc <- gcs
rm(gcs)

# frags
## EOF