library(docopt)
doc = 'column16_picardInsertSizeMetrics_v3.R

Usage:
  column16_picardInsertSizeMetrics_v3.R --input INPUT --orientation ORIENTATION --output OUTPUT [ --histogram HISTOGRAM ] [ --histogram_width HISTOGRAM_WIDTH ]

Options:
  --input INPUT             Input file path, column 16 of bedpe file (FR, RF, TANDEM orientation separated).
  --orientation ORIENTATION Read 1 and 2 orientation, FR, RF or TANDEM.
  --output OUTPUT           Output file path, wrangled picardInsertSizeMetrics output.
  --histogram HISTOGRAM     File to write insert size Histogram chart to.
  --histogram_width HISTOGRAM_WIDTH   Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. Also, when calculating mean and standard deviation, only bins <= Histogram_WIDTH will be included.
'

args = docopt(doc, version='picardInsertSizeMetricsBEDPE v0.2.0')

# args = list()
# args$histogram_width = 600
# args$input = "/Users/minghan/bioinfoproj/fragmentomics/column16/TGL48_0011_Ct_T_PE_325_WG.merged.sorted.sub1pct.coordSortd_coordSortd.F2316_TLENgt0.NOTANDEMS.FRs.col16"
# args$input = "/Users/minghan/bioinfoproj/fragmentomics/column16/TGL48_0011_Ct_T_PE_325_WG.merged.sorted.sub1pct.coordSortd_coordSortd.F2316_TLENgt0.NOTANDEMS.RFs.col16"
# args$input = "/Users/minghan/bioinfoproj/fragmentomics/column16/TGL48_0011_Ct_T_PE_325_WG.merged.sorted.sub1pct.coordSortd_coordSortd.F2316_TLENgt0.TANDEMS.col16"
# args$orientation = "FR"
# args$orientation = "RF"
# args$orientation = "TANDEM"
# args$output = "/Users/minghan/bioinfoproj/fragmentomics/column16/TGL48_0011_Ct_T_PE_325_WG.merged.sorted.sub1pct.coordSortd_coordSortd.F2316_TLENgt0.NOTANDEMS.FRs.col16.insertSizeMetrics"
# args$output = "/Users/minghan/bioinfoproj/fragmentomics/column16/TGL48_0011_Ct_T_PE_325_WG.merged.sorted.sub1pct.coordSortd_coordSortd.F2316_TLENgt0.NOTANDEMS.RFs.col16.insertSizeMetrics"
# args$output = "/Users/minghan/bioinfoproj/fragmentomics/column16/TGL48_0011_Ct_T_PE_325_WG.merged.sorted.sub1pct.coordSortd_coordSortd.F2316_TLENgt0.TANDEMS.col16.insertSizeMetrics"

## picard CollectInsertSizeMetrics parameters
DEVIATION = 10
MINIMUM_PCT = 0

## read in data
DATA = read.table(file = args$input)
INSERT_VEC = DATA$V1

## initial InsertSizeMetrics
MEDIAN = median(INSERT_VEC)
MAD = mad(INSERT_VEC, center = MEDIAN, constant = 1)
MIN = min(INSERT_VEC)
MAX = max(INSERT_VEC)
# MEDIAN; MAD; MIN; MAX

if (is.null(args$histogram_width)){
  ## CORE is default
  CORE_RANGE = c(MEDIAN - DEVIATION*MAD, MEDIAN + DEVIATION*MAD)
  MEAN = mean(INSERT_VEC[INSERT_VEC > CORE_RANGE[1] & INSERT_VEC < CORE_RANGE[2]])
  STDV = sd(INSERT_VEC[INSERT_VEC > CORE_RANGE[1] & INSERT_VEC < CORE_RANGE[2]])
} else {
  ## if HISTOGRAM_WIDTH specified, then mean and stdv calculated according to HISTOGRAM_WIDTH
  HISTOGRAM_WIDTH = as.integer(args$histogram_width)
  MEAN = mean(INSERT_VEC[INSERT_VEC <= HISTOGRAM_WIDTH])
  STDV = sd(INSERT_VEC[INSERT_VEC <= HISTOGRAM_WIDTH])
}
if (is.na(MEAN)){
  MEAN = 0
}
if (is.na(STDV)){
  STDV = 0
}
# MEAN; STDV

READ_PAIRS = length(INSERT_VEC)

## WIDTH_OF_X_PCTG ( uses ALL reads )
INSERT_VEC_LEN = length(INSERT_VEC)
WIDTH_OF_X_PCTG_LIST = list()
WIDTH_OF_X_PCTG_LIST["WIDTH_OF_10_PERCENT"] = 0
WIDTH_OF_X_PCTG_LIST["WIDTH_OF_20_PERCENT"] = 0
WIDTH_OF_X_PCTG_LIST["WIDTH_OF_30_PERCENT"] = 0
WIDTH_OF_X_PCTG_LIST["WIDTH_OF_40_PERCENT"] = 0
WIDTH_OF_X_PCTG_LIST["WIDTH_OF_50_PERCENT"] = 0
WIDTH_OF_X_PCTG_LIST["WIDTH_OF_60_PERCENT"] = 0
WIDTH_OF_X_PCTG_LIST["WIDTH_OF_70_PERCENT"] = 0
WIDTH_OF_X_PCTG_LIST["WIDTH_OF_80_PERCENT"] = 0
WIDTH_OF_X_PCTG_LIST["WIDTH_OF_90_PERCENT"] = 0
# WIDTH_OF_X_PCTG_LIST["WIDTH_OF_95_PERCENT"] = 0
WIDTH_OF_X_PCTG_LIST["WIDTH_OF_99_PERCENT"] = 0

INSERT_VEC.TABLE = as.data.frame(table(INSERT_VEC), stringsAsFactors = F)
INSERT_VEC.TABLE$INSERT_VEC = as.integer(INSERT_VEC.TABLE$INSERT_VEC)
# tail(INSERT_VEC.TABLE, 500)

## Picard way #######################################################
## going from median outwards one BASE at a time until max insert size

COVERED = 0
PCT_COVERED = 0
LOW = MEDIAN
HIGH = MEDIAN
# LOW;HIGH

while ((LOW >= MIN - 1 || HIGH <= MAX + 1) && PCT_COVERED < 0.9901 && LOW > 0){
  ## BIN is freq count
  LOWBIN = INSERT_VEC.TABLE[INSERT_VEC.TABLE$INSERT_VEC == LOW,]$Freq
  # print(paste0("LOWBIN: ", LOWBIN))
  if (length(LOWBIN) > 0){
    COVERED = COVERED + LOWBIN
    # print(paste0("COVERED_LOWBIN: ", COVERED))
  }
  ## so when initially LOW=HIGH=MEDIAN, only count MEDIAN once
  if (LOW != HIGH){
    HIGHBIN = INSERT_VEC.TABLE[INSERT_VEC.TABLE$INSERT_VEC == HIGH,]$Freq
    # print(paste0("HIGHBIN: ", HIGHBIN))
    if (length(HIGHBIN) > 0){
      COVERED = COVERED + HIGHBIN
    }
  }
  # print(paste0("HIGH AND LOW COVERED: ", COVERED))
  PCT_COVERED = COVERED / READ_PAIRS
  # print(paste0("PCT_COVERED: ", PCT_COVERED))
  DISTANCE = HIGH - LOW + 1
  # print(paste0("DISTANCE: ", DISTANCE))

  if (PCT_COVERED >= 0.1 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_10_PERCENT"] == 0){ WIDTH_OF_X_PCTG_LIST["WIDTH_OF_10_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.2 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_20_PERCENT"] == 0){ WIDTH_OF_X_PCTG_LIST["WIDTH_OF_20_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.3 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_30_PERCENT"] == 0){ WIDTH_OF_X_PCTG_LIST["WIDTH_OF_30_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.4 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_40_PERCENT"] == 0){ WIDTH_OF_X_PCTG_LIST["WIDTH_OF_40_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.5 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_50_PERCENT"] == 0){ WIDTH_OF_X_PCTG_LIST["WIDTH_OF_50_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.6 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_60_PERCENT"] == 0){ WIDTH_OF_X_PCTG_LIST["WIDTH_OF_60_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.7 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_70_PERCENT"] == 0){ WIDTH_OF_X_PCTG_LIST["WIDTH_OF_70_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.8 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_80_PERCENT"] == 0){ WIDTH_OF_X_PCTG_LIST["WIDTH_OF_80_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.9 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_90_PERCENT"] == 0){ WIDTH_OF_X_PCTG_LIST["WIDTH_OF_90_PERCENT"] = DISTANCE }
  # if (PCT_COVERED >= 0.95 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_95_PERCENT"] == 0){ WIDTH_OF_X_PCTG_LIST["WIDTH_OF_95_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.99 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_99_PERCENT"] == 0){ WIDTH_OF_X_PCTG_LIST["WIDTH_OF_99_PERCENT"] = DISTANCE }

  LOW = LOW - 1;
  # print(paste0("new LOW: ", LOW))
  HIGH = HIGH + 1;
  # print(paste0("new HIGH: ", HIGH))
}
# WIDTH_OF_X_PCTG_LIST

## when LOW becomes negative, will not have any more LOWBINs
## distance will now just = (insert size - median)*2 + 1
## no need to perform the while loop above for up to MAX insert size, will take forever if MAX insert size is > 1e6
# LOW;HIGH
HIGH_ONLY = HIGH
INSERT_VEC.TABLE.HIGH_ONLY = INSERT_VEC.TABLE[INSERT_VEC.TABLE$INSERT_VEC >= HIGH_ONLY,]
# INSERT_VEC.TABLE.HIGH_ONLY
NROW.INSERT_VEC.TABLE.HIGH_ONLY = nrow(INSERT_VEC.TABLE.HIGH_ONLY)
for (ROW in 1:NROW.INSERT_VEC.TABLE.HIGH_ONLY){
  # ROW = 1
  COVERED = COVERED + INSERT_VEC.TABLE.HIGH_ONLY[ROW,]$Freq
  PCT_COVERED = COVERED / READ_PAIRS

  if (PCT_COVERED >= 0.1 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_10_PERCENT"] == 0){
    DISTANCE = (INSERT_VEC.TABLE.HIGH_ONLY[ROW,]$INSERT_VEC - MEDIAN)*2 + 1
    WIDTH_OF_X_PCTG_LIST["WIDTH_OF_10_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.2 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_20_PERCENT"] == 0){
    DISTANCE = (INSERT_VEC.TABLE.HIGH_ONLY[ROW,]$INSERT_VEC - MEDIAN)*2 + 1
    WIDTH_OF_X_PCTG_LIST["WIDTH_OF_20_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.3 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_30_PERCENT"] == 0){
    DISTANCE = (INSERT_VEC.TABLE.HIGH_ONLY[ROW,]$INSERT_VEC - MEDIAN)*2 + 1
    WIDTH_OF_X_PCTG_LIST["WIDTH_OF_30_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.4 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_40_PERCENT"] == 0){
    DISTANCE = (INSERT_VEC.TABLE.HIGH_ONLY[ROW,]$INSERT_VEC - MEDIAN)*2 + 1
    WIDTH_OF_X_PCTG_LIST["WIDTH_OF_40_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.5 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_50_PERCENT"] == 0){
    DISTANCE = (INSERT_VEC.TABLE.HIGH_ONLY[ROW,]$INSERT_VEC - MEDIAN)*2 + 1
    WIDTH_OF_X_PCTG_LIST["WIDTH_OF_50_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.6 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_60_PERCENT"] == 0){
    DISTANCE = (INSERT_VEC.TABLE.HIGH_ONLY[ROW,]$INSERT_VEC - MEDIAN)*2 + 1
    WIDTH_OF_X_PCTG_LIST["WIDTH_OF_60_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.7 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_70_PERCENT"] == 0){
    DISTANCE = (INSERT_VEC.TABLE.HIGH_ONLY[ROW,]$INSERT_VEC - MEDIAN)*2 + 1
    WIDTH_OF_X_PCTG_LIST["WIDTH_OF_70_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.8 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_80_PERCENT"] == 0){
    DISTANCE = (INSERT_VEC.TABLE.HIGH_ONLY[ROW,]$INSERT_VEC - MEDIAN)*2 + 1
    WIDTH_OF_X_PCTG_LIST["WIDTH_OF_80_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.9 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_90_PERCENT"] == 0){
    DISTANCE = (INSERT_VEC.TABLE.HIGH_ONLY[ROW,]$INSERT_VEC - MEDIAN)*2 + 1
    WIDTH_OF_X_PCTG_LIST["WIDTH_OF_90_PERCENT"] = DISTANCE }
  # if (PCT_COVERED >= 0.95 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_95_PERCENT"] == 0){
  #   DISTANCE = (INSERT_VEC.TABLE.HIGH_ONLY[ROW,]$INSERT_VEC - MEDIAN)*2 + 1
  #   WIDTH_OF_X_PCTG_LIST["WIDTH_OF_95_PERCENT"] = DISTANCE }
  if (PCT_COVERED >= 0.99 && WIDTH_OF_X_PCTG_LIST["WIDTH_OF_99_PERCENT"] == 0){
    DISTANCE = (INSERT_VEC.TABLE.HIGH_ONLY[ROW,]$INSERT_VEC - MEDIAN)*2 + 1
    WIDTH_OF_X_PCTG_LIST["WIDTH_OF_99_PERCENT"] = DISTANCE }

}
# WIDTH_OF_X_PCTG_LIST

OUT.DF1 =
  data.frame(metrics = c("MEDIAN_INSERT_SIZE",
                         "MEDIAN_ABSOLUTE_DEVIATION",
                         "MIN_INSERT_SIZE",
                         "MAX_INSERT_SIZE",
                         "MEAN_INSERT_SIZE",
                         "STANDARD_DEVIATION",
                         "READ_PAIRS",
                         "PAIR_ORIENTATION"),
             values = c(MEDIAN,
                        MAD,
                        MIN,
                        MAX,
                        round(MEAN,5),
                        round(STDV,5),
                        READ_PAIRS,
                        args$orientation))

OUT.DF2 = data.frame(metrics = names(WIDTH_OF_X_PCTG_LIST),
                     values = do.call(rbind, WIDTH_OF_X_PCTG_LIST)[,1], row.names = NULL)

write.table(rbind(OUT.DF1, OUT.DF2), file = args$output,
            quote = F, sep = '\t', row.names = F, col.names = F)


# EOF
