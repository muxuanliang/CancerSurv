setwd("~/Documents/Research/Yingqi Zhao/CancerSuv")

library(CancerSurv)
library(readr)
library(dplyr)

data_raw <- read_excel_allsheets("PASSData.xlsx", tibble = FALSE)
data_old <- read_csv("PASS.csv")
data_old_stat <- data_old[,is.na(apply(data_old[data_old$pid==1,],2,sd))|(apply(data_old[data_old$pid==1,],2,sd)==0)]
colnamesToMerge <- c("pid", "ppt_cancerhx_agedx_value", "study_site_id","log_avg_size_fillin", "famhx_ever", "event_grade") # smoking status, bmi
data_old_stat <- data_old_stat[,colnamesToMerge]
data_tmp <- (distinct(data_old_stat))
if (length(unique(data_old_stat$pid))==NROW(data_tmp)){
  print("Dintinct rows")
} else {
  print("Duplicate rows. Need Check!")
}

# merge static tables
colnames(data_tmp)[1] <- "CISNET_ID"
data_static <- data_raw[[1]]
data_static <- data_static %>%  left_join(data_tmp, by = "CISNET_ID", copy=TRUE)

# merge time varying info
## create time x pid stamp
data_dynamic <- data.frame(CISNET_ID=c(data_raw[[2]]$CISNET_ID, data_raw[[3]]$CISNET_ID, data_raw[[4]]$CISNET_ID), TimeSince_Dx=c(data_raw[[2]]$TimeSince_Dx, data_raw[[3]]$TimeSince_Dx, data_raw[[4]]$TimeSince_Dx))
data_dynamic <- distinct(data_dynamic)
## start merge
data_dynamic <- data_dynamic %>%  left_join(data_raw[[2]], by = c("CISNET_ID", "TimeSince_Dx"), copy=TRUE)
data_dynamic <- data_dynamic %>%  left_join(data_raw[[3]], by = c("CISNET_ID", "TimeSince_Dx"), copy=TRUE)
data_dynamic <- data_dynamic %>%  left_join(data_raw[[4]], by = c("CISNET_ID", "TimeSince_Dx"), copy=TRUE)
## discard the data before the first biopsy
data_dynamic <- arrange(data_dynamic, CISNET_ID, TimeSince_Dx)
getPid <- unique(data_dynamic$CISNET_ID)
rowToDiscard <- NULL
for (id in getPid){
  start_row <- min(which(data_dynamic$CISNET_ID==id))
  nrows_to_discard <- which(data_dynamic$BxN[data_dynamic$CISNET_ID==id]==1)-1
  if(nrows_to_discard>0){
    rowToDiscard <- c(rowToDiscard, seq(start_row, start_row+nrows_to_discard-1, by=1))
  }
}
data_dynamic <- data_dynamic[-rowToDiscard,]

# Join with Static info
data <- data_dynamic %>% left_join(data_static, by = "CISNET_ID", copy=TRUE)

# write data
write_csv(
  data,
  "Merged_PASS.csv",
  na = "NA"
)

