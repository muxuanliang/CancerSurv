stan=function(coef){
  return(coef/abs(coef[1]))
}

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

time_standize <- function(vectorOfTime, standardizedTime){
  sapply(vectorOfTime, function(t){
    if (any(which((t-standardizedTime)>=0))){
      standardizedTime[max(which((t-standardizedTime)>=0))]
    } else {
      t
    }
  })
}

get_event_time <- function(pid, time, gleasonScores, treatYear){
  uniquePid <- unique(pid)
  eventTime <- array(0, c(length(uniquePid),1))
  eventIndicator <- array(0, c(length(uniquePid),1))
  gleasonScores[is.na(gleasonScores)] <- 3
  for (index in 1:length(uniquePid)){
    tmp <- uniquePid[index]
    if (NROW(gleasonScores[pid==tmp,,drop=FALSE])==1){
      progressionTime <- time[pid==tmp]
      treatYear_tmp <- treatYear[pid==tmp]
      if (is.na(treatYear_tmp)){
        eventTime[index] <- progressionTime
        eventIndicator[index] <- ifelse(sum(gleasonScores[pid==tmp,])>6,1,0)
      } else {
        if (progressionTime <= treatYear_tmp){
          eventTime[index] <- progressionTime
          eventIndicator[index] <- 1
        } else {
          eventTime[index] <- treatYear_tmp
          eventIndicator[index] <- 0
        }
      }
    } else {
      progressionTime <- ifelse(any(apply(gleasonScores[pid==tmp,],1,sum)>6),
                                time[pid==tmp][min(which(apply(gleasonScores[pid==tmp,],1,sum)>6))],
                                max(time[pid==tmp]))
      treatYear_tmp <- min(unique(treatYear[pid==tmp]))
      if (is.na(treatYear_tmp)){
        eventTime[index] <- progressionTime
        eventIndicator[index] <- ifelse(any(apply(gleasonScores[pid==tmp,],1,sum)>6),1,0)
      } else {
        if (progressionTime <= treatYear_tmp){
          eventTime[index] <- progressionTime
          eventIndicator[index] <- 1
        } else {
          eventTime[index] <- treatYear_tmp
          eventIndicator[index] <- 0
        }
      }
    }
  }
  return(data.frame(pid = uniquePid, eventTime = eventTime, eventIndicator = eventIndicator))
}

get_max_core_ratio <- function(pid,cores_ratio){
  max_cores_ratio <- array(0, c(length(cores_ratio),1))
  uniquePid <- unique(pid)
  for(tmp in uniquePid){
    max_cores_ratio[pid==tmp] <- get_max_core_ratio_per_patient(cores_ratio[pid==tmp])
  }
  max_cores_ratio
}

get_max_core_ratio_per_patient <- function(cores_ratio){
  max_cores_ratio <- array(0, c(length(cores_ratio),1))
  for (index in 1:length(cores_ratio)){
    max_cores_ratio[index] <- max(c(cores_ratio[1:index], 0), na.rm = TRUE)
  }
  max_cores_ratio
}

get_tpr_tnr <- function(fit=NULL, subjectId=NULL, timeToEvent=NULL, measureTime=NULL, covariate=NULL, censoringIndicator=NULL, tau0=6, timePoint, methodFitCensoring="km"){
  numberCovariate <- NCOL(covariate)
  if (is.null(colnames(covariate))){
    colnames(covariate) <- paste0("v", 1:(NCOL(covariate)))
  }
  dataFrame <- data.frame(covariate=covariate, subjectId=subjectId, timeToEvent=timeToEvent, measureTime=measureTime, censoringIndicator=censoringIndicator)

  # select the covariate to predict the decision at timePoint (recent covariate for each subject)
  uniquePid <- unique(dataFrame$subjectId)
  dataSelect <- NULL
  for (pid in uniquePid){
    idx_tmp <- max(which(dataFrame$measureTime[dataFrame$subjectId==pid]-timePoint<= tau0))
    cdata_tmp <- dataFrame[dataFrame$subjectId==pid,][idx_tmp,]
    dataSelect <- rbind(dataSelect, data_tmp)
  }
  estimatedDecision <- dataSelect$covariate %*% fit$coef[-1] + fit$coef[1] > 0


  dataSelect$indicator <- ifelse(dataSelect$timeToEvent-timePoint<=tau0, 1,0)
  ipc <- NULL
  sub_data=dataSelect[dataSelect$measureTime>=timePoint,]
  sub_data$ipcw=rep(0, nrow(sub_data))
  censored= 1-sub_data$censoringIndicator
  if (methodFitCensoring=="km"){
    km <- survival::survfit(survival::Surv(sub_data$timeToEvent,censored)~1)
    survest <- stepfun(km$time, c(1, km$surv))
    sub_data$ipcw[sub_data$indicator==1]=survest(sub_data$timeToEvent[sub_data$indicator==1])
    sub_data$ipcw[sub_data$indicator==0]=survest(tau0+timePoint)
  } else if (methodFitCensoring=="cox"){
    model=paste0("coxfit <- survival::coxph(survival::Surv(sub_data$timeToEvent,censored)~", paste(colnames(covariate), collapse=" + "),", data=sub_data)" )
    eval(parse(text=model))
    sub_data$ipcw[sub_data$indicator==1]=predict(coxfit, type="survival")[sub_data$indicator==1]
    sub_data$ipcw[sub_data$indicator==0]=predict(coxfit, newdata=data.frame(sub_data[sub_data$indicator==0, 1:numberCovariate],timeToEvent=rep(tau0+timePoint, times=sum(sub_data$indicator==0))), type="survival")
  }
  sub_data$ipc <- sub_data$ipcw

  denomm <- NULL
  idx <- which(sub_data$censoringIndicator==1|(sub_data$censoringIndicator==0&sub_data$indicator==0))
  sub_data <- sub_data[idx,]

  w_positive <- ifelse((sub_data$timeToEvent-timePoint)<=tau0, 1, 0)/sub_data$ipc
  w_negative <- ifelse((sub_data$timeToEvent-timePoint)>tau0, 1, 0)/sub_data$ipc

  tpr <-


  return(res)
}

