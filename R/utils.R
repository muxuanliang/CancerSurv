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

get_event_time <- function(pid, time, gleasonScores, treatYear, timeOnAS){
  uniquePid <- unique(pid)
  eventTime <- array(0, c(length(uniquePid),1))
  eventIndicator <- array(0, c(length(uniquePid),1))
  gleasonScores[is.na(gleasonScores)] <- 3
  for (index in 1:length(uniquePid)){
    tmp <- uniquePid[index]
    if (NROW(gleasonScores[pid==tmp,,drop=FALSE])==1){
      progressionTime <- unique(timeOnAS[pid==tmp])
      treatYear_tmp <- treatYear[pid==tmp]
      if (is.na(treatYear_tmp)){
        eventTime[index] <- progressionTime
        eventIndicator[index] <- ifelse(sum(gleasonScores[pid==tmp,])>6,1,0)
        if (eventIndicator[index]==0){
          if (unique(timeOnAS[pid==tmp])-max(time[pid==tmp])>2){
            eventTime[index] <- max(time[pid==tmp])+2
          }
        }
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
                                unique(timeOnAS[pid==tmp]))
      treatYear_tmp <- min(unique(treatYear[pid==tmp]))
      if (is.na(treatYear_tmp)){
        eventTime[index] <- progressionTime
        eventIndicator[index] <- ifelse(any(apply(gleasonScores[pid==tmp,],1,sum)>6),1,0)
        if (eventIndicator[index]==0){
          if (unique(timeOnAS[pid==tmp])-max(time[pid==tmp])>2){
            eventTime[index] <- max(time[pid==tmp])+2
          }
        }
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

get_max_value <- function(pid,value,value_baseline = NA){
  max_value <- array(0, c(length(value),1))
  uniquePid <- unique(pid)
  for(tmp in uniquePid){
    max_value[pid==tmp] <- get_max_value_per_patient(value[pid==tmp], value_baseline[pid==tmp])
  }
  max_value
}

get_max_value_per_patient <- function(value, value_baseline = NA){
  max_value <- array(0, c(length(value),1))
  for (index in 1:length(value)){
    max_value[index] <- max(c(value[1:index], 0, value_baseline), na.rm = TRUE)
  }
  max_value
}

get_recent_value <- function(pid,value,value_baseline = NA){
  recent_value <- rep(0, by=length(value))
  uniquePid <- unique(pid)
  for(tmp in uniquePid){
    if (anyNA(value_baseline[pid==tmp])){
      value_baseline[pid==tmp] <- get_baseline_value_per_patient(value[pid==tmp])
    }
    if (anyNA(value_baseline[pid==tmp])){
      print(tmp)
      next
    }
    recent_value[pid==tmp] <- get_recent_value_per_patient(value[pid==tmp], value_baseline[pid==tmp])
  }
  recent_value
}

get_baseline_value_per_patient <- function(value){
  value_vec_rm_na <- value[which(!is.na(value))]
  baseline_value <- rep(value_vec_rm_na[1], by=length(value))
}

get_recent_value_per_patient <- function(value, value_baseline = NA){
  recent_value <- array(0, c(length(value),1))
  for (index in 1:length(value)){
    value_vec <- c(value_baseline, value[1:index])
    value_vec_rm_na <- value_vec[which(!is.na(value_vec))]
    recent_value[index] <- value_vec_rm_na[length(value_vec_rm_na)]
  }
  recent_value
}

get_tpr_tnr <-
  function(fit = NULL,
           subjectId = NULL,
           timeToEvent = NULL,
           measureTime = NULL,
           measureTimeDiscrete = NULL,
           covariate = NULL,
           censoringIndicator = NULL,
           tau0 = 6,
           timePoints,
           methodFitCensoring = "km") {
    numberCovariate <- NCOL(covariate)
    if (is.null(colnames(covariate))) {
      colnames(covariate) <- paste0("v", 1:(NCOL(covariate)))
    }
    dataFrame <-
      data.frame(
        covariate = covariate,
        subjectId = subjectId,
        timeToEvent = timeToEvent,
        measureTime = measureTime,
        measureTimeDiscrete = measureTimeDiscrete,
        censoringIndicator = censoringIndicator
      )
    dataFrame$estimatedDecision <-
      covariate %*% fit$coef[-1] + fit$coef[1] > 0

    res <- data.frame(tpr = rep(NA, rep=length(timePoints)),
                      tnr = rep(NA, rep=length(timePoints)),
                      timePoint = timePoints)
    for (timePoint in timePoints) {
      # select the covariate to predict the decision at timePoint (recent covariate for each subject)
      uniquePid <- unique(dataFrame$subjectId)
      dataSelect <- NULL
      for (pid in uniquePid) {
        if (any(dataFrame$measureTimeDiscrete[dataFrame$subjectId == pid] == timePoint)) {
          idx_tmp <-
            (which(dataFrame$measureTimeDiscrete[dataFrame$subjectId == pid] == timePoint))
          data_tmp <- dataFrame[dataFrame$subjectId == pid, ][idx_tmp, ]
          dataSelect <- rbind(dataSelect, data_tmp)
        }
      }

      dataSelect$indicator <-
        ifelse(dataSelect$timeToEvent - dataSelect$measureTime <= tau0, 1, 0)
      ipc <- NULL
      sub_data = dataSelect
      sub_data$ipcw = rep(0, nrow(sub_data))
      censored = 1 - sub_data$censoringIndicator
      if (methodFitCensoring == "km") {
        km <-
          survival::survfit(survival::Surv(sub_data$timeToEvent, censored) ~ 1)
        survest <- stepfun(km$time, c(1, km$surv))
        sub_data$ipcw[sub_data$indicator == 1] = survest(sub_data$timeToEvent[sub_data$indicator ==
                                                                                1])
        sub_data$ipcw[sub_data$indicator == 0] = survest(tau0 + dataSelect$measureTime[sub_data$indicator ==
                                                                                         0])
      } else if (methodFitCensoring == "cox") {
        model = paste0(
          "coxfit <- survival::coxph(survival::Surv(sub_data$timeToEvent,censored)~",
          paste(colnames(covariate), collapse = " + "),
          ", data=sub_data)"
        )
        eval(parse(text = model))
        sub_data$ipcw[sub_data$indicator == 1] = predict(coxfit, type = "survival")[sub_data$indicator ==
                                                                                      1]
        sub_data$ipcw[sub_data$indicator == 0] = predict(coxfit,
                                                         newdata = data.frame(sub_data[sub_data$indicator == 0, 1:numberCovariate], timeToEvent =
                                                                                rep(
                                                                                  tau0 + timePoint, times = sum(sub_data$indicator == 0)
                                                                                )),
                                                         type = "survival")
      }
      sub_data$ipc <- sub_data$ipcw

      denomm <- NULL
      idx <-
        which(
          sub_data$censoringIndicator == 1 |
            (sub_data$censoringIndicator == 0 & sub_data$indicator == 0)
        )
      sub_data <- sub_data[idx, ]

      w_positive <-
        ifelse((sub_data$timeToEvent - sub_data$measureTime) <= tau0, 1, 0) / sub_data$ipc
      w_negative <-
        ifelse((sub_data$timeToEvent - sub_data$measureTime) > tau0, 1, 0) / sub_data$ipc

      res$tpr[res$timePoint == timePoint] <-
        mean(w_positive * sub_data$estimatedDecision * sub_data$censoringIndicator) /
        mean(w_positive * sub_data$censoringIndicator)
      res$tnr[res$timePoint == timePoint] <-
        mean(w_negative * (1 - sub_data$estimatedDecision)) / mean(w_negative)
    }
    return(res)
  }

sampleBasedOnPid <- function(pidList, sampleRatio = 0.5){
  uniquePid <- unique(pidList)
  sampledPid <- sample(uniquePid, size=floor(sampleRatio*length(uniquePid)))

  return(list(sampledPid=sampledPid, restPid=uniquePid[!(uniquePid %in% sampledPid)]))
}
