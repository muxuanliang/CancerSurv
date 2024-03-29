stan=function(coef){
  return(coef/abs(coef[1]))
}

# loss gets the loss of the chosen loss
loss <- function(x, loss_type){
  switch(loss_type,
         logistic = log(1+exp(-x)),
         exponential = exp(-x),
         hinge = (1-x)*((1-x)>=0),
         zeroOne = (x<0)
  )
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
        if (progressionTime < treatYear_tmp){
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
        if (progressionTime < treatYear_tmp){
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
  function(fit = NULL, coef = NULL, threshold = NULL,
           subjectId = NULL,
           timeToEvent = NULL,
           measureTime = NULL,
           measureTimeDiscrete = NULL,
           covariate = NULL,
           censoringIndicator = NULL,
           tau0 = 6,
           timePoints,
           methodFitCensoring = "km", time_varying=FALSE) {
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
    if (!time_varying){
      if (is.null(fit)){
        dataFrame$estimatedDecision <-
          covariate %*% coef > threshold
      } else {
        dataFrame$estimatedDecision <-
          covariate %*% fit$coef[-1] + fit$coef[1] > 0
      }
    } else {
      for (i in 1:nrow(covariate)){
        index <- which(measureTimeDiscrete[i]==timePoints)
        dataFrame$estimatedDecision[i] <-
          covariate[i,] %*% fit[[index]]$coef[-1] + fit[[index]]$coef[1] > 0
      }
    }


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

get_tpr_tnr_new <-function(fit = NULL,
                           biopsyTable.subjectId, biopsyTable.measureTime, biopsyTable.eventTime,
                           PSATable.covariate, PSATable.subjectId, PSATable.measureTime, PSATable.measureTimeDiscrete=NULL,
                           tau0=0.5, time_varying=FALSE, timePoints=NULL, opt=NULL){
  numberCovariate <- NCOL(PSATable.covariate)
  if (is.null(colnames(PSATable.covariate))) {
    colnames(PSATable.covariate) <- paste0("v", 1:(NCOL(PSATable.covariate)))
  }
  dataFrame <-
    data.frame(
      covariate = PSATable.covariate,
      subjectId = PSATable.subjectId,
      measureTime = PSATable.measureTime
    )
  if(!is.null(opt)){
    PSATable.covariate <- t(apply(PSATable.covariate,1,function(t){(t-opt$mean)/opt$sd}))
  }

  if (!time_varying){
    dataFrame$estimatedDecision <-
      (PSATable.covariate %*% fit$coef[-1] + fit$coef[1]) > 0
  } else {
    for (i in 1:nrow(PSATable.covariate)){
      if (PSATable.measureTimeDiscrete[i]>max(timePoints)){
        PSATable.measureTimeDiscrete[i] <- max(timePoints)
      }
      if (PSATable.measureTimeDiscrete[i]< min(timePoints)){
        PSATable.measureTimeDiscrete[i] <- min(timePoints)
      }
      index <- which(PSATable.measureTimeDiscrete[i]==timePoints)
      dataFrame$estimatedDecision[i] <-
        (PSATable.covariate[i,] %*% fit[[index]]$coef[-1] + fit[[index]]$coef[1]) > 0
    }
  }

  # get biopsyTable
  biopsyTable <- data.frame(subjectId=biopsyTable.subjectId, measureTime=biopsyTable.measureTime,
                            is.positive=abs(biopsyTable.measureTime-biopsyTable.eventTime)<0.1,
                            discard=biopsyTable.measureTime > biopsyTable.eventTime,
                            decision=NA)

  tmp_biopsy_table <- biopsyTable
  TPR <- NULL
  TNR <- NULL
  count <- 0
  negative_predict <- 0
  for (time in timePoints){
    pidList <- unique(tmp_biopsy_table$subjectId)
    for (pid in pidList){
      boolVector <- (dataFrame$measureTime[dataFrame$subjectId==pid]-time)<=tau0 & (dataFrame$measureTime[dataFrame$subjectId==pid]-time)>0
      if(any(boolVector)){
        select.index <- min(which(boolVector))
        selected.predict <- dataFrame$estimatedDecision[dataFrame$subjectId==pid][select.index]
        tmp_biopsy_table$decision[tmp_biopsy_table$subjectId==pid] <- selected.predict
        negative_predict <- negative_predict+(!selected.predict)
        count <- count+1
      } else {
        tmp_biopsy_table$discard[tmp_biopsy_table$subjectId==pid] <- TRUE
        next
      }
    }
    tmp_biopsy_table <- tmp_biopsy_table[!tmp_biopsy_table$discard,]

    nom <- ks(tmp_biopsy_table$measureTime, as.numeric(tmp_biopsy_table$decision*(tmp_biopsy_table$is.positive)), time+tau0)-ks(tmp_biopsy_table$measureTime, as.numeric(tmp_biopsy_table$decision*(tmp_biopsy_table$is.positive)), time)
    denom <- ks(tmp_biopsy_table$measureTime, as.numeric(tmp_biopsy_table$is.positive), time+tau0)-ks(tmp_biopsy_table$measureTime, as.numeric(tmp_biopsy_table$is.positive), time)
    TPR <- c(TPR, nom/denom)

    nom <- negative_predict/count-ks(tmp_biopsy_table$measureTime, as.numeric((!tmp_biopsy_table$decision)*(tmp_biopsy_table$is.positive)), time+tau0)
    denom <- 1-ks(tmp_biopsy_table$measureTime, as.numeric(tmp_biopsy_table$is.positive), time+tau0)
    TNR <- c(TNR, nom/denom)
  }

  data.frame(time=timePoints, TPR=TPR, TNR=TNR)
}

# ks gets the kernel estimation
ks <- function(xx, yy, xx.test){
  nobs <- nrow(as.matrix(xx))
  nvars <- ncol(as.matrix(xx))
  hopt <- (4/(nvars+2))^(1/(nvars+4)) * (nobs^(-1/(nvars+4)))
  wm <- function(t){
    if (ncol(as.matrix(xx))==1){
      weight <- exp(-0.5 * (as.numeric(t)-xx)^2/(hopt^2)) * hopt
    } else {
      weight <- apply(xx,1,function(s){exp(-0.5 * sum((t-s)^2)/(hopt^2)) * hopt^(ncol(xx))})
    }
    weighted.mean(yy, weight)
  }
  if (nrow(as.matrix(xx.test))==1) {
    yy.test <- wm(xx.test)
  } else {
    if (ncol((as.matrix(xx.test)))==1){
      yy.test <- sapply(as.matrix(xx.test), function(t){
        wm(t)
      })
    } else {
      yy.test <- apply(as.matrix(xx.test),1,function(t){
        wm(t)
      })
    }
  }
  yy.test
}

ks_sum <- function(xx, yy, xx.test){
  nobs <- nrow(as.matrix(xx))
  nvars <- ncol(as.matrix(xx))
  hopt <- (4/(nvars+2))^(1/(nvars+4)) * (nobs^(-1/(nvars+4)))
  wm <- function(t){
    if (ncol(as.matrix(xx))==1){
      weight <- exp(-0.5 * (as.numeric(t)-xx)^2/(hopt^2)) * hopt
    } else {
      weight <- apply(xx,1,function(s){exp(-0.5 * sum((t-s)^2)/(hopt^2)) * hopt^(ncol(xx))})
    }
    sum(yy * weight)
  }
  if (nrow(as.matrix(xx.test))==1) {
    yy.test <- wm(xx.test)
  } else {
    if (ncol((as.matrix(xx.test)))==1){
      yy.test <- sapply(as.matrix(xx.test), function(t){
        wm(t)
      })
    } else {
      yy.test <- apply(as.matrix(xx.test),1,function(t){
        wm(t)
      })
    }
  }
  yy.test
}

# ks gets the kernel estimation
ks_weight <- function(xx, xx.test, h){
  nobs <- nrow(as.matrix(xx))
  nvars <- ncol(as.matrix(xx))
  hopt <- (4/(nvars+2))^(1/(nvars+4)) * (nobs^(-1/(nvars+2)))
  wm <- function(t){
    if (ncol(as.matrix(xx))==1){
      weight <- exp(-0.5 * (as.numeric(t)-xx)^2/(hopt^2)) * hopt
    } else {
      weight <- apply(xx,1,function(s){exp(-0.5 * sum((t-s)^2)/(hopt^2)) * hopt^(ncol(xx))})
    }
    weight
  }
  if (nrow(as.matrix(xx.test))==1) {
    weight <- wm(xx.test)
  } else {
    if (ncol((as.matrix(xx.test)))==1){
      weight <- sapply(as.matrix(xx.test), function(t){
        wm(t)
      })
    } else {
      weight <- apply(as.matrix(xx.test),1,function(t){
        wm(t)
      })
    }
  }
  weight
}

sampleBasedOnPid <- function(pidList, sampleRatio = 0.5){
  uniquePid <- unique(pidList)
  sampledPid <- sample(uniquePid, size=floor(sampleRatio*length(uniquePid)))

  return(list(sampledPid=sampledPid, restPid=uniquePid[!(uniquePid %in% sampledPid)]))
}

BootstrapBasedOnPid <- function(pidList){
  uniquePid <- unique(pidList)
  sampledPid <- sample(uniquePid, size=length(uniquePid), replace = TRUE)

  return(sampledPid)
}

get_saved_missed_biopsy <- function(fit = NULL,
                                    subjectId = NULL,
                                    timeToEvent = NULL,
                                    measureTime = NULL,
                                    measureTimeDiscrete = NULL,
                                    covariate = NULL,
                                    censoringIndicator = NULL,
                                    positive_biopsy = NULL,
                                    negative_biopsy = NULL,
                                    tau0 = 6,
                                    time_varying = FALSE,
                                    timePoints, opt=NULL){
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
      censoringIndicator = censoringIndicator,
      positive_biopsy = positive_biopsy,
      negative_biopsy = negative_biopsy
    )
  if(!is.null(opt)){
    covariate <- t(apply(covariate,1,function(t){(t-opt$mean)/opt$sd}))
  }
  if (!time_varying){
    dataFrame$estimatedDecision <-
      covariate %*% fit$coef[-1] + fit$coef[1] > 0
  } else {
    for (i in 1:nrow(covariate)){
      index <- which(measureTimeDiscrete[i]==timePoints)
      dataFrame$estimatedDecision[i] <-
        covariate[i,] %*% fit[[index]]$coef[-1] + fit[[index]]$coef[1] > 0
    }
  }

  count_saved_negative <- count_missed_positive <- positive_biopsy <- negative_biopsy <- 0
  data_weighted <- NULL

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

    km <-
      survival::survfit(survival::Surv(sub_data$timeToEvent, censored) ~ 1)
    survest <- stepfun(km$time, c(1, km$surv))
    sub_data$ipcw[sub_data$indicator == 1] = survest(sub_data$timeToEvent[sub_data$indicator ==
                                                                            1])
    sub_data$ipcw[sub_data$indicator == 0] = survest(tau0 + dataSelect$measureTime[sub_data$indicator ==
                                                                                     0])

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

    count_missed_positive <- count_missed_positive + sum(w_positive*(sub_data$estimatedDecision==FALSE))
    count_saved_negative <- count_saved_negative + sum(w_negative*(sub_data$estimatedDecision==FALSE))

    sub_data$missed_biopsy <- (sub_data$estimatedDecision==FALSE) & (ifelse((sub_data$timeToEvent - sub_data$measureTime) <= tau0, 1, 0)==1)
    sub_data$saved_biopsy <- (sub_data$estimatedDecision==FALSE) & (ifelse((sub_data$timeToEvent - sub_data$measureTime) > tau0, 1, 0)==1)

    positive_biopsy <- positive_biopsy + sum(w_positive)
    negative_biopsy <- negative_biopsy + sum(w_negative)

    data_weighted <- rbind(data_weighted, sub_data)
  }
  list(negative_biopsy=negative_biopsy, saved_biopsy=count_saved_negative, positive_biopsy=positive_biopsy, missed_biopsy=count_missed_positive, data_weighted=data_weighted)
}

get_any_biopsy <- function(data_select_discrete, data_biopsy, tau0=0.5){
  negative_biopsy <- positive_biopsy <- NULL
  for (i in 1:nrow(data_select_discrete)){
    pid <- data_select_discrete$CISNET_ID[i]
    biopsy_time <- data_biopsy$TimeSince_Dx[data_biopsy$CISNET_ID==pid]
    if (any((biopsy_time >=data_select_discrete$timePoint[i]) & (biopsy_time < data_select_discrete$timePoint[i]+tau0))){
      index <- min(which(biopsy_time >=data_select_discrete$timePoint[i]))
      negative_biopsy[i] <- data_biopsy$negative_biopsy[data_biopsy$CISNET_ID==pid][index]
      positive_biopsy[i] <- data_biopsy$positive_biopsy[data_biopsy$CISNET_ID==pid][index]
    } else {
      negative_biopsy[i] <- positive_biopsy[i] <- 0
    }
  }
  list(negative_biopsy=negative_biopsy, positive_biopsy=positive_biopsy)
}

get_mean_value <- function(pid, value){
  res <- array(0, length(pid))
  uniquePid <- unique(pid)
  for (id in uniquePid){
    res[pid==id] <- mean(value[pid==id], na.rm=TRUE)
  }
  res
}

get_saved_missed_biopsy_backward <- function(fit = NULL,
                                    biopsyTable.subjectId, biopsyTable.measureTime, biopsyTable.isPositive,
                                    PSATable.covariate, PSATable.subjectId, PSATable.measureTime, PSATable.measureTimeDiscrete=NULL, tau0=0.5, time_varying=FALSE, timePoints=NULL){
  numberCovariate <- NCOL(PSATable.covariate)
  if (is.null(colnames(PSATable.covariate))) {
    colnames(PSATable.covariate) <- paste0("v", 1:(NCOL(PSATable.covariate)))
  }
  dataFrame <-
    data.frame(
      covariate = PSATable.covariate,
      subjectId = PSATable.subjectId,
      measureTime = PSATable.measureTime
    )
  if (!time_varying){
    dataFrame$estimatedDecision <-
      (PSATable.covariate %*% fit$coef[-1] + fit$coef[1]) > 0
  } else {
    for (i in 1:nrow(PSATable.covariate)){
      if (PSATable.measureTimeDiscrete[i]>max(timePoints)){
        PSATable.measureTimeDiscrete[i] <- max(timePoints)
      }
      if (PSATable.measureTimeDiscrete[i]< min(timePoints)){
        PSATable.measureTimeDiscrete[i] <- min(timePoints)
      }
        index <- which(PSATable.measureTimeDiscrete[i]==timePoints)
        dataFrame$estimatedDecision[i] <-
          (PSATable.covariate[i,] %*% fit[[index]]$coef[-1] + fit[[index]]$coef[1]) > 0
    }
  }

  count_negative <- count_positive <- count_saved_negative <- count_missed_positive <- count_discard <- array(FALSE, c(length(biopsyTable.subjectId), 1))
  for (i in 1:length(biopsyTable.subjectId)){
    which.subject <- biopsyTable.subjectId[i]
    boolVector <- biopsyTable.measureTime[i]-dataFrame$measureTime[dataFrame$subjectId==which.subject]<=tau0 & biopsyTable.measureTime[i]-dataFrame$measureTime[dataFrame$subjectId==which.subject]>0
    if(any(boolVector)){
      select.index <- max(which(boolVector))
      selected.predict <- dataFrame$estimatedDecision[dataFrame$subjectId==which.subject][select.index]
    } else {
      count_discard[i] <- TRUE
      next
    }
    if(biopsyTable.isPositive[i]==FALSE){
      count_negative[i] <- TRUE
    }
    if(biopsyTable.isPositive[i]==TRUE){
      count_positive[i] <- TRUE
    }
    if(biopsyTable.isPositive[i]==FALSE & selected.predict==FALSE){
      count_saved_negative[i] <- TRUE
    }
    if(biopsyTable.isPositive[i]==TRUE & selected.predict==FALSE){
      count_missed_positive[i] <- TRUE
    }
  }
  list(negative_biopsy=count_negative, saved_biopsy=count_saved_negative, positive_biopsy=count_positive, missed_biopsy=count_missed_positive, discard_biopsy=count_discard)
}


get_conflict_decision <- function(fit = NULL,
                                    subjectId = NULL,
                                    timeToEvent = NULL,
                                    measureTime = NULL,
                                    measureTimeDiscrete = NULL,
                                    covariate = NULL,
                                    censoringIndicator = NULL,
                                    positive_biopsy = NULL,
                                    negative_biopsy = NULL,
                                    tau0 = 6,
                                    time_varying = FALSE,
                                    timePoints){
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
      censoringIndicator = censoringIndicator,
      positive_biopsy = positive_biopsy,
      negative_biopsy = negative_biopsy
    )
  if (!time_varying){
    dataFrame$estimatedDecision <-
      covariate %*% fit$coef[-1] + fit$coef[1] > 0
  } else {
    for (i in 1:nrow(covariate)){
      index <- which(measureTimeDiscrete[i]==timePoints)
      dataFrame$estimatedDecision[i] <-
        covariate[i,] %*% fit[[index]]$coef[-1] + fit[[index]]$coef[1] > 0
    }
  }

  dataFrame
}
# calculate_conflict counts two possible conflicts.
calculate_conflict <- function(conflictTable){
  conflict_type_1 <- conflict_type_2 <- total_type_1 <- total_type_2 <-
    missed_due_to_longer <- saved_due_to_longer <- missed_anyway <- saved_anyway <-
    compensate_due_to_longer <- extra_due_to_longer <- wont_miss_anyway <- cannot_save_anyway <- total_biopsy <- 0
  for (i in 1:NROW(conflictTable)){
    if (is.na(conflictTable$estimatedDecision.y[i])){
      next
    } else {
      if (conflictTable$estimatedDecision.y[i]==FALSE){
        time <- conflictTable$measureTimeDiscrete[i]
        pid <- conflictTable$subjectId[i]
        selected_decision <- conflictTable$estimatedDecision.x[conflictTable$subjectId==pid & conflictTable$measureTimeDiscrete %in% c(time,time+0.5)]
        selected_positive_biopsy <- conflictTable$positive_biopsy[conflictTable$subjectId==pid & conflictTable$measureTimeDiscrete %in% c(time,time+0.5)]
        selected_negative_biopsy <- conflictTable$negative_biopsy[conflictTable$subjectId==pid & conflictTable$measureTimeDiscrete %in% c(time,time+0.5)]
        selected_weight <- conflictTable$ipc.x[conflictTable$subjectId==pid & conflictTable$measureTimeDiscrete %in% c(time,time+0.5)]
        if (length(selected_decision)<2) {
          next
        }
        if (any(selected_decision)){
          conflict_type_1 <- conflict_type_1+1
          if (sum(selected_positive_biopsy)>0){
            missed_due_to_longer <- missed_due_to_longer+1/min(selected_weight[selected_positive_biopsy==1])
          }
          if (sum(selected_negative_biopsy)>0){
            saved_due_to_longer <- saved_due_to_longer+1/min(selected_weight[selected_negative_biopsy==1])
          }
        } else{
          if (sum(selected_positive_biopsy)>0){
            missed_anyway <- missed_anyway+1/min(selected_weight[selected_positive_biopsy==1])
          }
          if (sum(selected_negative_biopsy)>0){
            saved_anyway <- saved_anyway+1/min(selected_weight[selected_negative_biopsy==1])
          }
        }
        total_type_1 <- total_type_1+1
      }
      if (conflictTable$estimatedDecision.y[i]==TRUE){
        time <- conflictTable$measureTimeDiscrete[i]
        pid <- conflictTable$subjectId[i]
        selected_decision <- conflictTable$estimatedDecision.x[conflictTable$subjectId==pid & conflictTable$measureTimeDiscrete %in% c(time,time+0.5)]
        selected_positive_biopsy <- conflictTable$positive_biopsy[conflictTable$subjectId==pid & conflictTable$measureTimeDiscrete %in% c(time,time+0.5)]
        selected_negative_biopsy <- conflictTable$negative_biopsy[conflictTable$subjectId==pid & conflictTable$measureTimeDiscrete %in% c(time,time+0.5)]
        selected_weight <- conflictTable$ipc.x[conflictTable$subjectId==pid & conflictTable$measureTimeDiscrete %in% c(time,time+0.5)]
        if (length(selected_decision)<2) {
          next
        }
        if (!any(selected_decision)){
          conflict_type_2 <- conflict_type_2+1
          if (sum(selected_positive_biopsy)>0){
            compensate_due_to_longer <- compensate_due_to_longer+1/min(selected_weight[selected_positive_biopsy==1])
          }
          if (sum(selected_negative_biopsy)>0){
            extra_due_to_longer <- extra_due_to_longer+1/min(selected_weight[selected_negative_biopsy==1])
          }
        } else {
          if (sum(selected_positive_biopsy)>0){
            wont_miss_anyway <- wont_miss_anyway+1/min(selected_weight[selected_positive_biopsy==1])
          }
          if (sum(selected_negative_biopsy)>0){
            cannot_save_anyway <- cannot_save_anyway+1/min(selected_weight[selected_negative_biopsy==1])
          }
        }
        total_type_2 <- total_type_2+1
      }
    }
  }
  return(list(conflict_type_1=conflict_type_1, conflict_type_2=conflict_type_2,
              total_type_1=total_type_1, total_type_2=total_type_2,
              missed_due_to_longer=missed_due_to_longer, saved_due_to_longer=saved_due_to_longer,
              missed_anyway=missed_anyway, saved_anyway=saved_anyway,
              compensate_due_to_longer=compensate_due_to_longer, extra_due_to_longer=extra_due_to_longer,
              wont_miss_anyway=wont_miss_anyway, cannot_save_anyway=cannot_save_anyway,
              total_negative =sum(conflictTable$negative_biopsy), total_positive=sum(conflictTable$positive_biopsy)))
}
