cancerSurvInterval <-function(biopsyTable.subjectId, biopsyTable.measureTime, biopsyTable.eventTime,
                              biopsyTable.primary_gleason, biopsyTable.secondary_gleason,
                           PSATable.covariate, PSATable.subjectId, PSATable.measureTime,
                           PSATable.measureTimeDiscrete=NULL, tau0=0.5, time_varying=FALSE,
                           timePoints=NULL, xi=10, lambdaSeq=NULL, nlambda=100, n.fold=2,
                           standardize = FALSE, hopt = 0.2*sd(biopsyTable$measureTime) *
                             (length(unique(PSATable.subjectId)))^(-1/5),
                           shift_fraction = 0.5, prevalence=NULL){


  # get prevalence
  if (is.null(prevalence)){
    prevalence <- get_prevalence(biopsyTable.subjectId,
                                 biopsyTable.measureTime,
                                 biopsyTable.eventTime,
                                 biopsyTable.primary_gleason,
                                 biopsyTable.secondary_gleason,
                                 tau0=tau0,
                                 timePoints = timePoints)
  }
  # fit
  numberCovariate <- NCOL(PSATable.covariate)
  if (is.null(colnames(PSATable.covariate))) {
    colnames(PSATable.covariate) <- paste0("v", 1:(NCOL(PSATable.covariate)))
  }
  time_sd <- sd(PSATable.measureTime)
  if (standardize){
    covariate.mean <- apply(PSATable.covariate,2,mean)
    covariate.sd <- apply(PSATable.covariate,2,sd)
    PSATable.covariate <- apply(PSATable.covariate,2,function(t){(t-mean(t))/sd(t)})
  } else {
    covariate.mean <- apply(PSATable.covariate,2,mean)*0
    covariate.sd <- apply(PSATable.covariate,2,sd)*0+1
  }

  dataFrame <-
    data.frame(
      covariate = PSATable.covariate,
      subjectId = PSATable.subjectId,
      measureTime = PSATable.measureTime
    )

  # get biopsyTable
  biopsyTable <- data.frame(subjectId=biopsyTable.subjectId, measureTime=biopsyTable.measureTime,
                            is.positive=(biopsyTable.primary_gleason+biopsyTable.secondary_gleason)>6,
                            discard=biopsyTable.measureTime > biopsyTable.eventTime,
                            decision=NA)

  pseudo_covariate <- pseudo_time <- pseudo_positive_weight <- pseudo_negative_weight <- trade_off_weight <- NULL
  for (time in timePoints){
    tmp_biopsy_table <- biopsyTable
    tmp_biopsy_table$last_biopsy_time <- NULL
    tmp_biopsy_table$final_biopsy <- FALSE
    for (index in 1:NROW(tmp_biopsy_table)){
      tmp_pid <- tmp_biopsy_table$subjectId[index]
      if (index==1){
        tmp_biopsy_table$last_biopsy_time[index] <- 0
      } else if (tmp_biopsy_table$subjectId[index-1]!=tmp_pid){
        tmp_biopsy_table$last_biopsy_time[index] <- 0
      } else {
        tmp_biopsy_table$last_biopsy_time[index] <- tmp_biopsy_table$measureTime[index-1]
      }

      if (tmp_biopsy_table$measureTime[index]==max(tmp_biopsy_table$measureTime[tmp_biopsy_table$subjectId==tmp_pid])){
        tmp_biopsy_table$final_biopsy[index] <- TRUE
      }
    }

    sd_scale <- sd(tmp_biopsy_table$measureTime)
    pidList <- unique(tmp_biopsy_table$subjectId)

    for (pid in pidList){
      boolVector <- (dataFrame$measureTime[dataFrame$subjectId==pid]-time)<=tau0 & (dataFrame$measureTime[dataFrame$subjectId==pid]-time)>0
      if(any(boolVector)){
        select.index <- min(which(boolVector))
        selected.covariate <- PSATable.covariate[dataFrame$subjectId==pid,,drop=FALSE][select.index,,drop=FALSE]
      } else {
        next
      }

    select_biopsy_table <- tmp_biopsy_table[tmp_biopsy_table$subjectId==pid & tmp_biopsy_table$is.positive==TRUE,]

    sd_scale <- sd(tmp_biopsy_table$measureTime)
    hopt_p <- sd_scale* length(unique(tmp_biopsy_table$subjectId))^{-1/6}*100 #6
    positive_weight <- ks_sum(cbind(select_biopsy_table$last_biopsy_time,select_biopsy_table$measureTime), as.numeric(select_biopsy_table$is.positive), cbind(time,time+tau0), hopt = hopt_p)/ks_sum(cbind(tmp_biopsy_table$last_biopsy_time,tmp_biopsy_table$measureTime), as.numeric(tmp_biopsy_table$is.positive), cbind(time,time+tau0), hopt = hopt_p)

    select_biopsy_table <- tmp_biopsy_table[tmp_biopsy_table$subjectId==pid & tmp_biopsy_table$is.positive==FALSE,]
    hopt_n <- sd_scale* length(unique(tmp_biopsy_table$subjectId))^{-1/5}*100 # 5
    negative_weight <- ks_sum(select_biopsy_table$measureTime, as.numeric(!select_biopsy_table$is.positive), time+tau0, hopt_n)/ks_sum(tmp_biopsy_table$measureTime, as.numeric(!tmp_biopsy_table$is.positive), time+tau0, hopt_n)

    pseudo_covariate <- rbind(pseudo_covariate,selected.covariate)
    pseudo_time <- c(pseudo_time, time)
    pseudo_positive_weight <- c(pseudo_positive_weight, positive_weight)
    pseudo_negative_weight <- c(pseudo_negative_weight, negative_weight)

    est_prevalence <- prevalence[which(time==timePoints)]
    trade_off_weight <- c(trade_off_weight, (1-est_prevalence)/est_prevalence*xi)
    }
  }

  # pseudo_data
  pseudo_data <- list(label=sign(pseudo_positive_weight-trade_off_weight*pseudo_negative_weight),
                      covariate=pseudo_covariate,
                      weight=abs(pseudo_positive_weight-trade_off_weight*pseudo_negative_weight))

  shift <- shift_fraction*(pseudo_positive_weight+trade_off_weight*pseudo_negative_weight)
  pseudo_positive_weight_shifted <- (pseudo_positive_weight-shift) * (pseudo_positive_weight-shift >= 0) -
    (trade_off_weight*pseudo_negative_weight-shift) * (trade_off_weight*pseudo_negative_weight-shift <= 0)
  pseudo_negative_weight_shifted <- -(pseudo_positive_weight-shift) * (pseudo_positive_weight-shift <= 0) +
    (trade_off_weight*pseudo_negative_weight-shift) * (trade_off_weight*pseudo_negative_weight-shift >= 0)
  pseudo_data <- list(label=c(array(1, length(pseudo_positive_weight_shifted)), array(-1, length(trade_off_weight*pseudo_negative_weight_shifted))),
                      covariate=rbind(pseudo_covariate, pseudo_covariate),
                      weight=c(pseudo_positive_weight_shifted, trade_off_weight*pseudo_negative_weight_shifted))
  covariate=pseudo_data$covariate
  outt=cbind(pseudo_data$weight, pseudo_data$label, covariate)
  colnames(outt)=c("outcome", "A", paste0("v", 1:(NCOL(covariate))))
  outt=as.data.frame(outt)
  outt[,2]=as.character(pseudo_data$label)
  moPropen <- modelObj::buildModelObj(model = ~1,solver.method = 'glm',
                                      solver.args = list('family'='binomial'),
                                      predict.method = 'predict.glm',
                                      predict.args = list(type='response'))

  tmp_text <- "v1"
  if (NCOL(covariate)>=2){
    for (iter in 2:(NCOL(covariate))){
      tmp_text <- paste0(tmp_text, "+v", iter)
    }
  }
  expr <- paste0("earlRes=try(DynTxRegime::owl(moPropen = moPropen,
                    data = outt, reward=outt$outcome/mean(outt$outcome), txName = 'A',
                    regime = ~ ",tmp_text,",lambdas =10^(seq(-2,1,length.out=11)), cvFolds = 2,
                    kernel = 'linear', surrogate = 'logit',kparam = NULL), TRUE)")
  eval(expr = parse(text=expr))
  coef_osf = DynTxRegime::regimeCoef(earlRes)

  return(list(fit=coef_osf,opt=list(mean=covariate.mean, sd=covariate.sd), prevalence=prevalence))
}
