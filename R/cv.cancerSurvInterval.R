cv.cancerSurvInterval <-function(biopsyTable.subjectId, biopsyTable.measureTime, biopsyTable.eventTime,
                              biopsyTable.primary_gleason, biopsyTable.secondary_gleason,
                              PSATable.covariate, PSATable.subjectId, PSATable.measureTime,
                              PSATable.measureTimeDiscrete=NULL, tau0=0.5, time_varying=FALSE,
                              timePoints=NULL, xi=10, lambdaSeq=NULL, nlambda=100, n.fold=2,
                              standardize = FALSE, hopt = 0.2*sd(biopsyTable$measureTime) *
                                (length(unique(PSATable.subjectId)))^(-1/5),
                              number_shift_fractions = 10){

  # cv select shift_fraction
  shift_fraction_seq <- seq(0, 0.5, length.out=number_shift_fractions)
  split_index <- sampleBasedOnPid(PSATable.subjectId)
  fit <- NULL
  val <- NULL
  for (shift_fraction_tmp in shift_fraction_seq){
    fit[[1]] <- cancerSurvInterval(biopsyTable.subjectId[biopsyTable.subjectId %in% split_index$sampledPid],
                                   biopsyTable.measureTime[biopsyTable.subjectId %in% split_index$sampledPid],
                                   biopsyTable.eventTime[biopsyTable.subjectId %in% split_index$sampledPid],
                                   biopsyTable.primary_gleason[biopsyTable.subjectId %in% split_index$sampledPid],
                                   biopsyTable.secondary_gleason[biopsyTable.subjectId %in% split_index$sampledPid],
                                   PSATable.covariate[PSATable.subjectId %in% split_index$sampledPid,],
                                   PSATable.subjectId[PSATable.subjectId %in% split_index$sampledPid],
                                   PSATable.measureTime[PSATable.subjectId %in% split_index$sampledPid],
                                   PSATable.measureTimeDiscrete[PSATable.subjectId %in% split_index$sampledPid],
                                   tau0=tau0, time_varying=time_varying,
                                   timePoints=timePoints, xi=xi, lambdaSeq=lambdaSeq,
                                   nlambda=nlambda, n.fold=n.fold,
                                   standardize = standardize, hopt = hopt,
                                   shift_fraction = shift_fraction_tmp)
    fit[[2]] <- cancerSurvInterval(biopsyTable.subjectId[biopsyTable.subjectId %in% split_index$restPid],
                                   biopsyTable.measureTime[biopsyTable.subjectId %in% split_index$restPid],
                                   biopsyTable.eventTime[biopsyTable.subjectId %in% split_index$restPid],
                                   biopsyTable.primary_gleason[biopsyTable.subjectId %in% split_index$restPid],
                                   biopsyTable.secondary_gleason[biopsyTable.subjectId %in% split_index$restPid],
                                   PSATable.covariate[PSATable.subjectId %in% split_index$restPid,],
                                   PSATable.subjectId[PSATable.subjectId %in% split_index$restPid],
                                   PSATable.measureTime[PSATable.subjectId %in% split_index$restPid],
                                   PSATable.measureTimeDiscrete[PSATable.subjectId %in% split_index$restPid],
                                   tau0=tau0, time_varying=time_varying,
                                   timePoints=timePoints, xi=xi, lambdaSeq=lambdaSeq,
                                   nlambda=nlambda, n.fold=n.fold,
                                   standardize = standardize, hopt = hopt,
                                   shift_fraction = shift_fraction_tmp)

    res_fit1_tmp <- get_tpr_tnr_proposed(fit=list(coef=fit[[1]]$fit), biopsyTable.subjectId[biopsyTable.subjectId %in% split_index$restPid],
                                         biopsyTable.measureTime[biopsyTable.subjectId %in% split_index$restPid],
                                         biopsyTable.eventTime[biopsyTable.subjectId %in% split_index$restPid],
                                         biopsyTable.primary_gleason[biopsyTable.subjectId %in% split_index$restPid],
                                         biopsyTable.secondary_gleason[biopsyTable.subjectId %in% split_index$restPid],
                                         PSATable.covariate[PSATable.subjectId %in% split_index$restPid,],
                                         PSATable.subjectId[PSATable.subjectId %in% split_index$restPid],
                                         PSATable.measureTime[PSATable.subjectId %in% split_index$restPid],
                                         PSATable.measureTimeDiscrete[PSATable.subjectId %in% split_index$restPid],
                                         tau0=tau0,
                                         timePoints = timePoints)
    res_fit2_tmp <- get_tpr_tnr_proposed(fit=list(coef=fit[[2]]$fit), biopsyTable.subjectId[biopsyTable.subjectId %in% split_index$sampledPid],
                                         biopsyTable.measureTime[biopsyTable.subjectId %in% split_index$sampledPid],
                                         biopsyTable.eventTime[biopsyTable.subjectId %in% split_index$sampledPid],
                                         biopsyTable.primary_gleason[biopsyTable.subjectId %in% split_index$sampledPid],
                                         biopsyTable.secondary_gleason[biopsyTable.subjectId %in% split_index$sampledPid],
                                         PSATable.covariate[PSATable.subjectId %in% split_index$sampledPid,],
                                         PSATable.subjectId[PSATable.subjectId %in% split_index$sampledPid],
                                         PSATable.measureTime[PSATable.subjectId %in% split_index$sampledPid],
                                         PSATable.measureTimeDiscrete[PSATable.subjectId %in% split_index$sampledPid],
                                         tau0=tau0,
                                         timePoints = timePoints)
    prevalence_tmp <- (res_fit1_tmp$prevalence+res_fit2_tmp$prevalence)/2
    val <- rbind(val, c(res=(res_fit1_tmp$TPR+xi*(1-prevalence_tmp)/prevalence_tmp*res_fit1_tmp$TNR+
                               res_fit2_tmp$TPR+xi*(1-prevalence_tmp)/prevalence_tmp*res_fit2_tmp$TNR)/2,
                        shift_fraction=shift_fraction_tmp))
  }
  val_sum <- apply(val[,1:length(timePoints)],1,mean)
  shift_fraction_select <- shift_fraction_seq[min(which.max(val_sum))]


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

  pseudo_covariate <- pseudo_time <- pseudo_positive_weight <- pseudo_negative_weight <- trade_off_weight <- prevalence <- NULL
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
    hopt_p <- sd_scale* length(unique(tmp_biopsy_table$subjectId))^{-1/6}
    hopt_n <- sd_scale* length(unique(tmp_biopsy_table$subjectId))^{-1/5}
    est_prevalence <- ks(cbind(tmp_biopsy_table$last_biopsy_time,tmp_biopsy_table$measureTime),
                         as.numeric(tmp_biopsy_table$is.positive), cbind(time,time+tau0), hopt = hopt_p)/
      (ks(cbind(tmp_biopsy_table$last_biopsy_time,tmp_biopsy_table$measureTime), as.numeric(tmp_biopsy_table$is.positive),
          cbind(time,time+tau0), hopt = hopt_p)+
         ks(tmp_biopsy_table$measureTime, as.numeric(!tmp_biopsy_table$is.positive), time+tau0, hopt_n))

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
      hopt_p <- sd_scale* length(unique(tmp_biopsy_table$subjectId))^{-1/6}
      positive_weight <- ks_sum(cbind(select_biopsy_table$last_biopsy_time,select_biopsy_table$measureTime), as.numeric(select_biopsy_table$is.positive), cbind(time,time+tau0), hopt = hopt_p)/ks_sum(cbind(tmp_biopsy_table$last_biopsy_time,tmp_biopsy_table$measureTime), as.numeric(tmp_biopsy_table$is.positive), cbind(time,time+tau0), hopt = hopt_p)

      select_biopsy_table <- tmp_biopsy_table[tmp_biopsy_table$subjectId==pid & tmp_biopsy_table$is.positive==FALSE,]
      hopt_n <- sd_scale* length(unique(tmp_biopsy_table$subjectId))^{-1/5}
      negative_weight <- ks_sum(select_biopsy_table$measureTime, as.numeric(!select_biopsy_table$is.positive), time+tau0, hopt_n)/ks_sum(tmp_biopsy_table$measureTime, as.numeric(!tmp_biopsy_table$is.positive), time+tau0, hopt_n)

      pseudo_covariate <- rbind(pseudo_covariate,selected.covariate)
      pseudo_time <- c(pseudo_time, time)
      pseudo_positive_weight <- c(pseudo_positive_weight, positive_weight)
      pseudo_negative_weight <- c(pseudo_negative_weight, negative_weight)

      trade_off_weight <- c(trade_off_weight, (1-est_prevalence)/est_prevalence*xi)
    }
    prevalence <- c(prevalence, est_prevalence)
  }

  shift <- shift_fraction_select*(pseudo_positive_weight+trade_off_weight*pseudo_negative_weight)
  pseudo_positive_weight_shifted <- (pseudo_positive_weight-shift) * (pseudo_positive_weight-shift >= 0) -
    (trade_off_weight*pseudo_negative_weight-shift) * (trade_off_weight*pseudo_negative_weight-shift <= 0)
  pseudo_negative_weight_shifted <- -(pseudo_positive_weight-shift) * (pseudo_positive_weight-shift <= 0) +
    (trade_off_weight*pseudo_negative_weight-shift) * (trade_off_weight*pseudo_negative_weight-shift >= 0)
  pseudo_data <- list(label=c(array(1, length(pseudo_positive_weight)), array(-1, length(trade_off_weight*pseudo_negative_weight))),
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
