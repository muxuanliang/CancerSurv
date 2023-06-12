cv.cancerSurvInterval <-function(biopsyTable.subjectId, biopsyTable.measureTime, biopsyTable.eventTime,
                              biopsyTable.primary_gleason, biopsyTable.secondary_gleason,
                              PSATable.covariate, PSATable.subjectId, PSATable.measureTime,
                              PSATable.measureTimeDiscrete=NULL, tau0=0.5, time_varying=FALSE,
                              timePoints=NULL, xi=10, lambdaSeq=NULL, nlambda=100, n.fold=2,
                              standardize = FALSE, hopt = 0.2*sd(biopsyTable$measureTime) *
                                (length(unique(PSATable.subjectId)))^(-1/5),
                              number_shift_fractions = 20, prevalence=NULL){

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

  # cv select shift_fraction
  shift_fraction_seq <- seq(0, 1, length.out=number_shift_fractions)
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
                                   shift_fraction = shift_fraction_tmp, prevalence=prevalence)
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
                                   shift_fraction = shift_fraction_tmp, prevalence=prevalence)

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
    val <- rbind(val, c(res=(res_fit1_tmp$TPR+xi*(1-prevalence)/prevalence*res_fit1_tmp$TNR+
                               res_fit2_tmp$TPR+xi*(1-prevalence)/prevalence*res_fit2_tmp$TNR)/2,
                        shift_fraction=shift_fraction_tmp))
  }
  val_sum <- apply(val[,1:length(timePoints)],1,mean)
  shift_fraction_select <- shift_fraction_seq[min(which.max(val_sum))]


  fit_select <- cancerSurvInterval(biopsyTable.subjectId,
                                 biopsyTable.measureTime,
                                 biopsyTable.eventTime,
                                 biopsyTable.primary_gleason,
                                 biopsyTable.secondary_gleason,
                                 PSATable.covariate,
                                 PSATable.subjectId,
                                 PSATable.measureTime,
                                 PSATable.measureTimeDiscrete,
                                 tau0=tau0, time_varying=time_varying,
                                 timePoints=timePoints, xi=xi, lambdaSeq=lambdaSeq,
                                 nlambda=nlambda, n.fold=n.fold,
                                 standardize = standardize, hopt = hopt,
                                 shift_fraction = shift_fraction_select, prevalence=prevalence)

  return(list(fit=fit_select$fit,opt=fit_select$opt, prevalence=prevalence, shift_fraction=shift_fraction_select, cv.value=val_sum))
}
