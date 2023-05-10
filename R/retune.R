retune <-function(fit = NULL,biopsyTable.subjectId, biopsyTable.measureTime, biopsyTable.eventTime,
                  biopsyTable.primary_gleason, biopsyTable.secondary_gleason,
                  PSATable.covariate, PSATable.subjectId, PSATable.measureTime, PSATable.measureTimeDiscrete,
                  tau0=0.5, time_varying=FALSE, timePoints=NULL, opt=NULL, xi=0.5){


  # get candidate
  numberCovariate <- NCOL(PSATable.covariate)
  if (is.null(colnames(PSATable.covariate))) {
    colnames(PSATable.covariate) <- paste0("v", 1:(NCOL(PSATable.covariate)))
  }
  dataFrame <-
    data.frame(
      covariate = PSATable.covariate,
      subjectId = PSATable.subjectId,
      measureTime = PSATable.measureTime,
      measureTimeDiscrete = PSATable.measureTimeDiscrete
    )
  if(!is.null(opt)){
    PSATable.covariate <- t(apply(PSATable.covariate,1,function(t){(t-opt$mean)/opt$sd}))
  }
  cutoff_seq <- sort(unique(PSATable.covariate %*% fit$coef[-1]))


  # try each cutoff
  value <- array(NA, dim=length(cutoff_seq))
  for (cutoff in cutoff_seq){
    res_tmp <- get_tpr_tnr_proposed(fit = list(coef=c(cutoff, fit$coef[-1])),
      biopsyTable.subjectId, biopsyTable.measureTime, biopsyTable.eventTime,
      biopsyTable.primary_gleason, biopsyTable.secondary_gleason,
      PSATable.covariate, PSATable.subjectId, PSATable.measureTime, PSATable.measureTimeDiscrete,
      tau0=tau0, time_varying=FALSE, timePoints=timePoints, opt=opt)

    value[cutoff==cutoff_seq] <- mean(res_tmp$TPR + xi * (1-res_tmp$prevalence)/res_tmp$prevalence * res_tmp$TNR)
  }
  list(coef=c(cutoff_seq[min(which.max(value))], fit$coef[-1]), value_opt=value[which.max(value)])
}
