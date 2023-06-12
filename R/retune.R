retune <-function(fit = NULL,biopsyTable.subjectId, biopsyTable.measureTime, biopsyTable.eventTime,
                  biopsyTable.primary_gleason, biopsyTable.secondary_gleason,
                  PSATable.covariate, PSATable.subjectId, PSATable.measureTime, PSATable.measureTimeDiscrete,
                  tau0=0.5, time_varying=FALSE, timePoints=NULL, opt=NULL, xi=0.5, prevalence=NULL){


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
  cutoff_seq <- seq(-min(unique(PSATable.covariate %*% fit$coef[-1])),
                    -max(unique(PSATable.covariate %*% fit$coef[-1])),
                    length.out=100)


  # try each cutoff
  value <- sapply(cutoff_seq, function(t){
    res_tmp <- get_tpr_tnr_proposed(fit = list(coef=c(t, fit$coef[-1])),
                                    biopsyTable.subjectId, biopsyTable.measureTime, biopsyTable.eventTime,
                                    biopsyTable.primary_gleason, biopsyTable.secondary_gleason,
                                    PSATable.covariate, PSATable.subjectId, PSATable.measureTime, PSATable.measureTimeDiscrete,
                                    tau0=tau0, time_varying=FALSE, timePoints=timePoints, opt=opt)

    if (is.null(prevalence)){
      value <- mean(res_tmp$TPR + xi * (1-res_tmp$prevalence)/res_tmp$prevalence * res_tmp$TNR)
    } else {
      value <- mean(res_tmp$TPR + xi * (1-prevalence)/prevalence * res_tmp$TNR)
    }
    value
  })

  list(coef=c(cutoff_seq[min(which.max(value))], fit$coef[-1]), value_opt=value[which.max(value)])
}
