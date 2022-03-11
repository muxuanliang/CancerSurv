cancerSurvInterval <-function(biopsyTable.subjectId, biopsyTable.measureTime, biopsyTable.eventTime,
                           PSATable.covariate, PSATable.subjectId, PSATable.measureTime,
                           PSATable.measureTimeDiscrete=NULL, tau0=0.5, time_varying=FALSE,
                           timePoints=NULL, xi=10, lambdaSeq=NULL, nlambda=100, n.fold=2, method='logit', standardize = FALSE){
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
                            is.positive=abs(biopsyTable.measureTime-biopsyTable.eventTime)<0.1,
                            discard=biopsyTable.measureTime > biopsyTable.eventTime)

  covariate.train <- NULL

  # get biopsyTable
  biopsyTable <- data.frame(subjectId=biopsyTable.subjectId, measureTime=biopsyTable.measureTime,
                            is.positive=abs(biopsyTable.measureTime-biopsyTable.eventTime)<0.1,
                            discard=biopsyTable.measureTime > biopsyTable.eventTime,
                            decision=NA)

  pseudo_covariate <- pseudo_time <- pseudo_positive_weight <- pseudo_negative_weight <- NULL
  for (time in timePoints){
    tmp_biopsy_table <- biopsyTable
    pidList <- unique(tmp_biopsy_table$subjectId)
    for (pid in pidList){
      boolVector <- (dataFrame$measureTime[dataFrame$subjectId==pid]-time)<=tau0 & (dataFrame$measureTime[dataFrame$subjectId==pid]-time)>0
      if(!any(boolVector)){
        tmp_biopsy_table$discard[tmp_biopsy_table$subjectId==pid] <- TRUE
      }
    }
    tmp_biopsy_table <- tmp_biopsy_table[!tmp_biopsy_table$discard,]

    for (pid in pidList){
      boolVector <- (dataFrame$measureTime[dataFrame$subjectId==pid]-time)<=tau0 & (dataFrame$measureTime[dataFrame$subjectId==pid]-time)>0
      if(any(boolVector)){
        select.index <- min(which(boolVector))
        selected.covariate <- PSATable.covariate[dataFrame$subjectId==pid,,drop=FALSE][select.index,,drop=FALSE]
      } else {
        next
      }

    select_biopsy_table <- tmp_biopsy_table[tmp_biopsy_table$subjectId==pid,]

    tmp_positive_weight <- ks_sum(select_biopsy_table$measureTime, as.numeric(select_biopsy_table$is.positive), time+tau0)/ks_sum(tmp_biopsy_table$measureTime, 1+0*as.numeric(tmp_biopsy_table$is.positive), time+tau0)
    denom_positive <- ks(tmp_biopsy_table$measureTime, as.numeric(tmp_biopsy_table$is.positive), time+tau0)-ks(tmp_biopsy_table$measureTime, as.numeric(tmp_biopsy_table$is.positive), time)
    positive_weight <-tmp_positive_weight/denom_positive

    tmp_negative_weight <- ks_sum(select_biopsy_table$measureTime, as.numeric(!select_biopsy_table$is.positive), time+tau0)/ks_sum(tmp_biopsy_table$measureTime, 1+0*as.numeric(tmp_biopsy_table$is.positive), time+tau0)
    denom_negative <- 1-ks(tmp_biopsy_table$measureTime, as.numeric(tmp_biopsy_table$is.positive), time+tau0)
    tmp_negative_weight2 <- ks_sum(select_biopsy_table$measureTime, as.numeric(select_biopsy_table$is.positive), time)/ks_sum(tmp_biopsy_table$measureTime, 1+0*as.numeric(tmp_biopsy_table$is.positive), time)
    negative_weight <- tmp_negative_weight/denom_negative+tmp_negative_weight2/denom_positive

    pseudo_covariate <- rbind(pseudo_covariate,selected.covariate)
    pseudo_time <- c(pseudo_time, time)
    pseudo_positive_weight <- c(pseudo_positive_weight, positive_weight)
    pseudo_negative_weight <- c(pseudo_negative_weight, negative_weight)
    }
  }

  # pseudo_data
  pseudo_data <- list(label=rep(c(1,-1), each=NROW(pseudo_covariate)),
                      covariate=rbind(pseudo_covariate,pseudo_covariate),
                      time=pseudo_time,
                      weight=c(pseudo_positive_weight, xi*c(pseudo_negative_weight)))

  if (method=='logit'){
    fit <- glmnet::cv.glmnet(x=cbind(1,pseudo_data$covariate), y=(pseudo_data$label>0), weights=pseudo_data$weight, family='binomial', intercept = FALSE, standardize = FALSE)
    return(list(fit=c(fit$glmnet.fit$a0[fit$lambda==fit$lambda.min],fit$glmnet.fit$beta[-1,fit$lambda==fit$lambda.min]), lambdaSeq = fit$lambda, lambda_opt=fit$lambda.min, cvm=fit$cvm, opt=list(mean=covariate.mean, sd=covariate.sd)))
  }

  if (is.null(lambdaSeq)) {
    kkt <- apply(pseudo_data$covariate, 2, function(t){sum(t*pseudo_data$weight*pseudo_data$label)})
    lambdaInit <- max(abs(kkt))
    lambdaSeq <- lambdaInit * (10/11)^(0:(nlambda-1))
  }

  if(length(lambdaSeq)>1){
    pidList <- unique(pseudo_data$pid)
    fold.id <- sample(seq_len(5), length(pidList), replace = TRUE)
    cvm <- array(0, c(n.fold, length(lambdaSeq)))
    for (index in 1:n.fold){
      training.id <- pidList[fold.id==index]
      validation.id <- pidList[fold.id!=index]
      fit <- lapply(lambdaSeq, function(t){
        bmrm::svmLP(x=cbind(1,pseudo_data$covariate[pseudo_data$pid %in% training.id,]), y=(pseudo_data$label[pseudo_data$pid %in% training.id]>0), LAMBDA = t, loss.weights = pseudo_data$weight[pseudo_data$pid %in% training.id])
      })

      validate <- sapply(fit, function(t){
        predict(t, x = cbind(1,pseudo_data$covariate[pseudo_data$pid %in% validation.id,]))
      })
      cvm[index,] <- apply(validate, 2, function(t){
        score <- mean(pseudo_data$weight[pseudo_data$pid %in% validation.id]*loss(pseudo_data$label[pseudo_data$pid %in% validation.id] * (t-0.5), loss_type='zeroOne'))
        score
      })
    }
    avecvm <- apply(cvm,2,mean)
    lambda_opt <- lambdaSeq[which.min(avecvm)]
  } else {
    lambda_opt <- lambdaSeq
  }

  fit <- bmrm::svmLP(x=cbind(1,pseudo_data$covariate), y=(pseudo_data$label>0), LAMBDA = lambda_opt, loss.weights = pseudo_data$weight)


  list(fit=fit, lambdaSeq=lambdaSeq, lambda_opt=lambda_opt, cvm=cvm, opt=list(mean=covariate.mean, sd=covariate.sd))
}
