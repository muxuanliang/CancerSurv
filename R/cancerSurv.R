cancerSurv=function(subjectId=NULL, timeToEvent=NULL, measureTime=NULL, measureTimeDiscrete=NULL, covariate=NULL, censoringIndicator=NULL, tau0=6, tradeoff=0.75, trade_off_type=c("standard", "cbr"), method=c("osf", "sm", "pc_glm", 'jama'), methodFitCensoring="km", time_varying = FALSE, coef_input = c(0.489,0.077, -0.968, 0,0.27, 0.476,-0.654,-1.715), prevalence=NULL){
  if (!time_varying){
  numberCovariate <- NCOL(covariate)
  if (is.null(colnames(covariate))){
    colnames(covariate) <- paste0("v", 1:(NCOL(covariate)))
  }
  dataFrame = data.frame(covariate, subjectId=subjectId, timeToEvent=timeToEvent, measureTime=measureTime, measureTimeDiscrete=measureTimeDiscrete, censoringIndicator=censoringIndicator)

  times=as.numeric(names(table(dataFrame$measureTimeDiscrete)))
  stopifnot(all(dataFrame$timeToEvent>=dataFrame$measureTime))

  all_survival <- data.frame(subjectId=dataFrame$subjectId,
                                  timeToEvent=dataFrame$timeToEvent,
                                  censoringIndicator=dataFrame$censoringIndicator)
  all_survival <- dplyr::distinct(all_survival)
  km_all <- survival::survfit(survival::Surv(all_survival$timeToEvent,all_survival$censoringIndicator)~1)
  survest_all <- stepfun(km_all$time, c(1, km_all$surv))

  dataFrame$indicator=ifelse(dataFrame$timeToEvent-dataFrame$measureTime<=tau0, 1,0)
  dataFrame$ipc <- NULL
  dataFrame$trade_off_weight <- trade_off_weight <-  NULL
  for(i in 1:length(times)){
    sub_data=dataFrame[dataFrame$measureTimeDiscrete==times[i],]
    sub_data$ipcw=rep(0, nrow(sub_data))
    sub_data$indicator=ifelse(sub_data$timeToEvent-sub_data$measureTime<=tau0, 1,0)
    sub_data_survival <- data.frame(subjectId=dataFrame$subjectId[dataFrame$measureTimeDiscrete>=times[i]],
                                    timeToEvent=dataFrame$timeToEvent[dataFrame$measureTimeDiscrete>=times[i]],
                                    censoringIndicator=dataFrame$censoringIndicator[dataFrame$measureTimeDiscrete>=times[i]])
    sub_data_survival <- dplyr::distinct(sub_data_survival)
    censored= 1-sub_data_survival$censoringIndicator
    if (methodFitCensoring=="km"){
      km <- survival::survfit(survival::Surv(sub_data_survival$timeToEvent,censored)~1)
      survest <- stepfun(km$time, c(1, km$surv))
      sub_data$ipcw[sub_data$indicator==1]=survest(sub_data$timeToEvent[sub_data$indicator==1])
      sub_data$ipcw[sub_data$indicator==0]=survest(tau0+sub_data$measureTime[sub_data$indicator==0])
    } else if (methodFitCensoring=="cox"){
      model=paste0("coxfit <- survival::coxph(survival::Surv(sub_data$timeToEvent,censored)~", paste(colnames(covariate), collapse=" + "),", data=sub_data)" )
      eval(parse(text=model))
      sub_data$ipcw[sub_data$indicator==1]=predict(coxfit, type="survival")[sub_data$indicator==1]
      sub_data$ipcw[sub_data$indicator==0]=predict(coxfit, newdata=data.frame(sub_data[sub_data$indicator==0, 1:numberCovariate],timeToEvent=rep(tau0+times[i], times=sum(sub_data$indicator==0))), type="survival")
    }
    if (is.null(prevalence)){
      prevalence <- max((survest_all(times[i])-survest_all(tau0+times[i]))/survest_all(times[i]),0.01)
      dataFrame$trade_off_weight[dataFrame$measureTimeDiscrete==times[i]] <- (1-prevalence)/prevalence
      trade_off_weight[i] <- (1-prevalence)/prevalence
    }
    dataFrame$trade_off_weight[dataFrame$measureTimeDiscrete==times[i]] <- (1-prevalence[i])/prevalence[i]
    trade_off_weight[i] <- (1-prevalence[i])/prevalence[i]
    for (j in 1:NROW(sub_data)){
      dataFrame$ipc[(dataFrame$subjectId==sub_data$subjectId[j])&(dataFrame$measureTimeDiscrete==times[i])] <- sub_data$ipcw[j]
    }
  }

  denomm <- NULL
  for(i in 1:length(times)){
    sub_data=dataFrame[dataFrame$measureTimeDiscrete==times[i],]
    idx=which(sub_data$censoringIndicator==1|(sub_data$censoringIndicator==0&sub_data$indicator==0))
    sub_data=sub_data[idx,]
    denomm[i]=sum(ifelse(sub_data$timeToEvent-sub_data$measureTime<=tau0, 1,0)/sub_data$ipc)/sum(1/sub_data$ipc)
  }

  if (trade_off_type=="standard"){
    dataFrame$trade_off_weight <- 1
    trade_off_weight <- 0*trade_off_weight+1
  }

  if (method=="osf"){
    weighted=weight=denom=NA
    for(i in 1:nrow(dataFrame)){
      sub_data=dataFrame[dataFrame$measureTimeDiscrete==dataFrame$measureTimeDiscrete[i],]
      idd=dataFrame$measureTimeDiscrete[i]/tau0 - min(dataFrame$measureTimeDiscrete)/tau0 + 1
      denom[i]=denomm[idd]
      w1=ifelse((dataFrame$timeToEvent[i]-dataFrame$measureTime[i])<=tau0, 1, 0)/denom[i]
      w2=ifelse((dataFrame$timeToEvent[i]-dataFrame$measureTime[i])>tau0, 1, 0)/(1-denom[i])
      if(dataFrame$indicator[i]==0|dataFrame$censoringIndicator[i]==1){
        weight[i]=(w1-(dataFrame$trade_off_weight[i]*tradeoff*w2))/dataFrame$ipc[i]
      } else {
        weight[i]=0
      }
      weighted[i]=weight[i]/sum(1/sub_data$ipc)
    }

    Idx=which(dataFrame$censoringIndicator==1|(dataFrame$censoringIndicator==0&dataFrame$indicator==0))
    dataFrame=dataFrame[Idx,]
    dataFrame$weight=weight[Idx]
    A=sign(weighted)[Idx]
    pr=rep(0, length(A))
    pr[A==1]=sum(A==1)/length(A)
    pr[A==-1]=sum(A==-1)/length(A)
    dataFrame$pr=pr

    covariate=dataFrame[1:numberCovariate]
    outt=cbind(abs(weighted)[Idx]*pr, A, covariate)
    colnames(outt)=c("outcome", "A", paste0("v", 1:(NCOL(covariate))))
    outt=as.data.frame(outt)
    outt[,2]=as.character(A)
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

    res = list(coef=stan(coef_osf))
  } else if (method=="sm") {
    covariate=covariate[which(dataFrame$censoringIndicator==1|(dataFrame$censoringIndicator==0&dataFrame$indicator==0)),]
    dataFrame=dataFrame[which(dataFrame$censoringIndicator==1|(dataFrame$censoringIndicator==0&dataFrame$indicator==0)),]
    glm_c=list()
    for(i in 1:length(times)){
      sub_data = dataFrame[dataFrame$measureTimeDiscrete==times[i],]
      model=paste0("glm(indicator~", paste(colnames(covariate), collapse=" + "),",data=sub_data,weight=1/ipc,family='binomial')" )
      glm_c[[i]]=(eval(parse(text=model)))
    }

    response=NA
    for(i in 1:nrow(dataFrame)){
      Id=which(dataFrame$measureTimeDiscrete[i]==times)
      response[i]= predict(glm_c[[Id]], newdata=as.data.frame(t(as.data.frame(covariate[i,]))), type="response")
    }
    dataFrame$response=response

    weight_glm = NULL
    visit_index = dataFrame$measureTimeDiscrete/tau0-min(dataFrame$measureTimeDiscrete)/tau0 + 1
    weight_glm = dataFrame$response/ denomm[visit_index]-trade_off_weight[visit_index]*tradeoff*((1-dataFrame$response)/(1- denomm[visit_index]))
    dataFrame$weight_glm = weight_glm
    model=paste0("geepack::geeglm(ifelse(weight_glm>0, 1, 0)~", paste(colnames(covariate), collapse=" + "),", id=subjectId, corstr='independence', data=dataFrame,family='binomial')" )
    coef_sm=coef(eval(parse(text=model)))

    res = list(coef=stan(coef_sm))
  } else if (method=="pc_glm") {
    id=(unique(dataFrame$subjectId))
    idx=NA
    for(i in 1:length(id)){
      idx[i]=which(dataFrame$subjectId==id[i])[1]
    }
    sub_data=dataFrame[idx,]
    sub_data$ipcw=rep(0, nrow(sub_data))
    sub_data$indicator=ifelse(sub_data$timeToEvent-sub_data$measureTime<=tau0, 1,0)

    censored= 1-sub_data$censoringIndicator
    km <- survival::survfit(survival::Surv(sub_data$time,censored)~1)
    survest <- stepfun(km$time, c(1, km$surv))

    dataFrame$ipcc <- rep(0, times=NROW(dataFrame))
    for(i in 1:length(times)){
      sub_data=dataFrame[dataFrame$measureTimeDiscrete==times[i],]
      sub_data$ipcw=rep(0, nrow(sub_data))
      sub_data$ipcw[which(sub_data$indicator==1)]=survest(sub_data$timeToEvent[which(sub_data$indicator==1)])
      sub_data$ipcw[which(sub_data$indicator==0)]=survest(tau0+sub_data$measureTime[which(sub_data$indicator==0)])
      idx=NA
      for(k in 1:nrow(sub_data)){
        idx[k]=which(sub_data$subjectId[k]==dataFrame$subjectId & dataFrame$measureTimeDiscrete==times[i])
      }
      dataFrame$ipcc[idx]=sub_data$ipcw
    }

    covariate=covariate[which(dataFrame$censoringIndicator==1|(dataFrame$censoringIndicator==0&dataFrame$indicator==0)),]
    dataFrame=dataFrame[which(dataFrame$censoringIndicator==1|(dataFrame$censoringIndicator==0&dataFrame$indicator==0)),]

    model=paste0("glm(indicator~", paste(colnames(covariate), collapse=" + "),",data=dataFrame,weight=1/ipcc,family='binomial')" )
    pc_glm=(eval(parse(text=model)))

    prob=predict(pc_glm, type="response", newdata=data.frame(covariate))
    dataFrame$risk_score=prob

    Cutoff=seq(1:100)*0.01
    Phi=Tpf=Fpf=NA
    visit_index=dataFrame$measureTimeDiscrete/tau0 - min(dataFrame$measureTimeDiscrete)/tau0 + 1
    for(kkl in 1:length(Cutoff)){
      cutoff=Cutoff[kkl]
      cutoff=rep(cutoff, length(times))
      rule=NA
      for(i in 1:nrow(dataFrame)){
        rule[i]=ifelse(dataFrame$risk_score[i]>cutoff[visit_index[i]], 1, -1)
      }

      dataFrame$rule=rule
      phi=tpf=fpf=NA
      for(i in 1:length(times)){
        sub_data=dataFrame[dataFrame$measureTimeDiscrete==times[i],]
        idx=which(sub_data$censoringIndicator==1|(sub_data$censoringIndicator==0&sub_data$indicator==0))
        sub_data=sub_data[idx,]
        sub1=sub_data[(sub_data$timeToEvent-sub_data$measureTime)<=tau0,]
        sub2=sub_data[(sub_data$timeToEvent-sub_data$measureTime)>tau0,]

        phi[i]=sum( ifelse(sub1$rule==1, 1, 0)/sub1$ipc)/sum(1/sub1$ipc)-trade_off_weight[i]*tradeoff*(sum( ifelse(sub2$rule==1, 1, 0)/sub2$ipc)/sum(1/sub2$ipc))
        tpf[i]=sum( ifelse(sub1$rule==1, 1, 0)/sub1$ipc)/sum(1/sub1$ipc)
        fpf[i]=sum( ifelse(sub2$rule==1, 1, 0)/sub2$ipc)/sum(1/sub2$ipc)
      }

      Phi[kkl]=mean(phi)
      Tpf[kkl]=mean(tpf)
      Fpf[kkl]=mean(fpf)
    }
    idd=which(Phi==max(Phi))
    cutoff=min(max(Cutoff[idd]),0.9)
    PC_GLM=c(coef(pc_glm)[1]-log(cutoff/(1-cutoff)), coef(pc_glm)[-1])

    res = list(coef=stan(PC_GLM))
  } else if (method=="jama") {
    coef <- coef_input

    covariate=covariate[which(dataFrame$censoringIndicator==1|(dataFrame$censoringIndicator==0&dataFrame$indicator==0)),]
    dataFrame=dataFrame[which(dataFrame$censoringIndicator==1|(dataFrame$censoringIndicator==0&dataFrame$indicator==0)),]

    dataFrame$risk_score=covariate %*%coef

    Cutoff=unique(dataFrame$risk_score)
    Phi=Tpf=Fpf=NA
    visit_index=dataFrame$measureTimeDiscrete/tau0 - min(dataFrame$measureTimeDiscrete)/tau0 + 1
    for(kkl in 1:length(Cutoff)){
      cutoff=Cutoff[kkl]
      cutoff=rep(cutoff, length(times))
      rule=NA
      for(i in 1:nrow(dataFrame)){
        rule[i]=ifelse(dataFrame$risk_score[i]>cutoff[visit_index[i]], 1, -1)
      }

      dataFrame$rule=rule
      phi=tpf=fpf=NA
      for(i in 1:length(times)){
        sub_data=dataFrame[dataFrame$measureTimeDiscrete==times[i],]
        idx=which(sub_data$censoringIndicator==1|(sub_data$censoringIndicator==0&sub_data$indicator==0))
        sub_data=sub_data[idx,]
        sub1=sub_data[(sub_data$timeToEvent-sub_data$measureTime)<=tau0,]
        sub2=sub_data[(sub_data$timeToEvent-sub_data$measureTime)>tau0,]

        phi[i]=sum( ifelse(sub1$rule==1, 1, 0)/sub1$ipc)/sum(1/sub1$ipc)-trade_off_weight[i]*tradeoff*(sum( ifelse(sub2$rule==1, 1, 0)/sub2$ipc)/sum(1/sub2$ipc))
        tpf[i]=sum( ifelse(sub1$rule==1, 1, 0)/sub1$ipc)/sum(1/sub1$ipc)
        fpf[i]=sum( ifelse(sub2$rule==1, 1, 0)/sub2$ipc)/sum(1/sub2$ipc)
      }

      Phi[kkl]=mean(phi)
      Tpf[kkl]=mean(tpf)
      Fpf[kkl]=mean(fpf)
    }
    idd=which(Phi==max(Phi))
    cutoff=max(Cutoff[idd])
    JAMA=c(-cutoff, coef)

    res = list(coef=stan(JAMA))
  }
  return(res)
  } else {
    timePoint <- unique(measureTimeDiscrete)
    res <- NULL
    index<-1
    for (time_point in timePoint){
      res[[index]] <- cancerSurv(subjectId=subjectId[measureTimeDiscrete==time_point], timeToEvent=timeToEvent[measureTimeDiscrete==time_point],
                 measureTime=measureTime[measureTimeDiscrete==time_point], measureTimeDiscrete=measureTimeDiscrete[measureTimeDiscrete==time_point],
                 covariate=covariate[measureTimeDiscrete==time_point,], censoringIndicator=censoringIndicator[measureTimeDiscrete==time_point],
                 tau0=tau0, tradeoff=tradeoff, trade_off_type = trade_off_type, method=method, methodFitCensoring=methodFitCensoring, time_varying = FALSE)
      index <- index+1
    }
    return(res)
  }
}
