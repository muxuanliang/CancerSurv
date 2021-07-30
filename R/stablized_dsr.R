stablized_dsr=function(dataFrame=NULL, subjectId=NULL, timeToEvent=NULL, measureTime=NULL, covariate=NULL, censoringIndicator=NULL, tau0=6, tradeoff=0.75, method=c("osf", "sm", "pc_glm"), methodFitCensoring="km"){
  if (is.null(dataFrame)){
    dataFrame = data.frame(subjectId=subjectId, timeToEvent=timeToEvent, measureTime=measureTime, covariate=covariate, censoringIndicator=censoringIndicator)
  }

  times=as.numeric(names(table(dataFrame$measureTime)))
  stopifnot(all(dataFrame$timeToEvent>=dataFrame$measureTime))

  dataFrame$indicator=ifelse(dataFrame$timeToEvent-dataFrame$measureTime<=tau0, 1,0)

  if (method != "pc_glm"){
    for(i in 1:length(times)){
      sub_data=dataFrame[dataFrame$measureTime==times[i],]
      sub_data$ipcw=rep(0, nrow(sub_data))
      sub_data$indicator=ifelse(sub_data$timeToEvent-sub_data$measureTime<=tau0, 1,0)
      denomm[i]=sum(sub_data$timeToEvent-sub_data$measureTime<=tau0)/nrow(sub_data)
      censored= 1-sub_data$censoringIndicator
      if (methodFitCensoring=="km"){
        km <- survival::survfit(survival::Surv(sub_data$timeToEvent,censored)~1)
        survest <- stepfun(km$time, c(1, km$surv))
        sub_data$ipcw[sub_data$indicator==1]=survest(sub_data$timeToEvent[sub_data$indicator==1])
        sub_data$ipcw[sub_data$indicator==0]=survest(u+times[i])
      } else if (methodFitCensoring=="cox"){
        coxfit <- survival::coxph(survival::Surv(sub_data$timeToEvent,censored)~sub_data$covariate)
        sub_data$ipcw[sub_data$indicator==1]=predict(coxfit, newdata=sub_data$timeToEvent[sub_data$indicator==1], type="survival")
        sub_data$ipcw[sub_data$indicator==0]=predict(coxfit, newdata=u+times[i], type="survival")
      }
      idx=NA
      for(k in 1:nrow(sub_data)){
        idx[k]=which(sub_data$subjectId[k]==dataFrame$subjectId)[i]
      }
      ipc[idx]=sub_data$ipcw
    }
    dataFrame$ipc=ipc
  }

  if (method=="osf"){
    for(i in 1:length(times)){
      sub_data=dataFrame[dataFrame$measureTime==times[i],]
      idx=which(sub_data$censoringIndicator==1|(sub_data$censoringIndicator==0&sub_data$indicator==0))
      sub_data=sub_data[idx,]
      denomm[i]=sum(ifelse(sub_data$timeToEvent-sub_data$measureTime<=tau0, 1,0)/sub_data$ipc)/sum(1/sub_data$ipc)
    }

    weighted=weight=denom=NA
    for(i in 1:nrow(dataFrame)){
      sub_data=dataFrame[dataFrame$measureTime==dataFrame$measureTime[i],]
      idd=dataFrame$measureTime[i]/6+1
      denom[i]=denomm[idd]
      w1=ifelse((dataFrame$timeToEvent[i]-dataFrame$measureTime[i])<=tau0, 1, 0)/denom[i]
      w2=ifelse((dataFrame$timeToEvent[i]-dataFrame$measureTime[i])>tau0, 1, 0)/(1-denom[i])
      if(dataFrame$indicator[i]==0|dataFrame$censoringIndicator[i]==1){
        weight[i]=(w1-(tradeoff*w2))/dataFrame$ipc[i]
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

    covariate=dataFrame$covariate
    outt=cbind(abs(weighted)[Idx]*pr, A, covariate)
    colnames(outt)=c("outcome", "A", paste0("v", 1:(NCOL(covariate))))
    outt=as.data.frame(outt)
    outt[,2]=as.character(A)
    moPropen <- buildModelObj(model = ~1,solver.method = 'glm',
                              solver.args = list('family'='binomial'),
                              predict.method = 'predict.glm',
                              predict.args = list(type='response'))


    earlRes=try(owl(moPropen = moPropen,
                    data = outt, reward=outt$outcome, txName = "A",
                    regime = ~ v1+v2,
                    lambdas =c(1,5, 10), cvFolds = 2,
                    kernel = "linear",#lambdas =5,cvFolds = 0L,
                    surrogate = 'hinge',kparam = NULL), TRUE)


    coef_osf = regimeCoef(earlRes)
    res = list(coef=stan(coef_osf))
    } else if (method=="sm") {

    } else if (method=="pc_glm") {
      for(i in 1:length(times)){
        sub_data=dataFrame[dataFrame$measureTime==times[i],]
        sub_data$ipcw=rep(0, nrow(sub_data))
        sub_data$ipcw[which(sub_data$indicator==1)]=rms::survest(sub_data$timeToEvent[which(sub_data$indicator==1)])
        sub_data$ipcw[which(sub_data$indicator==0)]=rms::survest(tau0+times[i])
        idx=NA
        for(k in 1:nrow(sub_data)){
          idx[k]=which(sub_data$subjectId[k]==dataFrame$subjectId & dataFrame$measureTime==times[i])
        }
        ipc[idx]=sub_data$ipcw
      }
      dataFrame$ipcc=ipc
    }




  }
















  glm_c=list()
  # for(i in 1:nrow(pc_data)){
  #  sub_data=pc_data[pc_data$meas.time==pc_data$meas.time[i],]
  pc_data$indicator=ifelse((pc_data$time-pc_data$meas.time)<=tau0, 1,0)
  #}
  # pc_data$JJ= relevel(as.factor(pc_data$meas.time) , ref="0")
  for(i in 1:length(times)){
    sub_data=pc_data[pc_data$meas.time==times[i],]
    glm_c[[i]]=glm(indicator~marker_1+marker_2, data=sub_data,weight=1/ipc,family="binomial")
  }


  response=NA

  for(i in 1:nrow(pc_data)){
    #  sub_data=pc_data[pc_data$meas.time==pc_data$meas.time[i],]
    Id=which(pc_data$meas.time[i]==times)
    response[i]= predict(glm_c[[Id]], newdata=pc_data[i,4:5], type="response")
  }
  #coef_glm=summary(lm(response~pc_data$marker_1+pc_data$marker_2))[[4]][,1]
  pc_data$response=response
  weight_glm=NA

  visit_index=pc_data$meas.time/tau0+1


  weight_glm=pc_data$response/ denomm[visit_index]-epsilon*((1-pc_data$response)/(1- denomm[visit_index]))

  # coef_glm=summary(glm(weight_glm~pc_data$marker_1+pc_data$marker_2))[[4]][,1]
  # weighted=weight=NA
  pc_data$weight_glm=weight_glm

  model=paste0("geeglm(ifelse(weight_glm>0, 1, 0)~", paste(colnames(covariate), collapse=" + "),", id=sub.id, corstr='independence', data=pc_data,family='binomial')" )

  coef_sm=coef(eval(parse(text=model)))

  #coef(glm(ifelse(weight_glm>0, 1, 0)~marker_1+marker_2, data=pc_data, family="binomial")


  model=paste0("glm(indicator~", paste(colnames(covariate), collapse=" + "),",data=pc_data,weight=1/ipcc,family='binomial')" )
  #model=paste0("geeglm(indicator~", paste(colnames(covariate), collapse=" + "),", id=sub.id,data=pc_data,family='binomial')" )
  pc_glm=(eval(parse(text=model)))

  prob=predict(pc_glm, type="response", newdata=pc_data)
  pc_data$risk_2=prob

  Cutoff=seq(1:100)*0.01
  Phi=Tpf=Fpf=NA
  for(kkl in 1:length(Cutoff)){
    cutoff=Cutoff[kkl]
    cutoff=rep(cutoff, 5)
    rule=NA
    for(i in 1:nrow(pc_data)){
      rule[i]=ifelse(pc_data$risk_2[i]>cutoff[visit_index], 1, -1)
    }

    pc_data$rule=rule
    phi=tpf=fpf=NA
    for(i in 1:length(times)){
      sub_data=pc_data[pc_data$meas.time==times[i],]

      #denom=sum((sub_data$time-sub_data$meas.time)<=tau0)/nrow(sub_data)
      idx=which(sub_data$status==1|(sub_data$status==0&sub_data$indicator==0))
      sub_data=sub_data[idx,]
      sub1=sub_data[(sub_data$time-sub_data$meas.time)<=tau0,]

      sub2=sub_data[(sub_data$time-sub_data$meas.time)>tau0,]

      phi[i]=sum(sub1$rule==1)/nrow(sub1)-epsilon*sum(sub2$rule==1)/nrow(sub2)
      phi[i]=sum( ifelse(sub1$rule==1, 1, 0)/sub1$ipc)/sum(1/sub1$ipc)-epsilon*(sum( ifelse(sub2$rule==1, 1, 0)/sub2$ipc)/sum(1/sub2$ipc))
      tpf[i]=sum( ifelse(sub1$rule==1, 1, 0)/sub1$ipc)/sum(1/sub1$ipc)
      fpf[i]=sum( ifelse(sub2$rule==1, 1, 0)/sub2$ipc)/sum(1/sub2$ipc)
    }

    Phi[kkl]=mean(phi)
    Tpf[kkl]=mean(tpf)
    Fpf[kkl]=mean(fpf)
  }
  idd=which(Phi==max(Phi))
  cutoff=Cutoff[idd]
  PC_GLM=c(coef(pc_glm)[1]-log(cutoff/(1-cutoff)), coef(pc_glm)[-1])



  result= list(stan(coef_sm),stan(coef_osf), stan(PC_GLM) )
  return(result)
}
