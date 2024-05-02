

brms_mediation <- function(model.m, model.y, model.my = NULL,
                           mediator = "med.name(s)", treat = "treat.name", outcome = NULL,
                           control.value = 0, treat.value = 1, 
                           dpar.outcome = "mu", dpar.mediator = NULL, mediator0 = NULL,
                           type = c("response","lp","HR","logHR"))
{
  if (is.null(mediator0)) mediator0 = mediator
  value = c(control.value, treat.value)
  type = type[1]
  if (!is.null(model.my)) 
  {
    model.m = model.my
    model.y = model.my
    if (is.null(outcome)) stop("'outcome' should be given")
  }
  
  dat.new = model.y$data
  outcome.lin = array(NA, dim = c(2, ndraws(model.y), nrow(dat.new)) )
  for (i in 1:2){
    dat.new[,treat] = value[i]
    outcome.lin[i,,] = posterior_linpred(model.y, newdata=dat.new, resp=outcome, dpar=dpar.outcome)
  }
  dat.new = model.m$data
  mediator.pred = array(NA, dim = c(2, length(mediator), ndraws(model.m), nrow(dat.new)) )
  for (i in 1:2){
    dat.new[,treat] = value[i]
    for (m in 1:length(mediator))
      mediator.pred[i,m,,] = posterior_epred(model.m, newdata=dat.new, resp=mediator0[m], dpar=dpar.mediator)
  }
  
  if (is.null(model.my))
  {
    dat = model.y$data
    X = standata(model.y)[["X"]]
    if (!is.null(dpar.outcome))
      if (dpar.outcome!="mu") 
        X = standata(model.y)[[paste("X_", dpar.outcome, sep="")]]
    vars = colnames(X)
    var = vars[grepl(treat,vars)]
    for (m in 1:length(mediator))
      var = c(var, vars[grepl(mediator[m],vars)])
    var = unique(var)
    b = as.matrix(model.y)[, c(paste("b_", var, sep=""))] 
    if (!is.null(dpar.outcome))
      if (dpar.outcome!="mu")
        b = as.matrix(model.y)[, c(paste("b_", dpar.outcome, "_", var, sep=""))] 
    colnames(b) = var
  }
  if (!is.null(model.my))
  {
    dat = model.y$data
    X = standata(model.y)[[paste("X_", outcome, sep="")]]
    if (!is.null(dpar.outcome))
      if (dpar.outcome!="mu") 
        X = standata(model.y)[[paste("X_", dpar.outcome, "_", outcome, sep="")]]
    vars = colnames(X)
    var = vars[grepl(treat,vars)]
    for (m in 1:length(mediator))
      var = c(var, vars[grepl(mediator[m],vars)])
    var = unique(var)
    b = as.matrix(model.y)[, c(paste("b_", outcome, "_", var, sep=""))] 
    if (!is.null(dpar.outcome))
      if (dpar.outcome!="mu")
        b = as.matrix(model.y)[, c(paste("b_", dpar.outcome, "_", outcome, "_", var, sep=""))] 
    colnames(b) = var
  }
  
  if (is.factor(dat[,treat])) {
    nam = paste(treat, levels(dat[,treat])[-1], sep="")
    x.treat = array(0, dim=c(2, nrow(dat), nlevels(dat[,treat])-1), dimnames=list(NULL,NULL,nam))
    for (i in 1:2) 
      x.treat[i,,levels(dat[,treat])[-1]==value[i]] = 1
  } else   {
    x.treat = array(0, dim=c(2, nrow(dat), 1), dimnames=list(NULL,NULL,treat))
    for (i in 1:2) x.treat[i,,] = value[i]
  }
  
  
  dat.new = dat
  for (j in 1:ncol(dat.new))  { 
    if (is.factor(dat.new[,j])) dat.new[,j] = levels(dat.new[,j])[1]
    else dat.new[,j] = 0
  }
  intercept = posterior_linpred(model.y, newdata=dat.new, resp=outcome, dpar=dpar.outcome)
  xb = outcome.lin
  cov = unique(unlist(strsplit(var,split=":",fixed=F)))
  inc = rep(TRUE, length(cov))
  for (k in 1:length(cov))  {
    if (grepl(treat,cov[k])) inc[k] = FALSE
    for (m in 1:length(mediator))
      if (grepl(mediator[m],cov[k])) inc[k] =FALSE
  }
  cov = cov[inc]
  for (i in 1:2)  {
    dat.new[, treat] = value[i]
    dat.new[, mediator] = dat[, mediator]
    if (length(cov) > 0) dat.new[,cov] = dat[,cov]
    xb[i,,] = posterior_linpred(model.y, newdata=dat.new, resp=outcome, dpar=dpar.outcome) - intercept
  }
  
  
  outcome.linpred = tm.lin = array(NA, dim = c(2, 2, ndraws(model.y), nrow(dat)))
  for (i in 1:2)
    for (j in 1:2)
    {
      tm.lin[i,j,,] = 0
      for (k in 1:length(var))
      {
        xx = strsplit(var[k], split=":", fixed=TRUE)[[1]]
        if (length(xx)==1)
        {
          if (grepl(treat,var[k]))
            tm.lin[i,j,,] = tm.lin[i,j,,] + b[,k,drop=F] %*% x.treat[i,,colnames(b)[k]]
          for (m in 1:length(mediator))
            if (grepl(mediator[m],var[k]))
              tm.lin[i,j,,] = tm.lin[i,j,,] + mediator.pred[j,m,,] * b[,k]
        }
        if (length(xx)==2)
        {
          for (m in 1:length(mediator))
            if (grepl(treat,var[k]) & grepl(mediator[m],var[k])) 
            {
              for (h in 1:dim(x.treat)[3]) 
                if(grepl(dimnames(x.treat)[[3]][h],var[k])) kk = h
              tm.lin[i,j,,] = tm.lin[i,j,,] + (x.treat[i,,colnames(b)[kk]] * mediator.pred[j,m,,]) * b[,k]
            }
          if (grepl(treat,var[k])) 
          {
            exc = NULL
            for (m in 1:length(mediator)) exc = c(exc,!grepl(mediator[m],var[k]))
            if (all(exc))
            {
              x.cov = X[,xx[!grepl(treat,xx)]]
              for (h in 1:dim(x.treat)[3]) 
                if(grepl(dimnames(x.treat)[[3]][h],var[k])) kk = h
              tm.lin[i,j,,] = tm.lin[i,j,,] + b[,k,drop=F] %*% (x.treat[i,,colnames(b)[kk]] * x.cov) 
            }
          }
        }
        if (length(xx) > 2) (stop("not allow three-way interactions"))
      }
      outcome.linpred[i,j,,] = outcome.lin[i,,] - xb[i,,] + tm.lin[i,j,,]
    }
  
  if (type=="response"|type=="HR")
  {
    fam = family(model.y, resp=outcome)
    outcome.pred = fam$linkinv(outcome.linpred)
    if (grepl("zero",fam$family) | grepl("hurdle",fam$family))
      if (dpar.outcome != "mu") outcome.pred = exp(outcome.linpred)/(1+exp(outcome.linpred))
  }
  if (type=="lp"|type=="logHR") outcome.pred = outcome.linpred
  
  res = cal.effects(outcome.pred)
  
  list(effects=res, outcome.pred=outcome.pred, outcome.linpred=outcome.linpred)
}


cal.effects <- function(outcome.pred)
{
  # direct effect: Y(1,M(t)) - Y(0,M(t))
  # control: Y(1,M(0)) - Y(0,M(0))
  direct_control = outcome.pred[2,1,,] - outcome.pred[1,1,,]
  # treated: Y(1,M(1)) - Y(0,M(1))
  direct_treated = outcome.pred[2,2,,] - outcome.pred[1,2,,]
  # mediation effect: Y(t,M(1)) - Y(t,M(0))
  # control: Y(0,M(1)) - Y(0,M(0))
  indirect_control = outcome.pred[1,2,,] - outcome.pred[1,1,,]
  # treated: Y(1,M(1)) - Y(1,M(0))
  indirect_treated = outcome.pred[2,2,,] - outcome.pred[2,1,,]
  # total effect: Y(1,M(1)) - Y(0,M(0))
  total = outcome.pred[2,2,,] - outcome.pred[1,1,,]
  direct = (direct_control + direct_treated)/2
  indirect = (indirect_control + indirect_treated)/2
  
  res = rbind(
    c(mean(indirect_control), median(indirect_control), sd(indirect_control),
      quantile(indirect_control, probs=c(0.025,0.975)),
      2*min(mean(indirect_control<0), mean(indirect_control>0))),
    
    c(mean(indirect_treated), median(indirect_treated), sd(indirect_treated), 
      quantile(indirect_treated, probs=c(0.025,0.975)),
      2*min(mean(indirect_treated<0), mean(indirect_treated>0))),
    
    c(mean(direct_control), median(direct_control), sd(direct_control), 
      quantile(direct_control, probs=c(0.025,0.975)),
      2*min(mean(direct_control<0), mean(direct_control>0))),
    
    c(mean(direct_treated), median(direct_treated), sd(direct_treated),
      quantile(direct_treated, probs=c(0.025,0.975)),
      2*min(mean(direct_treated<0), mean(direct_treated>0))),
    
    c(mean(total), median(total), sd(total), 
      quantile(total, probs=c(0.025,0.975)),
      2*min(mean(total<0), mean(total>0))), 
    
    c(mean(indirect), median(indirect), sd(indirect), 
      quantile(indirect, probs=c(0.025,0.975)),
      2*min(mean(indirect<0), mean(indirect>0))),
    
    c(mean(direct), median(direct), sd(direct), 
      quantile(direct, probs=c(0.025,0.975)),
      2*min(mean(direct<0), mean(direct>0)))
  ) # Bayes p-value: tail probability (see JMbayes), 2*min{pr(b<0), pr(b>0))}
  res[,1:5] = round(res[,1:5], digits=3)
  res[,6] = signif(res[,6], digits=2)
  rownames(res) = c("Indirect_control", "Indirect_treated", 
                    "Direct_control", "Direct_treated",
                    "Total Effect", "Indirect", "Direct")
  colnames(res) = c("Mean", "Median", "Sd", "l-95% CI", "u-95% CI", "Bayes_p")
  
  res
}



#************************************************************************************************

# incorrect interval estimates

brms_mediation1 <- function(model.m, model.y, 
                            mediator = "med.name", treat = "treat.name", 
                            control.value = 0, treat.value = 1)
{
  value = c(control.value, treat.value)
  mediator.pred = array(NA, dim = c(2, nrow(model.m$data)))
  for (i in 1:2)
  {
    dat = model.m$data
    dat[, treat] = value[i]
    mediator.pred[i,] = colMeans(posterior_epred(model.m, newdata=dat))
  }
  outcome.pred = array(NA, dim = c(2, 2, ndraws(model.y), nrow(model.y$data)) )
  for (i in 1:2)
    for (j in 1:2)
    {
      dat = model.y$data
      dat[, treat] = value[i]
      dat[, mediator] = mediator.pred[j,]
      outcome.pred[i,j,,] = posterior_epred(model.y, newdata=dat) 
    }
  
  res = cal.effects(outcome.pred)
  
  list(effects=res, mediator.pred=mediator.pred, outcome.pred=outcome.pred)
}

# too slow!

brms_mediation2 <- function(model.m, model.y, 
                            mediator = "med.name", treat = "treat.name", 
                            control.value = 0, treat.value = 1, 
                            ndraws = 1000)
{
  value = c(control.value, treat.value)
  mediator.pred = array(NA, dim = c(2, ndraws, nrow(model.m$data)))
  for (i in 1:2)
  {
    dat = model.m$data
    dat[, treat] = value[i]
    mediator.pred[i,,] = posterior_epred(model.m, newdata=dat, ndraws=ndraws)
  }
  outcome.pred = array(NA, dim = c(2, 2, ndraws, nrow(model.y$data)) )
  for (i in 1:2)
    for (j in 1:2)
      for (s in 1:ndraws) 
      {
        dat = model.y$data
        dat[, treat] = value[i]
        dat[, mediator] = mediator.pred[j,s,]
        outcome.pred[i,j,s,] = posterior_epred(model.y, newdata=dat, ndraws=1) 
      }
  
  res = cal.effects(outcome.pred)
  
  list(effects=res, mediator.pred=mediator.pred, outcome.pred=outcome.pred)
}

#******************************************************************************
