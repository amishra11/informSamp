##### Weighted Andersen-Gill Model ######
library(rootSolve)

### Score function for stratified, time-varying covariate version of Petra's model
#beta = parameter for predictors of interest and adjustment variables
#delta = stabilization parameter fo weights
# x = predictor of interest (e.g vaginal bleeding) (n x p)
# z = adjustment covariates (e.g. arm, age, baseline STI) (n x p)
# weight = IIRW weights (already inversed) (n x 1)
# nInd = indicator of event occurence (n x1)
# time = n x 1 vector of visit for each particiant
# strata = stratifcation variables (e.g. site )
weighted.AG.score <- function(beta,delta,x,z,weight,nInd,id,time,data,strata){
  beta <- as.matrix(beta)
  delta <- as.matrix(delta)
  
  design <- model.matrix( ~ x + z)
  if(ncol(design)>1) { design <- design[,2:ncol(design)] }#drop intercept column if there is one
  
  ## loop for stratum
  score.strat <- matrix(NA,nrow=length(levels(strata)), ncol=ncol(design))
  for (k in 1:nrow(score.strat)){
    #subset to stratum ubdex
    strLevel <- levels(strata)[k]
    strata.index <- strata==strLevel
    pred <- design[strata.index,]
    id.strat <- id[strata.index]
    time.strat <- time[strata.index]
    nInd.strat <- nInd[strata.index]
    weight.strat <- weight[strata.index]
    
    ##hazard for everyone in stratum at each visit 
    # if censored already removed from dataset (i.e. not in risk set)
    link <- data.frame("id"=id.strat,"time"=time.strat,"link"=exp(pred %*% beta)* exp(pred %*% delta))
    colnames(link) <- c("id","time","link")
    
    #At each visit sum hazard of everyone in risk set 
    #this will be the denominator for the fraction in the score function
    av2.dem <- aggregate(link$link,by=list(link$time),sum) 
    names(av2.dem) <- c("visit","dem")
  
    u <- unique(id.strat)
    ## this loop calculates the fraction in the score function for every individual 
    ## Each individual as an m x p matrix, where m = number of visits and p = number of predictors
    av2 <- matrix(0,nrow=nrow(av2.dem),ncol=ncol(pred))
    for(i in 1:length(u)){
      #getting subject specific info
      id.index <- id.strat==u[i]
      pred.tmp <- matrix(pred[id.index,],ncol=ncol(pred))
      link.tmp <- as.matrix(link$link[id.index])
    
      ## multiplies each column in x.tmp by link 
      av2.num.tmp <- sweep(pred.tmp,1,link.tmp,FUN = "*") #numerator of fraction (rows=visits, col=predictor)
      
      #if participant has fewer visits than end of study (33 visits) then add rows of 0 to bottom so 
      #I can easily add matrices together
      add.row <- length(av2.dem$visit) - length(link.tmp)
      av2.num.tmp2 <- rbind(av2.num.tmp,matrix(0,nrow = add.row,ncol=ncol(pred)) )
      
      av2.tmp <- sweep(av2.num.tmp2,1,av2.dem$dem, FUN = "/") #fraction in score function for person i
      av2 <- av2 + av2.tmp #fraction in score function (m x p matrix, where m is total visits in study)
    }
  
    ## second for loop to get the difference between xi(t) and fraction ###
    score <- rep(0,ncol=ncol(pred))
    for(j in 1:length(u)){
      id.index <- id.strat==u[j]
      pred.tmp <- matrix(pred[id.index,],ncol=ncol(pred))
    
      ## if not the same length as total number of visits then add some 0s to the end for easy addition/subtraction
      add.row <- length(av2.dem$visit) - nrow(pred.tmp)
      pred.tmp2 <- rbind(pred.tmp,matrix(0,nrow=add.row,ncol=ncol(pred)))
      event.tmp <- as.matrix(c(nInd.strat[id.index],rep(0,add.row)))
    
      #difference in score function
      diff <- pred.tmp2 - av2
    
      #multiple difference by weight at time t
      wt.tmp <- as.matrix(c(weight.strat[id.index],rep(0,add.row)))
      int.tmp <- sweep(diff,1,wt.tmp, FUN = "*")
    
      #only keep those observations at which event happend (i.e. dN_i(t)=1 ) 
      int.tmp2 <- matrix(int.tmp[event.tmp==1,],ncol=ncol(pred)) #this is the integral in the score function
      #it is an m x p matric (m= # visits, p = # predictors)
    
      #sum across time
      score.tmp <- apply(int.tmp2,2,sum) #score contribution for person i
      score <- score + score.tmp #scire contribution for stratum k
    }
    score.strat[k,] <- score #score for each stratum 
    #score.strat is a k x p matrix, k = numer of stratum, p = number of predictors
  }
  res <- apply(score.strat,2,sum) #sum across stratum
  # res is a px1 vector. Want to find beta such that res = 0 for all vector entries
  return(res)
}
