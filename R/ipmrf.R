ipmrf <-
function(marbolr,da,ntree){
#marbolr: randomforest obtained by Randomforest library
#da: data frame with the predictors only, not responses
#ntree: number of trees in the random forest
#if factors in da: levels should be numbered

mi=marbolr$inbag


if (packageVersion("randomForest") < "4.6-12") {
  warning("This function was not made for randomForest < 4.6-12.", call. = FALSE)
}

#in version of randomForest 4.6.12: The "inbag" component of the randomForest object now records the number of 
 # times an observation is included in the bootstrap sample, instead of just indicating whether it is included or not
 
 mi[mi>1]=1


#percentage of use per tree (if 2 times of a total of 2 splits is more than 2 times of a total of three)
dime=dim(da)
pup=matrix(0,nrow=dime[1],ncol=dime[2]) 




totob=ntree-rowSums(mi) #number of times out-of-bag

for (i in 1:ntree){
ar=getTree(marbolr,k=i) #sin label

#which individuals are out-of-bag
ib=mi[,i] #in bag
#out-of-bag
ob=which(ib==0)

#which variables used in prediction
for(j in 1:length(ob)){
da1=da[ob[j],]
wv=prevtree(ar,da1)
if (is.null(wv)){totob[ob[j]]=totob[ob[j]]-1}
else{

pupi=table(wv)/length(wv)
dond= unique(sort(wv))

pup[ob[j],dond]=pup[ob[j],dond]+ as.numeric(pupi)
}
}
}

pupf=pup
for(h in 1:dime[1]){
pupf[h,]=pup[h,]/totob[h]
}



dimnames(pupf)=dimnames(da)

return(pupf)
}
