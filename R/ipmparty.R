ipmparty <-
function(marbol,da,ntree){
#marbol: randomforest obtained by party
#da: data frame with the predictors only, not responses
#ntree: number of trees in the random forest


#percentage of use per tree (if 2 times of a total of 2 splits is more than 2 times of a total of three)
dime=dim(da)
pup=matrix(0,nrow=dime[1],ncol=dime[2]) 

#w=marbol@where #for version  1.0.25
w=marbol@weights #for version 1.2.3

mi=matrix(0,ncol=ntree,nrow=dim(da)[1])

for(i in 1:ntree){
mi[which( w[[i]]!=0),i]=1 # 0 when out-of-bag
}

totob=ntree-rowSums(mi) #number of times out-of-bag

for (i in 1:ntree){

ar=marbol@ensemble[[i]]

#which individuals are out-of-bag
ib=mi[,i] #in bag
#out-of-bag
ob=which(ib==0)




#which variables used for predicting
for(j in 1:length(ob)){
da1=da[ob[j],]
wv=prevtreeparty(ar,da1)

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
