ipmpartynew <-
function(marbol,da,ntree){
#marbol: randomforest obtained by party
#da: data frame with the predictors only, not responses
#ntree: number of trees in the random forest

da=as.data.frame(da)

#percentage of use per tree (if 2 times of a total of 2 splits is more than 2 times of a total of three)
dime=dim(da)
if(is.null(dime)){pup=matrix(0,nrow=1,ncol=length(da))
da=t(as.matrix(da))
dime=dim(da)
}else{
pup=matrix(0,nrow=dime[1],ncol=dime[2])} 


#w=marbol@get_where(newdata = as.data.frame(da))


#if(is.null(dime)){mi=matrix(0,ncol=ntree,nrow=1)
#}else{
#mi=matrix(0,ncol=ntree,nrow=dim(da)[1])}

#for(i in 1:ntree){
#mi[which( w[[i]]!=0),i]=1 # 0 when is out-of-bag
#}

#totob=ntree-rowSums(mi) #number of times out-of-bag

ob=dim(da)[1]

totob=rep(ntree,ob)


for (i in 1:ntree){

ar=marbol@ensemble[[i]]



#which variables used to predict
for(j in 1:ob){
#da1=as.matrix(da[j,],nrow=1,ncol=dim(da)[2]) #attention, depending in format of data, use this or line below (commenting)
da1=da[j,]
wv=prevtreeparty(ar,da1) 
if (is.null(wv)){totob[j]=totob[j]-1}
else{

pupi=table(wv)/length(wv)
dond= unique(sort(wv))

pup[j,dond]=pup[j,dond]+ as.numeric(pupi)
}
}
}


#pupf=pup/totob

pupf=pup
for(h in 1:dime[1]){
pupf[h,]=pup[h,]/totob[h]
}


dimnames(pupf)=dimnames(da)

return(pupf)

}
