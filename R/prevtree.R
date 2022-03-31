prevtree <-
function(ar,da1){
#return which variables used in prediction in tree with package randomForest

wv=c()
wrow=1
decision=ar[1,5]

while(decision != -1){

wv=c(wv,ar[wrow,3])

if (is.numeric(da1[,ar[wrow,3]])){#numeric

if(da1[ar[wrow,3]]<=ar[wrow,4]){
wrow=ar[wrow,1]
} else{
wrow=ar[wrow,2]
}
decision=ar[wrow,5]}
else{#factor
le=levels(da1[,ar[wrow,3]])
lel=length(le)
cod=number2binary(ar[wrow,4], lel)
if(cod[which(le== da1[,ar[wrow,3]])]){
wrow=ar[wrow,1]
} else{
wrow=ar[wrow,2]
}

decision=ar[wrow,5]
}


}
return(wv)
}
