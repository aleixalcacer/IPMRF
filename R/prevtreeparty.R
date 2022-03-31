prevtreeparty <-
function(x,da1){
#return which variables used in prediction in tree with package party

wv=c()

names(x) <- c("nodeID", "weights", "criterion", "terminal",
                  "psplit", "ssplits", "prediction", "left", "right")
                  

decision=x$terminal

while(decision==FALSE){

wv=c(wv,as.numeric(x$psplit[1]))

if (is.numeric(da1[,as.numeric(x$psplit[1])])){ #numeric
if(da1[,as.numeric(x$psplit[1])]<=x$psplit[3]){
x=x$left
} else{
x=x$right
}
}else{

if(x$psplit[2] == TRUE){ #ordered factor
if(da1[,as.numeric(x$psplit[1])]<=levels( da1[,as.numeric(x$psplit[1])])[as.numeric(x$psplit[3])]){
x=x$left
} else{
x=x$right
}
} else{ #nominal factor
if(da1[,as.numeric(x$psplit[1])] %in% levels( da1[,as.numeric(x$psplit[1])])[as.logical(x$psplit[3][[1]])]){
x=x$left
} else{
x=x$right
}
}
}



names(x) <- c("nodeID", "weights", "criterion", "terminal",
                  "psplit", "ssplits", "prediction", "left", "right")
decision=x$terminal
}
return(wv)
}
