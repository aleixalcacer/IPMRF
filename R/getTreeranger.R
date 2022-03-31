getTreeranger=function(rfobj, k=1){
  
  
  status=rep(1,length(rfobj$forest$child.nodeIDs[[k]][[1]]))  
  idx1=which(rfobj$forest$child.nodeIDs[[k]][[1]]==0)
  idx2=which(rfobj$forest$child.nodeIDs[[k]][[2]]==0)
  status[idx1]=-1
  status[idx2]=-1
  rfobj$forest$child.nodeIDs[[k]][[1]]=rfobj$forest$child.nodeIDs[[k]][[1]]+1
  rfobj$forest$child.nodeIDs[[k]][[1]][idx1]=0
  rfobj$forest$child.nodeIDs[[k]][[2]]=rfobj$forest$child.nodeIDs[[k]][[2]]+1
  rfobj$forest$child.nodeIDs[[k]][[2]][idx1]=0
  
  tree=cbind(rfobj$forest$child.nodeIDs[[k]][[1]], rfobj$forest$child.nodeIDs[[k]][[2]],
             rfobj$forest$split.varIDs[[k]],
             rfobj$forest$split.values[[k]],status)
  tree=as.matrix(tree)
  rownames(tree) <- 1:nrow(tree)
  
  colnames(tree)=c("left daughter", "right daughter", "split var", "split point","status")
  return(tree)
}
