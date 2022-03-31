ipmranger=function (marbolr, da, ntree) 
{
  mi=matrix(unlist(marbolr$inbag), ncol = ntree, byrow = TRUE)
  
  mi[mi > 1] = 1
  dime = dim(da)
  pup = matrix(0, nrow = dime[1], ncol = dime[2])
  totob = ntree - rowSums(mi)
  for (i in 1:ntree) {
    #ar = getTree(marbolr, k = i)
    ar = getTreeranger(marbolr, k = i)
    ib = mi[, i]
    ob = which(ib == 0)
    for (j in 1:length(ob)) {
      da1 = da[ob[j], ]
      wv = prevtree(ar, da1)
			if (is.null(wv)){totob[ob[j]]=totob[ob[j]]-1}
		else{

      pupi = table(wv)/length(wv)
      dond = unique(sort(wv))
      pup[ob[j], dond] = pup[ob[j], dond] + as.numeric(pupi)
			}
    }
  }
  pupf = pup
  for (h in 1:dime[1]) {
    pupf[h, ] = pup[h, ]/totob[h]
  }
  dimnames(pupf) = dimnames(da)
  return(pupf)
}
