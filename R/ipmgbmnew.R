ipmgbmnew=function (marbolr, da, ntree) 
{
  #idx.dep=which(colnames(da)==all.vars(formula)[1])
  da = as.data.frame(da)
  dime = dim(da)
  if (is.null(dime)) {
    pup = matrix(0, nrow = 1, ncol = length(da))
    da = t(as.matrix(da))
    dime = dim(da)
  }
  else {
    pup = matrix(0, nrow = dime[1], ncol = dime[2])
  }
  ob = dim(da)[1]
  totob = rep(ntree, ob)
  for (i in 1:ntree) {
   #ar = getTreeranger(marbolr, k = i)
    ar=getTreegbm(marbolr,da,k=1)
    for (j in 1:ob) {
      da1 = da[j, ]
      wv = prevtree(ar, da1)
      if (is.null(wv)) {
        totob[j] = totob[j] - 1
      }
      else {
        pupi = table(wv)/length(wv)
        dond = unique(sort(wv))
        pup[j, dond] = pup[j, dond] + as.numeric(pupi)
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
