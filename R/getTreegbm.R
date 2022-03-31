getTreegbm = function (gbm, X, k) {
  #library(gbm)
  list <- pretty.gbm.tree(gbm, i.tree = k)
  a <- list
  v2int <- function(v) {
    sum((-v + 1) / 2 * 2 ^ seq(0, (length(v) - 1), 1))
  }
  splitBin = sapply(gbm$c.splits, v2int)
  rownames(a) <- 1:nrow(a)
  a$status <- rep(1, length(a$SplitVar))
  
  
  a <- a[, c("LeftNode",
             "RightNode",
             "MissingNode",
             "SplitVar",
             "SplitCodePred",
             "status")]
  a[which(a[, "SplitVar"] >= 0), c("SplitVar", "LeftNode",
                                   "RightNode", "MissingNode")] <-
    a[which(a[, "SplitVar"] >=
              0), c("SplitVar", "LeftNode", "RightNode", "MissingNode")] +   1
  #ix <- a$MissingNode[which(a$SplitVar <0)]
  #if (length(ix) > 0)
  #  a$status[ix] <- -1
  a <-
    a[, c("LeftNode",
          "RightNode",
          "SplitVar",
          "SplitCodePred",
          "status")]
  cat <- which(sapply(X, is.factor) & !sapply(X, is.ordered))
  ix <- which(a[, "SplitVar"] %in% cat)
  for (i in ix)
    a[i, "SplitCodePred"] <- splitBin[a[i,
                                        "SplitCodePred"] + 1]
  a$LeftNode[which((a$SplitVar) == -1)] = 0
  a$RightNode[which((a$SplitVar) == -1)] = 0
  a$SplitVar[which((a$SplitVar) == -1)] = 0
  a$status[which((a$SplitVar) == 0)] = -1
  colnames(a) <-
    c("left daughter",
      "right daughter",
      "split var",
      "split point",
      "status")
  
  
  return(a)
}
