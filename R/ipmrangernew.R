#' IPM casewise with RF by \pkg{ranger} for new cases, whose responses do not need to be known
#' 
#' The IPM of a new case, i.e. one not used to grow the forest and whose true
#' response does not need to be known, is computed as follows. The new case is put
#' down each of the \emph{ntree} trees in the forest. For each tree, the case goes from the
#' root node to a leaf through a series of nodes. The variable split in these nodes
#' is recorded. The percentage of times a variable is selected along the case's way
#' from the root to the terminal node is calculated for each tree. Note that we do not
#' count the percentage of times a split occurred on variable k in tree t, but only the
#' variables that intervened in the prediction of the case. The IPM for this new case
#' is obtained by averaging those percentages over the \emph{ntree} trees.
#' 
#' The random forest is based on a fast implementation of CART.
#'
#'
#' @param marbolr 
#' Random forest obtained with \code{\link[ranger]{ranger}}.
#' Responses can be of the same type supported by \code{\link[ranger]{ranger}}.
#' Note that  not only numerical or nominal, but also  ordered responses,
#' censored response variables and multivariate responses can be considered with \code{ipmparty}.
#' 
#' @param da
#' Data frame with the predictors only, not responses, for the new cases.
#' Each row corresponds to an observation and each column corresponds to a predictor,
#' which obviously must be the same variables used as predictors in the training set.
#' 
#' @param ntree 
#' Number of trees in the random forest.
#'
#'
#' @return
#' It returns IPM for new cases. It is a matrix with as many rows as cases are in da,
#' and as many columns as predictors are in da.
#' 
#' @export
#' 
#' @details 
#' All details are given in Epifanio (2017).
#' 
#' @references 
#' Pierola, A. and Epifanio, I. and Alemany, S. (2016) An ensemble of ordered logistic regression
#' and random forest for child garment size matching. \emph{Computers & Industrial Engineering},
#' \bold{101}, 455--465.
#' 
#' Epifanio, I. (2017) Intervention in prediction measure: a new approach to assessing variable
#' importance for random forests. \emph{BMC Bioinformatics}, \bold{18}, 230.
#'
#' @author 
#' Stefano Nembrini, Irene Epifanio
#' 
#' @note 
#' See Epifanio (2017) about the parameters of RFs to be used, the advantages and limitations of
#' IPM, and in particular when CART is considered with predictors of different types.
#' 
#' @seealso 
#' \code{\link{ipmparty}},
#' \code{\link{ipmrf}},
#' \code{\link{ipmranger}},
#' \code{\link{ipmpartynew}},
#' \code{\link{ipmrfnew}},
#' \code{\link{ipmgbmnew}}
#' 
#' @examples
#' \dontrun{
#'     if (require("ranger")) {
#'         num.trees = 500
#'         rf <-
#'             ranger(Species ~ .,
#'                    data = iris,
#'                    keep.inbag = TRUE,
#'                    num.trees = num.trees)
#'         
#'         IPM_complete = apply(ipmrangernew(rf, iris[, -5], num.trees), FUN = mean, 2)
#'     }
#' }
ipmrangernew = function(marbolr, da, ntree) {
    #marbolr: randomforest obtained by ranger library
    #da: data frame with the predictors only, not responses
    #ntree: number of trees in the random forest
    
    
    da = as.data.frame(da)
    
    #percentage of use per tree (if 2 times of a total of 2 splits is more than 2 times of a total of three)
    dime = dim(da)
    if (is.null(dime)) {
        pup = matrix(0, nrow = 1, ncol = length(da))
        da = t(as.matrix(da))
        dime = dim(da)
    } else{
        pup = matrix(0, nrow = dime[1], ncol = dime[2])
    }
    
    
    
    #totob=ntree
    ob = dim(da)[1]
    
    totob = rep(ntree, ob)
    
    
    for (i in 1:ntree) {
        #ar=getTree(marbolr,k=i) #sin label
        ar = getTreeranger(marbolr, k = i)
        # ob=dim(da)[1]
        
        #which variables used in prediction
        for (j in 1:ob) {
            da1 = da[j, ]
            wv = prevtree(ar, da1)
            if (is.null(wv)) {
                totob[j] = totob[j] - 1
            }
            else{
                pupi = table(wv) / length(wv)
                dond = unique(sort(wv))
                
                pup[j, dond] = pup[j, dond] + as.numeric(pupi)
            }
        }
    }
    
    #pupf=pup/totob
    pupf = pup
    for (h in 1:dime[1]) {
        pupf[h, ] = pup[h, ] / totob[h]
    }
    
    
    
    
    dimnames(pupf) = dimnames(da)
    
    return(pupf)
}
