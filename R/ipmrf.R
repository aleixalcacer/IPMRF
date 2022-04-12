#' IPM casewise with CART-RF by \pkg{randomForest} for OOB samples
#' 
#' The IPM for a case in the training set is calculated by considering and averaging
#' over only the trees where the case belongs to the OOB set. The case is put
#' down each of the trees  where the case belongs to the OOB set.
#' For each tree, the case goes from the root node to a leaf through a series of nodes.
#' The variable split in these nodes is recorded. The percentage of times a variable is
#' selected along the case's way from the root to the terminal node is calculated for each tree.
#' Note that we do not count the percentage of times a split occurred on variable $k$ in tree $t$,
#' but only the variables that intervened in the prediction of the case. The IPM for this case is
#' obtained by averaging those percentages over only the trees where the case belongs to the
#' OOB set. 
#' 
#' The random forest is based on CART.
#'
#' @param marbolr 
#' Random forest obtained with \code{\link[randomForest]{randomForest}}.
#' Responses can be of the same type supported by \code{\link[randomForest]{randomForest}}.
#' Note that  not only numerical or nominal, but also  ordered responses, censored response
#' variables and multivariate responses can be considered with \code{ipmparty}.
#' 
#' @param da 
#' Data frame with the predictors only, not responses, of the training set used for computing
#' \emph{marbolr}. Each row corresponds to an observation and each column corresponds to a
#' predictor. Predictors can be numeric, nominal or an ordered factor.
#' 
#' @param ntree 
#' Number of trees in the random forest.
#'
#' @return
#' It returns IPM for cases in the training set. It is estimated when they are OOB observations.
#' It is a matrix with as many rows as cases are in da, and as many columns as predictors are in da.
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
#' Irene Epifanio
#' 
#' @note
#' See Epifanio (2017) about the parameters of RFs to be used, the advantages and limitations of
#' IPM, and in particular when CART is considered with predictors of different types.
#' 
#' @seealso 
#' \code{\link{ipmparty}},
#' \code{\link{ipmranger}},
#' \code{\link{ipmpartynew}},
#' \code{\link{ipmrfnew}},
#' \code{\link{ipmrangernew}},
#' \code{\link{ipmgbmnew}}
#' 
#' @examples
#' 
#' \dontrun{
#'     #Note: more examples can be found at
#'     #https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1650-8
#'     
#'     if (require("mlbench")) {
#'         #data used by Breiman, L.: Random forests. Machine Learning 45(1), 5--32 (2001)
#'         data(PimaIndiansDiabetes2)
#'         Diabetes <- na.omit(PimaIndiansDiabetes2)
#'         
#'         set.seed(2016)
#'         require(randomForest)
#'         ri <-
#'             randomForest(
#'                 diabetes  ~ .,
#'                 data = Diabetes,
#'                 ntree = 500,
#'                 importance = TRUE,
#'                 keep.inbag = TRUE,
#'                 replace = FALSE
#'             )
#'         
#'         #GVIM and PVIM (CART-RF)
#'         im = importance(ri)
#'         im
#'         #rank
#'         ii = apply(im, 2, rank)
#'         ii
#'         
#'         #IPM based on CART-RF (randomForest package)
#'         da = Diabetes[, 1:8]
#'         ntree = 500
#'         #IPM case-wise computed with OOB
#'         pupf = ipmrf(ri, da, ntree)
#'         
#'         #global IPM
#'         pua = apply(pupf, 2, mean)
#'         pua
#'         
#'         #IPM by classes
#'         attach(Diabetes)
#'         puac = matrix(0, nrow = 2, ncol = dim(da)[2])
#'         puac[1, ] = apply(pupf[diabetes == 'neg', ], 2, mean)
#'         puac[2, ] = apply(pupf[diabetes == 'pos', ], 2, mean)
#'         colnames(puac) = colnames(da)
#'         rownames(puac) = c('neg', 'pos')
#'         puac
#'         
#'         #rank IPM
#'         #global rank
#'         rank(pua)
#'         #rank by class
#'         apply(puac, 1, rank)
#'     }
#' }
ipmrf <-
    function(marbolr, da, ntree) {
        #marbolr: randomforest obtained by Randomforest library
        #da: data frame with the predictors only, not responses
        #ntree: number of trees in the random forest
        #if factors in da: levels should be numbered
        
        mi = marbolr$inbag
        
        
        if (utils::packageVersion("randomForest") < "4.6-12") {
            warning("This function was not made for randomForest < 4.6-12.",
                    call. = FALSE)
        }
        
        #in version of randomForest 4.6.12: The "inbag" component of the randomForest object now records the number of
        # times an observation is included in the bootstrap sample, instead of just indicating whether it is included or not
        
        mi[mi > 1] = 1
        
        
        #percentage of use per tree (if 2 times of a total of 2 splits is more than 2 times of a total of three)
        dime = dim(da)
        pup = matrix(0, nrow = dime[1], ncol = dime[2])
        
        
        
        
        totob = ntree - rowSums(mi) #number of times out-of-bag
        
        for (i in 1:ntree) {
            ar = randomForest::getTree(marbolr, k = i) #sin label
            
            #which individuals are out-of-bag
            ib = mi[, i] #in bag
            #out-of-bag
            ob = which(ib == 0)
            
            #which variables used in prediction
            for (j in 1:length(ob)) {
                da1 = da[ob[j], ]
                wv = prevtree(ar, da1)
                if (is.null(wv)) {
                    totob[ob[j]] = totob[ob[j]] - 1
                }
                else{
                    pupi = table(wv) / length(wv)
                    dond = unique(sort(wv))
                    
                    pup[ob[j], dond] = pup[ob[j], dond] + as.numeric(pupi)
                }
            }
        }
        
        pupf = pup
        for (h in 1:dime[1]) {
            pupf[h, ] = pup[h, ] / totob[h]
        }
        
        
        
        dimnames(pupf) = dimnames(da)
        
        return(pupf)
    }
