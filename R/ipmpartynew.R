#' IPM casewise with CIT-RF by \pkg{party} for new cases, whose responses do not need to be known
#'
#'The IPM of a new case, i.e. one not used to grow the forest and whose true
#'response does not need to be known, is computed as follows. The new case is put
#'down each of the \emph{ntree} trees in the forest. For each tree, the case goes from the
#'root node to a leaf through a series of nodes. The variable split in these nodes
#'is recorded. The percentage of times a variable is selected along the case's way
#'from the root to the terminal node is calculated for each tree. Note that we do not
#'count the percentage of times a split occurred on variable k in tree t, but only the
#'variables that intervened in the prediction of the case. The IPM for this new case
#'is obtained by averaging those percentages over the \emph{ntree} trees. 
#'The random forest is based on CIT (Conditional Inference Trees)
#'
#' @param marbol
#' Random forest obtained with \code{\link[party]{cforest}}.
#' Responses in the training set can be of the same type supported by
#' \code{\link[party]{cforest}}, not only numerical or nominal, but also
#' ordered responses, censored response variables and multivariate responses.
#' 
#' @param da
#' Data frame with the predictors only, not responses, for the new cases.
#' Each row corresponds to an observation and each column corresponds to a predictor,
#' which obviously must be the same variables used as predictors in the training set.
#' Predictors can be numeric, nominal or an ordered factor.
#' 
#' @param ntree
#' Number of trees in the random forest.
#'
#' @return
#' It returns IPM for new cases. It is a matrix with as many rows as cases are in da,
#' and as many columns as predictors are in da. IPM can be estimated for any kind of RF
#' computed by \code{\link[party]{cforest}}, including multivariate RF.
#' 
#' @export
#'
#' @details
#' All details are given in Epifanio (2017).
#'
#' @references
#' Pierola, A. and Epifanio, I. and Alemany, S. (2016) An ensemble of ordered
#' logistic regression and random forest for child garment size matching.
#' \emph{Computers & Industrial Engineering}, \bold{101}, 455--465.
#' 
#' Epifanio, I. (2017) Intervention in prediction measure: a new approach to assessing
#' variable importance for random forests. \emph{BMC Bioinformatics}, \bold{18}, 230.
#'
#' @author 
#' Irene Epifanio
#' 
#' @note 
#' See Epifanio (2017) about advantages and limitations of IPM, and about the
#' parameters to be used in \code{\link[party]{cforest}}.
#' 
#' @seealso 
#' \code{\link{ipmparty}},
#' \code{\link{ipmrf}},
#' \code{\link{ipmranger}},
#' \code{\link{ipmrfnew}},
#' \code{\link{ipmrangernew}},
#' \code{\link{ipmgbmnew}}
#' 
#' @examples
#' \dontrun{
#'     #Note: more examples can be found at
#'     #https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1650-8
#'     
#'     ## -------------------------------------------------------
#'     ## Example from \code{\link[party]{varimp}} in \pkg{party}
#'     ## Classification RF
#'     ## -------------------------------------------------------
#'     
#'     
#'     library(party)
#'     
#'     
#'     #IMP based on CIT-RF (party package)
#'     ntree = 50
#'     #readingSkills: data from party package
#'     da = readingSkills[, 1:3]
#'     set.seed(290875)
#'     readingSkills.cf3 <- cforest(score ~ .,
#'                                  data = readingSkills,
#'                                  control = cforest_unbiased(mtry = 3, ntree = 50))
#'     
#'     #new case
#'     nativeSpeaker = 'yes'
#'     age = 8
#'     shoeSize = 28
#'     da1 = data.frame(nativeSpeaker, age, shoeSize)
#'     
#'     #IPM case-wise computed for new cases for party package
#'     pupfn = ipmpartynew(readingSkills.cf3, da1, ntree)
#' pupfn
#' }
#' 
ipmpartynew <-
    function(marbol, da, ntree) {
        #marbol: randomforest obtained by party
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
        
        
        #w=marbol@get_where(newdata = as.data.frame(da))
        
        
        #if(is.null(dime)){mi=matrix(0,ncol=ntree,nrow=1)
        #}else{
        #mi=matrix(0,ncol=ntree,nrow=dim(da)[1])}
        
        #for(i in 1:ntree){
        #mi[which( w[[i]]!=0),i]=1 # 0 when is out-of-bag
        #}
        
        #totob=ntree-rowSums(mi) #number of times out-of-bag
        
        ob = dim(da)[1]
        
        totob = rep(ntree, ob)
        
        
        for (i in 1:ntree) {
            ar = marbol@ensemble[[i]]
            
            
            
            #which variables used to predict
            for (j in 1:ob) {
                #da1=as.matrix(da[j,],nrow=1,ncol=dim(da)[2]) #attention, depending in format of data, use this or line below (commenting)
                da1 = da[j, ]
                wv = prevtreeparty(ar, da1)
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
