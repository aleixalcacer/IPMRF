#' IPM casewise with CIT-RF by \pkg{party} for OOB samples
#' 
#' The IPM for a case in the training set is calculated by considering and averaging
#' over only the trees where the case belongs to the OOB set. The case is put
#' down each of the trees  where the case belongs to the OOB set. For each tree,
#' the case goes from the root node to a leaf through a series of nodes. The variable split
#' in these nodes is recorded. The percentage of times a variable is selected along the
#' case's way from the root to the terminal node is calculated for each tree. Note that we
#' do not count the percentage of times a split occurred on variable k in tree t, but only
#' the variables that intervened in the prediction of the case. The IPM for this case is
#' obtained by averaging those percentages over only the trees where the case belongs to the
#' OOB set. The random forest is based on CIT (Conditional Inference Trees).
#'
#' @param marbol
#' Random forest obtained with \code{\link[party]{cforest}}. Responses can be of the same
#' type supported by \code{\link[party]{cforest}}, not only numerical or nominal, but also
#' ordered responses, censored response variables and multivariate responses.
#' @param da
#' Data frame with the predictors only, not responses, of the training set
#' used for computing \emph{marbol}. Each row corresponds to an observation and each column
#' corresponds to a predictor. Predictors can be numeric, nominal or an ordered factor.
#' @param ntree 
#' Number of trees in the random forest.
#' 
#' @return
#' It returns IPM for cases in the training set. It is estimated when they are OOB observations.
#' It is a matrix with as many rows as cases are in da, and as many columns as predictors are in da.
#' IPM can be estimated for any kind of RF computed by \code{\link[party]{cforest}}, including
#' multivariate RF.
#' 
#' @export
#'
#' @details 
#' All details are given in Epifanio (2017).
#' 
#' @references 
#' Pierola, A. and Epifanio, I. and Alemany, S. (2016) An ensemble of ordered logistic regression and random forest for child garment size matching. \emph{Computers & Industrial Engineering}, \bold{101}, 455--465.
#' 
#' Epifanio, I. (2017) Intervention in prediction measure: a new approach to assessing variable importance for random forests. \emph{BMC Bioinformatics}, \bold{18}, 230.
#' 
#' @author
#' Irene Epifanio
#' 
#' @note 
#' See Epifanio (2017) about advantages and limitations of IPM, and about the parameters to be
#' used in \code{\link[party]{cforest}}.
#' 
#' @seealso
#' \code{\link{ipmpartynew}},
#' \code{\link{ipmrf}},
#' \code{\link{ipmranger}},
#' \code{\link{ipmrfnew}},
#' \code{\link{ipmrangernew}},
#' \code{\link{ipmgbmnew}}
#' 
#' @examples
#'
#' #Note: more examples can be found at 
#' #https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1650-8
#' 
#' \dontrun{
#'     ## -------------------------------------------------------
#'     ## Example from \code{\link[party]{varimp}} in \pkg{party}
#'     ## Classification RF
#'     ## -------------------------------------------------------
#' 
#'     library(party)
#'     
#'     #from help in varimp by party package
#'     set.seed(290875)
#'     readingSkills.cf <- cforest(score ~ ., data = readingSkills,
#'                                 control = cforest_unbiased(mtry = 2, ntree = 50))
#'     
#'     # standard importance
#'     varimp(readingSkills.cf)
#'     
#'     # the same modulo random variation
#'     varimp(readingSkills.cf, pre1.0_0 = TRUE)
#'     
#'     # conditional importance, may take a while...
#'     varimp(readingSkills.cf, conditional = TRUE)
#' }
#' 
#' #IMP based on CIT-RF (party package)
#' library(party)
#' 
#' ntree<-50
#' #readingSkills: data from party package
#' da<-readingSkills[,1:3] 
#' set.seed(290875)
#' readingSkills.cf3 <- cforest(score ~ ., data = readingSkills,
#'                              control = cforest_unbiased(mtry = 3, ntree = 50))
#' 
#' #IPM case-wise computed with OOB with party
#' pupf<-ipmparty(readingSkills.cf3 ,da,ntree)
#' 
#' #global IPM
#' pua<-apply(pupf,2,mean) 
#' pua
#' 
#' 
#' \dontrun{
#'     ## -------------------------------------------------------
#'     ## Example from \code{\link[randomForestSRC]{var.select}} in \pkg{randomForestSRC} 
#'     ## Multivariate mixed forests
#'     ## -------------------------------------------------------
#'     
#'     if(require("randomForestSRC")) {
#'     
#'         #from help in var.select by randomForestSRC package
#'         mtcars.new <- mtcars
#'         mtcars.new$cyl <- factor(mtcars.new$cyl)
#'         mtcars.new$carb <- factor(mtcars.new$carb, ordered = TRUE)
#'         mv.obj <- rfsrc(cbind(carb, mpg, cyl) ~., data = mtcars.new,
#'                         importance = TRUE)
#'         var.select(mv.obj, method = "vh.vimp", nrep = 10) 
#'         
#'         #different variables are selected if var.select is repeated
#'     }
#' }
#' 
#' #IMP based on CIT-RF (party package)
#' if(require("randomForestSRC")) {
#'     mtcars.new <- mtcars
#'     
#'     ntree<-500
#'     da<-mtcars.new[,3:10] 
#'     mc.cf <- cforest(carb+ mpg+ cyl ~., data = mtcars.new,
#'                      control = cforest_unbiased(mtry = 8, ntree = 500))
#'     
#'     #IPM case-wise computing with OOB with party
#'     pupf<-ipmparty(mc.cf ,da,ntree) 
#'     
#'     #global IPM
#'     pua<-apply(pupf,2,mean) 
#'     pua
#'     
#'     #disp and hp are consistently selected as more important if repeated
#' }
#' 
ipmparty <- function(marbol, da, ntree) {
    #marbol: randomforest obtained by party
    #da: data frame with the predictors only, not responses
    #ntree: number of trees in the random forest
    
    
    #percentage of use per tree (if 2 times of a total of 2 splits is more than 2 times of a total of three)
    dime = dim(da)
    pup = matrix(0, nrow = dime[1], ncol = dime[2])
    
    #w=marbol@where #for version  1.0.25
    w = marbol@weights #for version 1.2.3
    
    mi = matrix(0, ncol = ntree, nrow = dim(da)[1])
    
    for (i in 1:ntree) {
        mi[which(w[[i]] != 0), i] = 1 # 0 when out-of-bag
    }
    
    totob = ntree - rowSums(mi) #number of times out-of-bag
    
    for (i in 1:ntree) {
        ar = marbol@ensemble[[i]]
        
        #which individuals are out-of-bag
        ib = mi[, i] #in bag
        #out-of-bag
        ob = which(ib == 0)
        
        
        
        
        #which variables used for predicting
        for (j in 1:length(ob)) {
            da1 = da[ob[j],]
            wv = prevtreeparty(ar, da1)
            
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
        pupf[h,] = pup[h,] / totob[h]
    }
    
    
    dimnames(pupf) = dimnames(da)
    
    return(pupf)
    
}
