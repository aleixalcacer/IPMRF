\name{ipmrfnew}
\alias{ipmrfnew}
\title{
IPM casewise with CART-RF by \pkg{randomForest} for new cases, whose responses do not need to be known}
\description{
The IPM of a new case, i.e. one not used to grow the forest and whose true
response does not need to be known, is computed as follows. The new case is put
down each of the \emph{ntree} trees in the forest. For each tree, the case goes from the
root node to a leaf through a series of nodes. The variable split in these nodes
is recorded. The percentage of times a variable is selected along the case's way
from the root to the terminal node is calculated for each tree. Note that we do not
count the percentage of times a split occurred on variable k in tree t, but only the
variables that intervened in the prediction of the case. The IPM for this new case
is obtained by averaging those percentages over the \emph{ntree} trees.

The random forest is based on CART}
\usage{
ipmrfnew(marbolr, da, ntree)
}
\arguments{
  \item{marbolr}{
Random forest obtained with \code{\link[randomForest]{randomForest}}. Responses can be of the same type supported by \code{\link[randomForest]{randomForest}}. Note that  not only numerical or nominal, but also  ordered responses,  censored response variables and multivariate responses can be considered with \code{ipmparty}.}
  \item{da}{
Data frame with the predictors only, not responses, for the new cases. Each row corresponds to an observation and each column corresponds to a predictor, which obviously must be the same variables used as predictors in the training set.}
  \item{ntree}{
Number of trees in the random forest.}
}
\details{
All details are given in Epifanio (2017).
}
\value{
It returns IPM for new cases. It is a matrix with as many rows as cases are in da, and as many columns as predictors are in da.
}
\references{
Pierola, A. and Epifanio, I. and Alemany, S. (2016) An ensemble of ordered logistic regression and random forest for child garment size matching. \emph{Computers & Industrial Engineering}, \bold{101}, 455--465.

Epifanio, I. (2017) Intervention in prediction measure: a new approach to assessing variable importance for random forests. \emph{BMC Bioinformatics}, \bold{18}, 230.
}
\author{
Irene Epifanio
}
\note{
See Epifanio (2017) about the parameters of RFs to be used, the advantages and limitations of IPM, and in particular when CART is considered with predictors of different types.
}
\seealso{
\code{\link{ipmparty}}, \code{\link{ipmrf}}, \code{\link{ipmranger}}, \code{\link{ipmpartynew}}, \code{\link{ipmrangernew}},\code{\link{ipmgbmnew}}
}
\examples{
\dontrun{
#Note: more examples can be found at 
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1650-8


library(mlbench)
#data used by Breiman, L.: Random forests. Machine Learning 45(1), 5--32 (2001)
data(PimaIndiansDiabetes2) 
Diabetes <- na.omit(PimaIndiansDiabetes2)

set.seed(2016)
require(randomForest)
ri<- randomForest(diabetes  ~ ., data=Diabetes, ntree=500, importance=TRUE,
keep.inbag=TRUE,replace = FALSE) 

#new cases
da1=rbind(apply(Diabetes[Diabetes[,9]=='pos',1:8],2,mean),
apply(Diabetes[Diabetes[,9]=='neg',1:8],2,mean))


#IPM case-wise computed for new cases for randomForest package
ntree=500
pupfn=ipmrfnew(ri, as.data.frame(da1),ntree)
pupfn
}
}
\keyword{tree}
\keyword{multivariate}
