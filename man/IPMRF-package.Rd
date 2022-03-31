\name{IPMRF-package}
\alias{IPMRF-package}
\alias{IPMRF}
\docType{package}
\title{Intervention in Prediction Measure (IPM) for Random Forests}
\description{
It computes IPM for assessing variable importance for
random forests. See I. Epifanio (2017). Intervention in prediction measure: a new approach to assessing variable importance for random forests. BMC Bioinformatics.}
\details{
 \tabular{ll}{
        Package: \tab IPMRF\cr
        Type: \tab Package\cr
        Version: \tab 1.2\cr
        Date: \tab 2017-08-09\cr
    }
 Main Functions:
 \itemize{
  \item {ipmparty:} {IPM casewise with CIT-RF by \pkg{party} for OOB samples}
	\item {ipmpartynew:} {IPM casewise with CIT-RF by \pkg{party} for new samples}
\item {ipmrf:} {IPM casewise with CART-RF by \pkg{randomForest} for OOB samples}
	\item {ipmrfnew:} {IPM casewise with CART-RF by \pkg{randomForest} for new samples}
	\item {ipmranger:} {IPM casewise with RF by \pkg{ranger} for OOB samples}
	\item {ipmrangernew:} {IPM casewise with RF by \pkg{ranger} for new samples}
	\item {ipmgbmnew:} {IPM casewise with GBM by \pkg{gbm} for new samples}
 }
}
\author{Irene Epifanio, Stefano Nembrini}
\references{
Pierola, A. and Epifanio, I. and Alemany, S. (2016) An ensemble of ordered logistic regression and random forest for child garment size matching. \emph{Computers & Industrial Engineering}, \bold{101}, 455--465.

Epifanio, I. (2017) Intervention in prediction measure: a new approach to assessing variable importance for random forests. \emph{BMC Bioinformatics}, \bold{18}, 230.
}
\keyword{Variable importance measure}
\keyword{Random forest}
\keyword{Generalized Boosted Regression}
\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1650-8
}
