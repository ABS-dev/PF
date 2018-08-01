#' Package for PF.
#'
#' Includes functions related to prevented fraction. 
#'
#' \tabular{ll}{
#' Package: \tab pf-package\cr
#' Type: \tab Package\cr
#' Version: \tab 9.5.3\cr
#' Date: \tab 2018-08-01\cr
#' License: \tab MIT\cr
#' LazyLoad: \tab yes\cr
#' LazyData: \tab yes\cr
#' }
#'
#' @name PF-package
#' @aliases PF
#' @docType package
#' @author David Siev \email{David.Siev@@aphis.usda.gov}
#' @section Resources:
#' 
#' \itemize{
#' \item GUIDANCE: \url{https://www.aphis.usda.gov/aphis/ourfocus/animalhealth/veterinary-biologics/biologics-regulations-and-guidance/ct_vb_statwi}
#' \item QUICK START: \url{https://github.com/ABS-dev/PF/blob/master/README.md}
#' \item BUG REPORTS: \url{https://github.com/ABS-dev/PF/issues}
#' }
# @examples
# #---------------------------------------------
# # Checking PF package
# #---------------------------------------------
# example(RRsc)
# example(RRstr)
# example(RRmh)
# example(RRor)
# example(phiWt)
# example(tauWt)
# example(rsbWt)
# example(rsb)
# example(RRlsi)
# example(IDRsc)
# example(IDRlsi)
# #---------------------------------------------
# # The next two take a moment to run
# #---------------------------------------------
# example(RRtosst)
# example(RRotsst)
# #---------------------------------------------
# # End examples
# #---------------------------------------------
# invisible()
NA

#' @name New
#' @title New dataset
#' @description New dataset
#' @docType data
#' @format a data frame with 52 observations of the following 3 variables, no NAs
#' \describe{
#' \item{cage}{cage ID. 1 - 26}
#' \item{tx}{treatment. one of 'con' or 'vac'}
#' \item{pos}{numeric indicator of positive response. 0 = FALSE or 1 = TRUE}
#' }
#' @references We need some references
#' @keywords datasets
NA

#' @name Set1
#' @title Set1 dataset
#' @description Set1 dataset
#' @rdname dfSet1
#' @docType data
#' @format a data.frame with 6 observation of the following 4 variables, no NAs
#' \describe{
#' \item{y}{number positive}
#' \item{n}{total number in group \code{tx} x \code{clus}}
#' \item{tx}{treatment 'vac' or 'con'}
#' \item{clus}{cluster ID}
#' }
#' @references We need some references
#' @keywords datasets
NA

#' @name Table6
#' @title Table6 dataset
#' @description Table6 dataset
#' @rdname dfTable6
#' @docType data
#' @format a data.frame with 8 observations of the following 4 variables, no NAs
#' \describe{
#' \item{y}{number positive}
#' \item{n}{total number in group \code{tx} x \code{clus}}
#' \item{tx}{treatment 'a' or 'b'}
#' \item{clus}{cluster ID}
#' }
#' @references Table 1 from Gart (1985) 
#' @keywords datasets
NA

#' @name bird
#' @title bird dataset
#' @description bird dataset
#' @docType data
#' @format a data.frame with 6 observations of the following 4 variables, no NAs
#' \describe{
#' \item{y}{number positive}
#' \item{n}{total number in group \code{tx} x \code{all}}
#' \item{tx}{treatment 'vac' or 'con'}
#' \item{all}{all?}
#' }
#' @references we need some references
#' @keywords datasets
NA

#' @name birdm
#' @title birdm dataset
#' @description birdm dataset
#' @docType data
#' @format a data.frame with 6 observations of the following 4 variables, no NAs
#' \describe{
#' \item{y}{number positive}
#' \item{n}{total number in group \code{tx} x \code{all}}
#' \item{tx}{treatment 'vac' or 'con'}
#' \item{all}{all?}
#' }
#' @references we need some references
#' @keywords datasets
NA

#' @name set1
#' @title set1 dataset
#' @description set1 dataset
#' @format a 3 x 4 matrix of data in \code{\link{Set1}}
#' @references we need some references!
#' @keywords datasets
NA

#' @name table6
#' @description table6 dataset
#' @title table6 dataset
#' @format matrix for of data in \code{\link{Table6}}
#' @keywords datasets
NA

#' @name rat
#' @title rat dataset
#' @description rat dataset
#' @format a data.frame with 32 observations of the following 3 variables, no NAs
#' \describe{
#' \item{y}{number positive}
#' \item{n}{total number}
#' \item{group}{treatment group: 'control' or 'treated'}
#' }
#' @references Weil's rat data (Table 1 of Rao and Scott)
#' @keywords datasets
NA

