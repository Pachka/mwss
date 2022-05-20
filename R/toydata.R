#' Toy data for mwss examples
#'
#' Toy data gerenated for mwss examples contains baseline data for a 5-ward facility.
#'
#' @docType data
#'
#' @usage data(toydata)
#'
#' @format An object of class \code{"list"}.
#'
#' @keywords datasets
#'
#' @examples
#' data(toydata)
#' str(toydata)
#' nwards <- toydata$nwards
#' ward_names <- toydata$ward_names
#' \donttest{
#' list2env(toydata,envir=.GlobalEnv)
#' matContact <- randomContacts(pop_size_H, ward_names)$contactMat
#' plot_connectivity(matContact, pop_size_P, verbose = FALSE)
#' }

"toydata"
