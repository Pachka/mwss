#' Generate events
#'
#' @description Internal function used by \code{mwss} function to generate the event matrix required by SimInf package
#'
#' @usage build_E(compartments, SA)
#'
#' @return Matrix Empty E matrix to store programmed events
#'
#' @keywords internal
#'
#' @noRd

build_E <- function(compartments, SA){

### Build E
### Note: no event for now, entry and exit are handle by transition driven by average length of stay
### will be use to implement transfer of patients between wards
# sample nodes are specified here for events
# Entry
# Exit
# External transfer
# Internal transfer (change of status -vaccinated) but stay in the same ward

ncompartments <- compartments %>% length %>% as.integer

E <- structure(rep(rep(0, ncompartments), 3),
               .Dim = c(ncompartments, 3L),
               .Dimnames = list(compartments,
                                c("1", "2", "3")))
###
### Possible events:
###

## exit patient / vaccination Patient (internal trans)
E[compartments[grepl('PW', compartments)], 1] <- 1

## new patient
# 1- with screening area: entry in one of the PSA compartment
if(isTRUE(SA))
E[compartments[grepl('PSA', compartments)], 2] <- 1 else
# 2- without screening area entry in a PW compartment
E[compartments[grepl('PW', compartments)], 2] <- 1

return(E)
}
