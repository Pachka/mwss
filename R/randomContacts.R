#' Generate random matrix of contact
#'
#' @description \code{randomContacts} returns a list with random contact between wards.
#'
#' @param pop_size_H Numerical vector. Contains the number of professional affiliated to each ward in the same order than ward_names argument.
#' @param ward_names String vector. Contains the ward names in the same order than pop_size_H argument.
#'
#' @importFrom stats runif
#'
#' @return List of two data.frame.
#' sharedTime: proportion of time spent in each ward (row) for each professional (column) ;
#' contactMat: average time spent in each ward (column) by professionals of each ward (row).
#'
#' @examples
#' data("toydata")
#' list2env(toydata,envir=.GlobalEnv)
#'
#' randomContacts(pop_size_H, ward_names)
#'
#' @export

randomContacts <- function(pop_size_H, ward_names){

  nward <- length(pop_size_H)
  nHCWS <- sum(pop_size_H)

  contacts <- matrix(data = 0,
                     nrow = nward, ncol = nHCWS,
                     dimnames = list(ward_names,
                                     seq(nHCWS)))

  for(i in colnames(contacts)){

    pTime <- round((round(runif(nward),2) * 100), digits = -1)
    pTime %<>% .[cumsum(.)<= 100]
    if(length(pTime) < nward)
      pTime[(length(pTime)+1)] <- 100 - sum(pTime)
    if(length(pTime) < nward)
      pTime[(length(pTime)+1) : nward] <- 0

    from <- cumsum(c(1, pop_size_H)) %>% .[-length(.)]
    to <- cumsum(pop_size_H)

    attachedward <- which(as.numeric(i) >= from & as.numeric(i) <= to)

    contacts[attachedward, i] <- max(pTime)
    contacts[-attachedward, i] <- sample(x = pTime[-which.max(pTime)], (nward-1))

  }

  contacts %<>% data.frame

  contactMat <- timeShare(contacts, namesincol1 = FALSE)

  list(sharedTime = contacts,
       contactMat = contactMat)

}
