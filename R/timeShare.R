#' Built contact matrix based individual schedules of professionals.
#'
#' @description \code{timeShare} returns summaries the individual schedules of professionals at the ward scale to produce the contact matrix: average time spent in each ward (column) by professionals of each ward (row).
#' The function considers that a professional is affiliated to the ward where he spend the most part of its activity. In case of equality, the order of rows is used to affiliate professionals to wards.
#'
#' @usage timeShare(contacts, namesincol1 = TRUE)
#'
#' @param contacts data.frame. Contains the proportion of time spent by each professional (column) in each ward (rows).
#' @param namesincol1 logical. If TRUE, the first column of contacts must contain the ward ids.
#' If FALSE, the ward ids are rownames of the contact data.table. Default is TRUE.
#'
#' @examples
#' ## Two wards: "a" and "b" with three professionals in "a" and two in "b"
#' ward_names <- letters[1:2]
#' pop_size_H <- c(3, 2)
#'
#' contacts <- matrix(data = 0,
#'                    nrow = length(ward_names), ncol = sum(pop_size_H),
#'                    dimnames = list(ward_names,
#'                                    seq(sum(pop_size_H))))
#'
#' ## Individual schedules
#' contacts <- apply(contacts, 2, function(x){
#' # Random proportions
#' rp <- round(runif(1) * 100)
#' c(rp,100-rp)
#' })
#'
#' contacts <- as.data.frame(contacts)
#'
#' rownames(contacts) <- ward_names
#'
#' timeShare(contacts, namesincol1 = FALSE)
#'
#' @return data.frame
#'
#' @export


timeShare <- function(contacts, namesincol1 = TRUE){

  if(isTRUE(namesincol1)){
    rownames(contacts) <- contacts[,1]
    contacts <- contacts[,-1]
  }

  wards <- rownames(contacts)
  HCWSs <- colnames(contacts)
  nwards <- length(wards)



  # Each HCWS is affiliated to the ward where it spend more time
  affiliation <- sapply(HCWSs, function(HCWS)
    contacts[,HCWS] %>% which.max %>% wards[.]
    )

  contactMat <- lapply(wards, function(WD){
    nconnect <- HCWSs[affiliation == WD] %>% length

    if(nconnect == 0)
      return(rep(0, nwards))

    if(nconnect == 1)
    return(contacts[, HCWSs[affiliation == WD]] %>% divide_by(sum(.)) %>% round(.,2))

    if(length(HCWSs[affiliation == WD]) > 1)
    return(contacts[, HCWSs[affiliation == WD]] %>% rowMeans %>% round(., 2))

  }) %>% do.call(rbind, .)

  rownames(contactMat) <- colnames(contactMat) <- wards

  return(contactMat)
}
