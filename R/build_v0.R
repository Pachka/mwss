#' Build and initialize v0
#'
#' @description Internal function used by \code{mwss} function to generate the initial v0 matrix required by SimInf package.
#'  During simulation weighted proportion of infected in specific subpopulations will be stored in v0 at the end of each time step.
#'
#' @usage build_v0(u0, SA)
#'
#' @return Matrix. Empty v0 matrix.
#'
#' @keywords internal
#' @noRd

build_v0 <- function(u0, SA){

  if(isTRUE(SA))
    cols <- c("wpropInfHorig", "wpropInfHdest", "wpropInfPSAdest", "wpropInfPWdest",
              # admission counter per ward
              rbind(paste0("nadm", seq(nrow(u0))),
                    paste0("nadmInf", seq(nrow(u0))))) else
                      cols <- c("wpropInfHorig", "wpropInfHdest", "wpropInfPWdest")

    v0 <- matrix(data = 0,
                 nrow = nrow(u0),
                 ncol = length(cols),
                 dimnames = list(NULL, cols))

  v0 %<>% data.frame

  return(v0)
}
