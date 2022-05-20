#' Build the ward network
#'
#' @description Internal function used by the \code{ plot_connectivity} function.
#' \code{wardsNet} returns an network object.
#' Ward weights are defined by the size argument.
#' Edges between wards depends on values in the connections matrix.
#'
#' @usage wardsNet(matContact, size)
#'
#' @param connections Matrix. Define the proportion of time spent by professionals in the different wards.
#' @param size Numerical vector. Define the weight attribute of the wards/nodes. Using population or subpopulation sizes is recommended.
#'
#' @importFrom network set.vertex.attribute
#' @importFrom network as.network
#'
#' @return A network object
#'
#' @keywords internal
#' @noRd

wardsNet <- function(connections, size){

  if(sum(connections$nHCWS != 0) > 0){

    connections %<>% .[.$nHCWS!=0, ]

    g <- graph_from_data_frame(connections, directed = TRUE, vertices = data.frame(name = names(size), size = size))

    E(g)$weight <- connections$nHCWS %<>% divide_by(max(.))

    g %<>% intergraph::asNetwork(.)

  } else {
    num_nodes <- length(size)
    my_sociomatrix <- matrix(rep(0, num_nodes * num_nodes),
                             # edge values
                             nrow = num_nodes,
                             #nrow must be same as ncol
                             ncol = num_nodes)

    # diag(my_sociomatrix) <- 0

    g <- as.network(
      x = my_sociomatrix,
      # the network object
      directed = TRUE,
      # specify whether the network is directed
      loops = FALSE,
      # do we allow self ties (should not allow them)
      matrix.type = "adjacency" # the type of input
    )

    set.vertex.attribute(g, "vertex.names", names(size))

  }

  return(g)
}
