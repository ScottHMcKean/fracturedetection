#' Fixed Point Function
#'
#' @description This function runs a recursive fixed point approximation
#' to expand the strong edges to weak edges. Not the most efficient function
#' in the world.
#'
#' @param ws list of weak and strong edges
#' @return ws list of weak and strong edges
expandStrong <- function(ws) {
  overlap <- grow(ws$strong,2) & ws$weak
  ws$strong[overlap] <- TRUE
  ws$weak[overlap] <- FALSE
  ws
}
