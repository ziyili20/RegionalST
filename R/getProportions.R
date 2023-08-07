#' Define an accessor method for Proportion_CARD
#'
#' @param card A CARD object. 
#'
#' @return A matrix containing the spot-level cell type proportion information
#' @export
#' @examples
#' # getProportion(card)
#'

getProportion <- function(card) {
    
    if (!inherits(card, "CARD")) { # you should replace with the appropriate object
        stop("Input must be a CARD object")
    }
    return(card@Proportion_CARD)
}