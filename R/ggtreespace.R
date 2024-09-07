#' @title Plot phylomorphospace
#'
#' @description This function plots a phylomorphospace by mapping a tree
#' object onto a vector space like morphospace.
#' @param tr a tree object. This should be an object of class that is
#'           compatible with `ggtree`, typically an object of class
#'           `phylo` or `treedata`.
#' @param trait Trait data as a data frame or matrix, where each row
#' represents a tree tip or node.
#'
#'     For data matching the number of tips, ancestral traits are reconstructed
#'     for internal nodes.
#'
#'     For data equal to the total number of nodes, values are directly used as
#'     node coordinates.
#' @param mapping aesthetic mapping
#' @param ... additional parameters for customization with `ggtree`. Please
#' use `?ggtree::ggtree` for more information.
#' @return ggtreeSpace object
#' @importFrom ggtree ggtree
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 coord_cartesian
#' @examples
#' library(ggtree)
#' library(phytools)
#' library(ggtreeSpace)
#'
#' tr <- rtree(15)
#' td <- fastBM(tr, nsim = 2)
#' ggtreespace(tr, td) +
#'     geom_tippoint()
#'
#' @export
ggtreespace <- function(tr, trait, mapping = NULL, ...) {
  
    if (is.null(trait)) {
        stop("Traits data is required.")
    }
  
    if (!is.data.frame(trait) && !is.matrix(trait) && !is.vector(trait)) {
      stop("The input trait data must be a data frame, matrix, or vector 
           cotaining trait names.")
    }
  
    if (is.data.frame(trait) || is.matrix(trait)) {
      if (is.null(colnames(trait)) || length(colnames(trait)) == 0) {
        c <- c("x", "y")
      } else {
        c <- colnames(trait)
      }
    
    trd <- make_ts_data(tr, trait)

    p <- ggtree(trd,
        mapping = mapping,
        layout = "equal_angle",
        ...
    ) +
        theme_treespace() +
        labs(
            x = c[1],
            y = c[2]
        )

    suppressMessages(p <- p + coord_cartesian())

    class(p) <- c("ggtreeSpace", class(p))

    return(p)
    } 
  
    if (is.vector(trait)){
      
      if (!inherits(tr, "treedata")){
        stop("Only treedata object supports internal calling.")
      }
      
      if (length(trait) == 1){
        stop("More than one traits is needed.")
      }
      
      if (length(trait) > 2) {
        warning("Only the first 2 trait data names will be used.")
      }
      
      
      p <- ggtree(tr, mapping = mapping, layout = intern_call, 
                  layout.params = list(t = trait, as.graph = FALSE)) +
              theme_treespace() +
              labs(
                  x = trait[1],
                  y = trait[2]
                  )
      
      suppressMessages(p <- p + coord_cartesian())
      
      class(p) <- c("ggtreeSpace", class(p))
      
      return(p)
    } 
}


