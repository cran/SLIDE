#' Bootstrap Cutoff for SLIDE
#'
#' @param uninfected A dataframe of protein expression levels in an uninfected subset of cells. All columns must be numeric expression levels.
#' @param fraction Fraction of uninfected cells used for the SLIDE bootstrap. This value should reflect the size of the infected cell population relative to the uninfected cell population. Must be between >0 and <=0.2. Default is 0.1.
#' @param iter Iterations of the SLIDE bootstrap. Default is 20.
#' @param level Alpha level to determine cutoff for a significant mean distance ratio. Default is 0.05.
#' @return This function bootstraps a cutoff for the mean distance ratio from the uninfected cell population by running the SLIDE procedure on a random subset of uninfected cells a specified number of times (default 20x). Only to be used if our pre-determined value of 1.2 is insufficient. Depending on how large the uninfected dataset is, this function may take many hours to complete.
#' @examples 
#' cutoff <- bootstrap_cutoff(uninfected = UN_sig, fraction = (nrow(I_sig)/nrow(UN_sig)), iter=1)
#' cutoff
#' @references
#' Sen, N., Mukherjee, G., Sen, A., Bendall, S. C., Sung, P., Nolan, G. P., & Arvin, A. M. (2014). Single-Cell Mass Cytometry Analysis of Human Tonsil T Cell Remodeling by Varicella Zoster Virus. Cell Reports, 8(2), 633–645. https://doi.org/10.1016/j.celrep.2014.06.024
#' @export
bootstrap_cutoff <- function(uninfected, fraction=0.1, iter=20, level=0.05)
{
  # is size between 0 and 1?
  if (fraction > 0.2 | fraction <= 0)
    stop("Fraction of uninfected cells must be >0 and <=0.2")

  nullDist <- numeric(iter)
  size <- fraction*nrow(uninfected)

  for(i in 1:iter) {
    index <- sample(x = c(1:nrow(uninfected)),
                    size = nrow(uninfected),
                    replace = F)
    I <- uninfected[index[1:size],]
    UN <- uninfected[index[size:length(index)],]
    control_results <- slide(uninfected = UN, infected = I)
    nullDist[i] <- control_results$meanRatio
  }
  return(quantile(nullDist, 1-level))
}
