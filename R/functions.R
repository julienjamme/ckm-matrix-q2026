#' Compute several utility measures 
#'
#' @param freqs matrix of frequencies of the counts in the tabular data
#' @param transition_matrix transition matrix from ptable package or ckm package
#' @param precision integer, max loss to measure the utility
#'
#' @returns data.frame : i, pi, U1, L1, L2
#' 
#' @details
#' U1_i = P(|Z| <= precision | X = i), où d est la précision souhaitée
#' L1_i = E(|Z| | X = i)
#' L2 = E(|Z| | X > 0) = Somme_{i>0} P(X=i) * L1_i # P(X=i) estimée par les frequences
#' empiriques (hors 0)
#' 
#' U1, L1 are computed conditionally to i
#' L2 is computed for all the tabular data
#'  
#' @export
#'
#' @examples
assess_apriori_utility <- function(freqs, transition_matrix, precision = 3){
  
  require(dplyr)
  
  i_max <- max(transition_matrix@pTable$i)
  
  freqs |> 
    mutate(i = ifelse(i >= i_max, i_max, i)) |>
    filter(i > 0) |>
    group_by(i) |>
    summarise(pi = sum(p_hat), .groups = "drop") |>
    full_join(
      transition_matrix@pTable |> filter(i > 0),
      by = "i"
    ) |> 
    group_by(i, pi) |> 
    summarise(
      U1 = sum(p[abs(v) <= precision]),
      L1 = sum(abs(j-i)*p),
      .groups="drop"
    ) |>
    mutate(
      L2 = sum(pi * L1)
    )
  
}