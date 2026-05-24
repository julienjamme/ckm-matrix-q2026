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


#' Title
#'
#' @param D1 Deviation for small counts
#' @param V1 Variance for small counts
#' @param js1 js for small counts
#' @param D2 Deviation for last count
#' @param V2 Variance for last count
#' @param js2 js for last count
#'
#' @returns
#' @export
#'
#' @examples
build_stacked_matrix <- function(
    D1, V1, js1,
    D2, V2, js2=0
){
  
  matrice_transition1 <- ptable::create_cnt_ptable(D = D1, V =V1, js=js1)
  matrice_transition2 <- ptable::create_cnt_ptable(D = D2, V =V2, js=js2)
  
  matrice_transition <- matrice_transition1
  
  matrice_transition@pTable <- bind_rows(
    matrice_transition1@pTable |> filter(i != max(i)),
    matrice_transition1@pTable |> filter(i == max(i)) |> select(i, j) |> 
      bind_cols(
        matrice_transition2@pTable |> filter(i == max(i)) |> select(p:type)
      )
  )
  
  matrice_transition@empResults <- matrice_transition@pTable |> 
    group_by(i) |>
    summarise(
      p_mean = sum(v*p),
      p_var = sum(v^2*p) - sum(v*p)^2,
      p_sum = sum(p),
      .groups = "drop"
    ) |>
    left_join(matrice_transition@pTable |> filter(i==j) |> select(i, p_stay = p), by ="i") |>
    mutate(p_stay = ifelse(is.na(p_stay), 0, p_stay)) |>
    data.table::as.data.table()
  
  return(
    list(
      matrix_1 = matrice_transition1,
      matrix_2 = matrice_transition2,
      stacked = matrice_transition
    )
  )
}


