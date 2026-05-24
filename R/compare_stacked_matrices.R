library(ptable)
library(ckm)
library(purrr)
library(dplyr)
library(ggplot2)

source("R/theme_ggplot.R")
source("R/functions.R")
source("R/data_preparation.R")
source("R/parameters.R")


stacked_matrix <- build_stacked_matrix(
  D1=15,V1=35,js1=10,
  D2=15,V2=5,js2=0
)


# Assess risk - Utility

tabulars_risks_stacked <- imap(
  tabulars_frequencies,
  \(freq, name){
    
        js=10
        freq_threshold <- js+1
        risks <- ckm::assess_risk(stacked_matrix$stacked, freq, I = 1:js, J = 1:freq_threshold) |>
          tail(1) |>
          rename(risk = qij) |>
          mutate(tab=name, matrix="stacked")
        
        utility <- assess_apriori_utility(freq, stacked_matrix$stacked, precision = 3) |> tail(1) |>
          select(imax=i,U1,L1,L2)
        
        return(
          risks |> bind_cols(utility)
        )
      }
    ) |> list_rbind()


tabulars_risks_stacked

tabulars_risks |>
  filter(D==15,V==35,js==10) |> mutate(matrix="original") |>
  select(tab, matrix, risk, U1) |>
  bind_rows(
    tabulars_risks_stacked |>
      select(tab, matrix, risk, U1)
  ) |>
  arrange(tab) |>
  knitr::kable(format="latex", digits=3, caption="")

