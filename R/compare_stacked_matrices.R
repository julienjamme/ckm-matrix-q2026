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


stacked_matrix$stacked@empResults |>
  knitr::kable(format="latex",digits=3,caption="")

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

# Empirical perturbation CKM-like

set.seed(40889)

tabulars_data_noisy <- tabulars_data |>
  map(
    \(tab_data){
      N <- nrow(tab_data)
      keys <- runif(N)
      original <- tab_data |> mutate(CK = keys) |>
        apply_ckm(D=15,V=35,js=10, ck_var = "CK", cnt_var = "NB") 
      original <- original$tab |> rename(NB_mat_o = NB_ckm)
      
      stacked <- tab_data |> mutate(CK = keys) |>
        apply_ckm_with_stacked_mat(stacked_matrix$stacked) |>
        rename(NB_mat_s = NB_ckm)
      
      return(full_join(original, stacked))
    },.progress = TRUE
  )


# Outputs --------------------

tabulars_risks |>
  filter(D==15,V==35,js==10) |> mutate(matrix="original") |>
  select(tab, matrix, risk, U1) |>
  bind_rows(
    tabulars_risks_stacked |>
      select(tab, matrix, risk, U1)
  ) |>
  arrange(tab) |> 
  full_join(
    tabulars_data_noisy |> 
      imap(\(tab, name){
        tab |> 
          summarise(
            original = mean(abs(NB - NB_mat_o)),
            stacked = mean(abs(NB - NB_mat_s))
          ) |>
          tidyr::pivot_longer(c(original, stacked), names_to = "matrix", values_to = "MAE_all") |>
          full_join(
            tab |> 
              summarise(
                original = mean(abs(NB[NB>10] - NB_mat_o[NB>10])),
                stacked = mean(abs(NB[NB>10] - NB_mat_s[NB>10]))
              ) |>
              mutate(tab = name) |>
              tidyr::pivot_longer(c(original, stacked), names_to = "matrix", values_to = "MAE_11p")
          )
        
      }
      ) |> 
      list_rbind(), 
    by = c("tab", "matrix")) |>
  knitr::kable(format="latex", digits=3, caption="")


# empirical check of bias -------------

tabulars_data_noisy$tab_6 |> 
  group_by(NB) |>
  summarise(mean_o = mean(NB_mat_o), mean_s = mean(NB_mat_s), .groups="drop") |>
  mutate(bias_o = mean_o - NB, bias_s = mean_s - NB) |>
  View()





