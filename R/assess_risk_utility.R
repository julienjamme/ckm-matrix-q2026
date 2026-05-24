library(dplyr)
library(purrr)
library(ggplot2)

source("R/functions.R")


# Parameters tested ----------------

source("R/parameters.R")

# Build transition matrices -----------------------  
ckm_matrices <- pmap(ckm_parameters, ckm::create_transition_matrix, .progress=TRUE)

# Tabular data ------------------------------------
source("R/data_preparation.R")

# Compute Risk -----------------------------------

tabulars_risks <- imap(
  tabulars_frequencies,
  \(freq, name){
    map(
      ckm_matrices,
      \(mat){
        D <- mat@pParams@D
        V <- mat@pParams@V
        js <- mat@pParams@js
        freq_threshold <- js+1
        risks <- ckm::assess_risk(mat, freq, I = 1:js, J = 1:freq_threshold) |>
          tail(1) |>
          rename(risk = qij) |>
          mutate(D=D,V=V,js=js,tab=name) |>
          mutate(
            nb_cell_not_null = freq |> summarise(nb = sum(N[i > 0])) |> pull(nb),
            nb_cell_total = freq |> summarise(nb = sum(N)) |> pull(nb),
            sm_count = freq |> summarise(pc = sum(N[i %in% 1:js])/sum(N[i>0])*100) |> pull(pc))
        
        utility <- assess_apriori_utility(freq, mat, precision = 3) |> tail(1) |>
          select(imax=i,U1,L1,L2)
        
        return(
          risks |> bind_cols(utility)
        )
      }
    ) |> list_rbind()
  },
  .progress=TRUE
) |>
  list_rbind()


# Latex output ---------
tabulars_risks |>
  filter(js == 4) |>
  select(tab,V,js,sm_count,risk) |>
  arrange(sm_count) |>
  tidyr::pivot_wider(names_from = V, values_from = risk, names_prefix = "V=") |>
  knitr::kable(format="latex", digits=3,caption="")

tabulars_risks |>
  filter(js == 10) |>
  select(tab,V,js,sm_count,risk) |>
  arrange(sm_count) |>
  tidyr::pivot_wider(names_from = V, values_from = risk, names_prefix = "V=") |>
  knitr::kable(format="latex", digits=3,caption="")


# Graphical representation ------------------

source("R/theme_ggplot.R")

tabulars_risks |>
  filter(js==4) |>
  ggplot() +
  geom_point(
    aes(y = risk, x = sm_count, color = as.factor(V)),
    alpha = 0.65, size = 3
  ) +
  geom_hline(yintercept = 0.6, color = "grey55", linetype = "dotted") +
  geom_hline(yintercept = 0.8, color = "grey55", linetype = "dotted") +
  scale_x_continuous("Small counts > 0 et < 5\n(%)", expand = c(0,0), limits = c(0,100)) +
  scale_y_continuous("Risk", breaks = seq(0,1,0.2), expand = c(0,0), limits = c(0,1)) +
  scale_color_viridis_d("Variance", direction = -1) +
  # scale_color_brewer("Variance", type = "qual", palette = 7) +
  # ggtitle("Risque d'inférence sur les petits comptages en fonction\ndu taux de petits comptages",
  #         subtitle = "En fonction du niveau de variance V et pour D=15 et js=4\nChaque point d'une couleur donnée représente un tableau") +
  theme(legend.position = "inside", legend.position.inside = c(0.75, 0.25))

ggsave("outputs/risk_vs_small_counts_s5.png", device = "png", dpi = 150)


tabulars_risks |>
  filter(js==10) |>
  ggplot() +
  geom_point(
    aes(y = risk, x = sm_count, color = as.factor(V)),
    alpha = 0.65, size = 3
  ) +
  geom_hline(yintercept = 0.6, color = "grey55", linetype = "dotted") +
  geom_hline(yintercept = 0.8, color = "grey55", linetype = "dotted") +
  scale_x_continuous("Small counts > 0 et < 11\n(%)", expand = c(0,0), limits = c(0,100)) +
  scale_y_continuous("Risk", breaks = seq(0,1,0.2), expand = c(0,0), limits = c(0,1)) +
  scale_color_viridis_d("Variance", direction = -1) +
  # scale_color_brewer("Variance", type = "qual", palette = 7) +
  # ggtitle("Risque d'inférence sur les petits comptages en fonction\ndu taux de petits comptages",
  #         subtitle = "En fonction du niveau de variance V et pour D=15 et js=4\nChaque point d'une couleur donnée représente un tableau") +
  theme(legend.position = "inside", legend.position.inside = c(0.75, 0.25))

ggsave("outputs/risk_vs_small_counts_s10.png", device = "png", dpi = 150)







