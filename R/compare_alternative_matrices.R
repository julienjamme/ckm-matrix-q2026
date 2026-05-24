library(ptable)
library(ckm)
library(purrr)
library(dplyr)
library(ggplot2)

source("R/theme_ggplot.R")

# Example 

ptable::create_cnt_ptable(D=5,V=10,js=0,optim = 1)@empResults |> 
  select(i,p_mean,p_var,p_sum) |>
  knitr::kable(format="latex",digits=4,caption="original")

map(
  0:5,\(js){
    ptable::create_cnt_ptable(D=15,V=10,js=js,optim = 1)@empResults |> 
      filter(i %in% 1:2) |>
      select(i,p_var) |>
      mutate(js=js)
  }
) |>
  list_rbind() |>
  tidyr::pivot_wider(names_from = i, values_from = p_var, names_prefix = "i=") |>
  knitr::kable(format="latex",caption="",digits=4)


ptable::create_cnt_ptable(D=15,V=10,js=0,optim = 4)@empResults |> 
  select(i,p_mean,p_var,p_sum) |>
  knitr::kable(format="latex",digits=4,caption="optim 4") #error

ptable::create_cnt_ptable(D=15,V=10,js=0,optim = 5)@empResults |> 
  select(i,p_mean,p_var,p_sum) |>
  knitr::kable(format="latex",digits=4,caption="alternative")


mat_optim1 <- ptable::create_cnt_ptable(D=15,V=10,js=0,optim = 1)@pTable
mat_optim5 <- ptable::create_cnt_ptable(D=15,V=10,js=0,optim = 5)@pTable


mat_optim1 |> filter(i==1) |> mutate(type="original") |>
  bind_rows(
    mat_optim5 |> filter(i==1) |> mutate(type="biased")
  ) |>
  ggplot() +
  geom_bar(aes(x=v,y=p,fill=type), stat="identity", position="dodge") +
  scale_x_continuous("deviation", breaks=-15:15, expand = c(0,0)) +
  scale_y_continuous("Probability", breaks = seq(0,1,0.1), expand = c(0,0)) +
  scale_fill_brewer("", palette = 2, type = "qual") +
  theme(legend.position = "inside", legend.position.inside = c(0.7,0.7))
  # ggtitle("Distributions of P(X'|X=1) in the biased and unbiased versions")

ggsave("outputs/distribution_compared_i1.png", device = "png", dpi = 150)

mat_optim1 |> filter(i==max(i)) |> mutate(type="original") |>
  bind_rows(
    mat_optim5 |> filter(i==max(i)) |> mutate(type="biased")
  ) |>
  ggplot() +
  geom_bar(aes(x=v,y=p,fill=type), stat="identity", position="dodge") +
  scale_x_continuous("deviation", breaks=-15:15, expand = c(0,0)) +
  scale_y_continuous("Probability", breaks = seq(0,1,0.1), expand = c(0,0)) +
  scale_fill_brewer("", palette = 2, type = "qual") +
  theme(legend.position = "inside", legend.position.inside = c(0.7,0.7))

ggsave("outputs/distribution_compared_imax.png", device = "png", dpi = 150)


# Parameters tested ----------------

ckm_parameters <- ckm::build_parameters_table(
  Ds = 15,
  Vs = c(6.1, 10, 12.1, 15, 20.1, 25, 30.1),
  jss = c(4,6,8,10)
)|>
  filter(
    !(js == 6 & V < 12),
    !(js == 8 & V < 20),
    !(js == 10 & V < 30)
  )

# Build transition matrices -----------------------  

# ckm_matrices <- pmap(ckm_parameters, ptable::create_cnt_ptable, optim = 1, .progress=TRUE)
ckm_bench_matrix <- ptable::create_cnt_ptable(D = 15, V =30.1, js=10, optim = 1)
ckm_alter_matrices <- pmap(ckm_parameters, ptable::create_cnt_ptable, optim = 5, .progress=TRUE)

# Frequencies ----------------------------

source("R/data_preparation.R") #=> to get tabulars_frequencies

uniform_frequencies <- data.frame(
  i = 0:100,
  N = 10,
  p_hat = 10/(10*101)
)

# Assess risk with alternative biased version -------------------------

# original_risks <- 
alternative_risks <- map(
  
  ckm_alter_matrices, \(mat){
    D <- mat@pParams@D
    V <- mat@pParams@V
    js <- mat@pParams@js
    bias_mean <- mean(mat@empResults$p_mean)
    bias_max <- max(mat@empResults$p_mean)
    bias_big_counts <- mat@empResults |> tail(1) |> pull(p_mean)
    
    tabular3_freqs <- assess_risk(mat, tabulars_frequencies$tab_3, I=1:10,J=1:11) |>
      tail(1) |>
      rename(risk = qij) |>
      mutate(tab="tab3") |>
      bind_cols(
        assess_apriori_utility(tabulars_frequencies$tab_3, mat, precision = 3) |> 
          tail(1) |>
          select(imax=i,U1,L1,L2)
      )
    
    uniform_freqs <- assess_risk(mat, uniform_frequencies, I=1:10,J=1:11) |>
      tail(1) |>
      rename(risk = qij) |>
      mutate(tab="uniform")  |>
      bind_cols(
        assess_apriori_utility(uniform_frequencies, mat, precision = 3) |> 
          tail(1) |>
          select(imax=i,U1,L1,L2)
      )
    
    return(
      tabular3_freqs |> 
        bind_rows(uniform_freqs) |>
        mutate(D=D,V=V,js=js,
               bias_mean=bias_mean,
               bias_max=bias_max,
               bias_big_counts=bias_big_counts)
    )
    
  },.progress = TRUE
) |>
  list_rbind()



benchmark_risk <- bind_rows(
  tabular3_freqs <- assess_risk(ckm_bench_matrix, tabulars_frequencies$tab_3, I=1:10,J=1:11) |>
    tail(1) |>
    rename(risk = qij) |>
    mutate(tab="tab3") |>
    bind_cols(
      assess_apriori_utility(tabulars_frequencies$tab_3, ckm_bench_matrix, precision = 3) |> 
        tail(1) |>
        select(imax=i,U1,L1,L2)
    ),
  uniform_freqs <- assess_risk(ckm_bench_matrix, uniform_frequencies, I=1:10,J=1:11) |>
    tail(1) |>
    rename(risk = qij) |>
    mutate(tab="uniform") |>
    bind_cols(
      assess_apriori_utility(uniform_frequencies, ckm_bench_matrix, precision = 3) |> 
        tail(1) |>
        select(imax=i,U1,L1,L2)
    )
) |>
  mutate(
    D=ckm_bench_matrix@pParams@D,
    V=ckm_bench_matrix@pParams@V,
    js=ckm_bench_matrix@pParams@js,
    bias_mean = mean(ckm_bench_matrix@empResults$p_mean),
    bias_max = max(ckm_bench_matrix@empResults$p_mean),
    bias_big_counts = ckm_bench_matrix@empResults |> tail(1) |> pull(p_mean)
  )

alternative_risks |> bind_rows(benchmark_risk) |>
  filter(tab == "tab3") |>
  select(D,V,js,starts_with("bias"),risk, U1) |>
  arrange(risk) |>
  knitr::kable(format="latex",digits=3,caption="")

# alternative_risks |> bind_rows(benchmark_risk) |>
#   filter(tab == "uniform") |>
#   select(D,V,js,starts_with("bias"),risk, U1) |>
#   arrange(risk)

