library(ckm)
library(dplyr)

# CKM PARAMETERS ------------------------------

ckm_parameters <- ckm::build_parameters_table(
  Ds = 15,
  Vs = c(6.1, 10, 15, 20.1, 25, 30.1, 35, 40, 45, 50),
  jss = c(4,10)
) |>
  filter(
    !(js == 4 & V > 25),
    !(js == 10 & V < 30)
  )