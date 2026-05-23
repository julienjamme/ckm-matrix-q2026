library(dplyr)
library(purrr)

# Download detailed tabular data for French Census 
# And Build agregated data with different level of details
#----------------------------------

temp <- tempfile()
download.file("https://www.insee.fr/fr/statistiques/fichier/8582668/TD_FOR2_2022_csv.zip",temp)
census_tab <- readr::read_delim(
  unz(temp, "TD_FOR2_2022.csv"), 
  delim = ";", 
  col_types = readr::cols(NB = "n", .default="c")
) |>
  mutate(DEP = substr(CODGEO, 1, 2)) |>
  mutate(NB = as.integer(ifelse(is.na(NB),0,round(NB, digits = 0)))) |>
  filter(NB > 0)
unlink(temp)

# Description of the table ----
str(census_tab)

nrow(census_tab) # More than 3 million cells
census_tab |> filter(NB < 1) |> nrow() 
census_tab |> filter(NB == 1) |> nrow() 
quantile(census_tab$NB, probs = seq(0,1,0.1))


# Description of the tabular data to build ----

tabulars_BK <- list(
   # c("DEP","SEXE","AGEQ65", "DIPL_19"),
   # c("CODGEO","SEXE"),
   c("CODGEO","DIPL_19"),
   c("CODGEO","AGEQ65"),
   c("CODGEO","SEXE","AGEQ65"),
   c("CODGEO","SEXE","DIPL_19"),
   c("CODGEO","AGEQ65","DIPL_19"),
   c("CODGEO","SEXE","AGEQ65", "DIPL_19")
)
names(tabulars_BK) <- paste0("tab_", seq_along(tabulars_BK))

tabulars_data <- imap(tabulars_BK, 
    \(cats,tb){
      print(tb)
      res <- census_tab |> 
        group_by(across(all_of(cats))) |>
        summarise(NB = sum(NB), .groups = "drop") |>
        mutate(
          is_small_5 = NB > 0 & NB < 5,
          is_small_11 = NB > 0 & NB < 11
        )
      count(res, is_small_5) |> print()
      count(res, is_small_11) |> print()
      return(res)
    }
)
        
str(tabulars_data)

tabulars_data |> map(\(t) count(t, is_small_5) |> mutate(pc = n/sum(n)*100) |> filter(is_small_5 == TRUE) |> pull(pc))
tabulars_data |> map(\(t) count(t, is_small_11) |> mutate(pc = n/sum(n)*100) |> filter(is_small_11 == TRUE) |> pull(pc))


# Compute frequencies distribution --------------
tabulars_frequencies <- tabulars_data |>
  imap(
    \(tab,name){
    tab |> rename(nb_obs = NB) |> ckm::compute_frequencies(cat_vars = tabulars_BK[[name]])
  })


# purrr::iwalk(
#   list(freq_census_tab1 = freq_cens_tab, freq_census_tab2 = freq_cens_tab2),  
#   \(d,n) write.csv(d, file = file.path("data", paste0(n,".csv")), row.names = FALSE)
# )