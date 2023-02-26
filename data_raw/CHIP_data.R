library(dplyr)
library(tidyr)
chip_data <- haven::read_dta('./data_raw/CHIP2018.dta')
df <- chip_data |> select(
    urban,
    wage,
    edu_year,
    experience,
    gender,
    married,
    minzu,
    hukou_local,
    # hukou_ag, # too many NAs
    employ,
    public
) |> filter(
    !is.na(wage) & wage > 0,
    employ == 1
) |>  mutate(
    log_wage = log(wage), .after = wage,
) |> mutate(
    edu_group = as.numeric(edu_year > 12), .before = log_wage # = 1 if postsecondary
) |> mutate (
    exper = experience,
    exper2 = experience^2, .after = edu_year
) |> drop_na() |> select(
    -c(wage, experience, employ)
)

chip_data <- df |> filter(
    urban == 1 # Urban
) |> select(
    -c(urban)
)
usethis::use_data(chip_data, overwrite = TRUE)
