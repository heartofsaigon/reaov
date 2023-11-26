## code to prepare `DATASET` dataset goes here

if(!library(faraway, logical.return = T)) install.packages("faraway")
if(!library(tidyverse, logical.return = T)) install.packages("faraway")

mydata1<- faraway::gala|> as_tibble()|> select(-Endemics)
mydata2<- faraway::odor|> as_tibble()


usethis::use_data(mydata1, overwrite = TRUE)
usethis::use_data(mydata2, overwrite = TRUE)
