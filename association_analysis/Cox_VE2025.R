
devtools::install("V:/Somapanel_GNF/Brynjolfur/hvtools")
library(hvtools)
library(usethis)
library(devtools)
library(dplyr)
library(tidyverse)
library(tidymodels)
library(xlsx)

#install 7K protein data ("V:/Somapanel_GNF/HVdata/")
library(HVdata)

#Phenotype data
pheno <-  data.table::fread("V:/Somapanel_GNF/Rdata/RANDIDS_AGES_GNF_PHENOTYPES.csv")

#pheno %>% select(starts_with("MI"))

#Protein data
proteins = readSoma7k("boxcox","excluded", "a1", humanonly = T, exclude_sample_outliers = T)  %>%
  mutate(S_ID = as.numeric(S_ID))

pheno_data = pheno %>%
  mutate(CHDEVENTA1 = pheno$CHDEVENTA & pheno$CHDEVENTB==0,
         CHFREGEVENTA1 = pheno$CHFREGEVENTA & pheno$CHFREGEVENTB==0,
         MIEVENTA1 = pheno$MIEVENTA & pheno$MIEVENTB==0,
         STROKEEVENTA1 = pheno$STROKEEVENTA & pheno$STROKEEVENTB==0)


pheno = pheno %>% inner_join(proteins)

pheno_data$y <- (pheno_data$MIEVENTA1)

protein_data<- proteins #|>


fmla <- survival::Surv(MIFU, y) ~ AGE + SEX

results <- cox_single_point(fmla, pheno_data, protein_data = protein_data)

results_with_names <- results |>
  inner_join(
    read_protein_ids(),
    by = "SEQ_ID"
  ) |>
  select(SEQ_ID, 'Target',  EntrezGeneSymbol,  everything()) |>
  arrange(p.value_adj)

results_with_names |>
  head(100)

View(results_with_names)


write.xlsx(results_with_names, file = "name.xlsx")

