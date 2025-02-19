


library(HVdata)
library(HVtools)
library(xlsx)
library(dplyr)
library(stringr)



# Data --------------------------------------------------------------------

# AGES phenotypes
pheno <-  data.table::fread("V:/Somapanel_GNF/Rdata/RANDIDS_AGES_GNF_PHENOTYPES.csv")

pheno = pheno %>% mutate(CORNCAL1_mod = log(pheno$CORNCAL1 + 1))
pheno = pheno %>% mutate(TG_mod = log(pheno$S_TG))
pheno = pheno %>%
  mutate(S_LDL_adjmed = ifelse(STATINS==1, S_LDL/0.7, S_LDL),
         S_CHOL_adjmed = ifelse(STATINS==1, S_CHOL/0.8, S_LDL),
         SYS_adjmed = ifelse(HTNMED==1, SYS+15, SYS),
         DID_adjmed = ifelse(HTNMED==1, DID+10, DID),
         PP_adjmed = SYS_adjmed - DID_adjmed)

pheno = pheno %>% mutate(CORNCAL1_mod = log(pheno$CORNCAL1 + 1))

# 7k soma data
xpsct = HVdata::readSoma7k("boxcox","excluded", "a1", humanonly = F, exclude_sample_outliers = T)  %>%
  mutate(S_ID = as.numeric(S_ID))

# protein annotation
protann = hvtools::read_protein_ids()

# join phenotypes and protein data in one data frame
data1 = pheno %>% inner_join(xpsct)


outcome = "CORNCAL1_mod"

# subset data on non-missing for outcome, other criteria could be added
subdata = data1 %>%
  filter(!is.na(!!rlang::parse_expr(outcome)))

# fit a base model with the outcome and covariates of interest
fit0 <- lm(as.formula(paste0(outcome, " ~ AGE + SEX")), data=subdata) 

# Perform protein associations
for (p in 3:ncol(xpsct))
{ prot = names(xpsct)[p]  
model <- update(fit0, as.formula(paste0(". ~ .+ ",prot)), data=subdata)

newline = data.frame(Somamer = prot,
                     Beta = coef(summary(model))[prot,1],
                     SE = coef(summary(model))[prot,2],
                     P = coef(summary(model))[prot,4],
                     R2_adj = summary(model)$adj.r.squared)

if (exists("res")){
  res = bind_rows(res, newline)
} else {
  res = newline
}
}

res$FDR <- p.adjust(res$P,method = "fdr")
res$Bonf <- p.adjust(res$P,method = "bonferroni")

totallist = res %>%
  left_join(protann, by=c("Somamer"="SEQ_ID")) %>%
  select(Somamer, Target, EntrezGeneSymbol, Beta, SE,  P, FDR, Bonf, R2_adj, everything()) %>%
  arrange(P)

View(totallist)
library(xlsx)

write.xlsx(totallist, file=paste0("linear_", outcome,".xlsx"), row.names=F)



