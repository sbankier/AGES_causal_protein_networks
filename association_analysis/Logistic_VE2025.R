library(HVdata)
library(hvtools)
library(xlsx)
if(!require(rms)){
  install.packages("rms")
  library(rms)
} 



# Data --------------------------------------------------------------------

# AGES phenotypes
pheno <-  data.table::fread("V:/Somapanel_GNF/Rdata/RANDIDS_AGES_GNF_PHENOTYPES.csv")
#?readSoma7k

# 7k soma data
xpsct = readSoma7k("boxcox","excluded", "a1", humanonly = T, exclude_sample_outliers = T)  %>%
  mutate(S_ID = as.numeric(S_ID))

# protein annotation
protann = read_protein_ids()

# Define outcome specifically
pheno = pheno %>%
  mutate(prev3 = pheno$prev1 & pheno$prev2)

# join phenotypes and protein data in one data frame
data1 = pheno %>% inner_join(xpsct)

# Define base model for both sexes -------------------------------------------------------

# define outcome
outcome = "prev3"

# subset data on non-missing for outcome, other criteria could be added
subdata = data1 %>%
  filter(!is.na(!!rlang::parse_expr(outcome)))

# fit a base model with the outcome and covariates of interest
fit0 <- lrm(as.formula(paste0(outcome, " ~ AGE + SEX + EGFR")), data=subdata) 

X0 <- model.matrix(fit0) 
xpindex <- ncol(X0)+1

# Define base model for sexes separately -------------------------------------------------------
subdata$SEX <- factor(subdata$SEX, levels = c(1, 2), labels = c("Male", "Female"))
# Define the outcome variable
outcome <- "prev3"
# Fit logistic regression models separately for each sex. SMOKINGSTATUS + BMI + ALCOHOLGWEEK
fits_by_sex <- list()
sexes <- levels(subdata$SEX)

for (sex in sexes) {
  subdata_sex <- subset(subdata, SEX == "Male")
  fit0 <- lrm(as.formula(paste0(outcome, " ~ AGE")), data = subdata_sex)
}
X0 <- model.matrix(fit0) 
xpindex <- ncol(X0)+1

# Protein associations ----------------------------------------------------


# create data frame to collect association results
resframe <- data.frame(matrix(data = NA, ncol = 7, nrow = ncol(xpsct))) 
names(resframe) <- c("Somamer", "P", "R2.adj", "AUC","Beta","Brier","SE")
names(resframe) <- c("Somamer", "Beta","SE", "P",  "R2.adj", "AUC","Brier")

# loop over all proteins and test each for association with outcome
for (i in 3:ncol(xpsct)){
  cat (paste0(i-2, "..."))
  prot = names(xpsct)[i]
  fiti <- update(fit0, as.formula(paste0(". ~ .+ ",prot)), data=subdata)
  resframe[i-2,1] <- names(xpsct)[i]
  resframe[i-2,2] <- fiti$coefficients[xpindex]
  resframe[i-2,3] <- sqrt(diag(vcov(fiti)))[xpindex]
  resframe[i-2,4] <- 2*pnorm(abs(fiti$coefficients[xpindex] / sqrt(diag(fiti$var))[xpindex]),lower.tail = F)
  resframe[i-2,5] <- fiti$stats[10]
  resframe[i-2,6] <- fiti$stats[6]
  resframe[i-2,7] <- fiti$stats[15]
}
resframe$OR <- exp(resframe$Beta)
resframe$P.adj_bonf <- p.adjust(resframe$P, method = "bonferroni")
resframe$P.adj_BH <- p.adjust(resframe$P, method = "BH")

totallist = resframe %>%
  left_join(protann, by=c("Somamer"="SEQ_ID")) %>%
  select(Somamer, Target, EntrezGeneSymbol, Beta, SE, OR,  P, P.adj_BH, P.adj_BH, everything()) %>%
  arrange(P)

View(totallist)
# save output as csv or excel
# write.csv(totallist, file=paste0("logistic_", outcome,".csv"), row.names=F)
write.xlsx(totallist, file=paste0("logistic_", outcome,".xlsx"), row.names=F)
