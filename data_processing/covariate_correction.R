library(dplyr)
library(data.table)

#import script parameters from config
source("scripts/config/data_processing/covariate_correction_config.R")

#import protein data in wide and remove Nan
targets <- fread(targets_fn, data.table=FALSE)
targets_na <- na.omit(targets)

#import covariates
cov <- fread("data/proteins/AGES_agesex.txt", data.table=FALSE)
targets_cov <- left_join(targets_na, cov, by='S_ID')

#remove labels
targets_filt <- subset(targets_cov, select = -c(S_ID, PROJ, AGE, SEX))

#fit linear model with covariates
dat_adj <- function(df) {
	fit <- lm(df ~ AGE + SEX, data = targets_cov)
	return(resid(fit))
}

#obtain model residuals
adj_res <- lapply(targets_filt, dat_adj)
df_adj <- data.frame(adj_res)
cols <- subset(targets_na, select = c(S_ID, PROJ))
out <- cbind(cols, df_adj)

#export adjusted proteins
fwrite(out, out_fn, row.names=FALSE, sep="\t", quote = FALSE)