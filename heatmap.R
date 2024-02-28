library("dplyr")
library("gplots")
library("readr")
library("tidyr")
library("ggplot2")

# Ancestral Lineage

covar_al_df <- read_csv("/Users/cindy/projects/covariability/doublesResult_AncestralLineage.csv",
                       col_types=cols(posI=col_factor(), posJ=col_factor()))

covar_al_mpbon_df <- covar_al_df %>% select(posI, posJ, Mp.bon)

options(repr.plot.width=10, repr.plot.height=10)

covar_al_mpbon_df$mpbon_cats <- cut(covar_al_mpbon_df$Mp.bon, breaks=c(.001, .002, .05, .5, 1))

ggplot(data=covar_al_mpbon_df,
       mapping=aes(x=posJ, y=posI)) +
  geom_tile(aes(fill=mpbon_cats), color="white") +
  theme_grey(base_size=8) +
  labs(x="", y="", title="Mp.bon - Ancestral Lineage") +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(drop=FALSE, values=colorRampPalette(c("navy","mediumorchid2"))(4)) +
  guides(fill=guide_legend(title="Mp.bon"))

# Lambda

covar_lambda_df <- read_csv("/Users/cindy/projects/covariability/doublesResult_Lambda.csv",
                            col_types=cols(posI=col_factor(), posJ=col_factor()))

covar_lambda_mpbon_df <- covar_lambda_df %>% select(posI, posJ, Mp.bon)

options(repr.plot.width=10, repr.plot.height=10)

covar_lambda_mpbon_df$mpbon_cats <- cut(covar_lambda_mpbon_df$Mp.bon, breaks=c(.001, .002, .05, .5, 1))

ggplot(data=covar_lambda_mpbon_df,
       mapping=aes(x=posJ, y=posI)) +
  geom_tile(aes(fill=mpbon_cats), color="white") +
  theme_grey(base_size=8) +
  labs(x="", y="", title="Mp.bon - Lambda Variant") +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(drop=FALSE, values=colorRampPalette(c("green4","lightpink"))(4)) +
  guides(fill=guide_legend(title="Mp.bon"))

