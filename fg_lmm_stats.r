require(lme4)
require(lmerTest)
require(emmeans)
require(magrittr)
require(ggplot2)
phases = c('baseline','acquisition','extinction')
cons = c('CS+','CS-')
groups = c('healthy','ptsd')

#behavior example
setwd('/Users/ach3377/Documents/gPPI')
df <- read.csv('pfc_ers_cleaned.csv')
df2 <- df[which(df$phase != 'baseline'),]
ers <- lmer(rsa~phase*condition*roi*group + (1|subject),data=df,REML=FALSE)
anova(ers)
ers.emm <- emmeans(ers,list(pairwise ~ phase:condition:roi:group),adjust="None")
beh.emmc <- broom::tidy(beh.emm$`pairwise differences of encode_phase, condition, response_phase`)

w_group <- lmer(rsa~phase*condition*roi*group + (1|subject),data=df,REML=FALSE)
anova(w_group)
wo_group <- lmer(rsa~phase*condition*roi + (1|subject),data=df,REML=FALSE)
anova(wo_group)
lrtest(wo_group,w_group)

#univariate
betas <- read.csv('subcortical_betas.csv')
bmod <- lmer(beta~phase*condition*roi*group + (1|subject),data=betas,REML=FALSE)
anova(bmod)

#uni to ers
df <- read.csv('ers_subcort_betas_full.csv')
mod <- lmer(dACC_ers~hc_tail*phase*condition*group + (1|subject),data=df,REML=FALSE)
anova(mod)

#more specific models
cdf = df[which(df$group == 'healthy' & df$phase != 'baseline'),]
pdf = df[which(df$group == 'ptsd' & df$phase != 'baseline'),]

cmod <- lmer(dACC_ers~amyg_cem*phase*condition + (1|subject),data=cdf,REML=FALSE)
anova(cmod)
pmod <- lmer(dACC_ers~amyg_cem*phase*condition + (1|subject),data=pdf,REML=FALSE)
anova(pmod)

