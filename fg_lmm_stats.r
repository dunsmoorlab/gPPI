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
mod <- lmer(vmPFC_ers~amyg_cem*phase*condition*group + (1|subject),data=df,REML=FALSE)
anova(mod)

#more specific models
cdf = df[which(df$group == 'healthy' & df$phase != 'baseline'),]
pdf = df[which(df$group == 'ptsd' & df$phase != 'baseline'),]

cmod <- lmer(dACC_ers~amyg_cem*phase*condition + (1|subject),data=cdf,REML=FALSE)
anova(cmod)
pmod <- lmer(dACC_ers~amyg_cem*phase*condition + (1|subject),data=pdf,REML=FALSE)
anova(pmod)



#CS+ > CS- differential data 
df1 <- read.csv('ers_subcort_betas_diff.csv')
df2 <- df[which(df$phase == 'extinction'),]
df2$subject <- factor(df2$subject)
df2$vmPFC_ers <- ave(df2$vmPFC_ers,df2$group,FUN=scale)
df2$hc_head <- ave(df2$hc_head,df2$group,FUN=scale)
mod <- lm(vmPFC_ers~hc_head+group,data=df2)#+ (1|subject),data=df2,REML=FALSE)
summary(mod)


df <- read.csv('subcort_betas_lmm.csv')
edf <- df[which(df$phase == 'extinction' & df$group == 'ptsd'),]
amyg_mod <- lmer(dACC_ers~amyg_bla*condition + (1|subject), data=edf, REML=FALSE)
amyg_mod2 <- lmer(dACC_ers~amyg_bla*amyg_cem*condition*group + (1|subject),data=edf,REML=FALSE)
anova(amyg_mod)

summary(amyg_mod)

bla_mod <- lmer(vmPFC_ers~condition*phase*group + (1|subject),data=df,REML=FALSE)
anova(bla_mod)

ldf <- df[which(df$phase != 'baseline'),]
mod <- glmer(factor(phasecon) ~ dACC_ers * group + (1|subject),data=df,family="binomial")
anova(mod)
summary(mod)



#quick connectivity stuff
df <- read.csv('quick-conn.csv')
df <- df[which(df$group == 'healthy'),]
mod <- lmer(conn~target*cope + (1|subject),data=df,REML=FALSE)
anova(mod)
summary(mod)
