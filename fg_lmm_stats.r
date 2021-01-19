require(lme4)
require(lmerTest)
require(emmeans)
require(magrittr)
require(ggplot2)
require(effects)
require(interactions)
require(jtools)
require(dplyr)
require(reghelper)
require(RColorBrewer)
require(dotwhisker)
afex::set_sum_contrasts()

phases = c('baseline','acquisition','extinction')
cons = c('CS+','CS-')
groups = c('healthy','ptsd')
group_pal <- colors()[c(565,555)]
phase_pal <- colors()[c(454,614)]
con_pal <- colors()[c(17,92)]
#setwd('/Users/ach3377/Documents/gPPI')
setwd('C:\\Users\\ACH\\Documents\\gPPI')
df <- read.csv('all_data_lmm.csv')
ers <- df[which(df$phase %in% c('baseline','acquisition','extinction')),]
emo <- df[which(df$phase %in% c('acquisition','extinction')),]
emoz <- read.csv('emo_data_lmm_zscore.csv')

emm_options(lmer.df="satterthwaite")
emm_options(lmerTest.limit = 15000)

#PFC reinstatement stuff----------------------------------------------------------------------
ers <- read.csv('pfc_ers_cleaned_lmm.csv')
ers <- ers[which(ers$phase %in% c('acquisition','extinction')),]
ers <- ers %>% 
  mutate(condition = recode(condition, 
                            "CS-" = "0", 
                            "CS+" = "1"),
         phase = recode(phase,
                        "baseline" = "a",
                        "acquisition" = "b",
                        "extinction" = "c"))
mod <- lmer(ers~phase*condition*roi*group + (1|stimulus) + (1|subject),data=ers,REML=FALSE)
summary(mod)
anova(mod)

condif <- emmeans(mod, revpairwise ~ condition|phase*group*roi, adjust='None', lmer.df="satterthwaite")
condif$contrasts

condif.emm <- emmeans(mod, ~condition*phase*group*roi, adjust='None', lmer.df="satterthwaite")

healthy_fear_dACC <- c(-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
healthy_fear_vmPFC <- c(0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0)

healthy_ext_dACC <- c(0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0)
healthy_ext_vmPFC <- c(0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0)

ptsd_fear_dACC <- c(0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0)
ptsd_fear_vmPFC <- c(0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0)

ptsd_ext_dACC <- c(0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0)
ptsd_ext_vmPFC <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1)


contrast(condif.emm, method=list(healthy_fear_dACC - healthy_fear_vmPFC), adjust="None")
contrast(condif.emm, method=list(healthy_ext_vmPFC - healthy_ext_dACC), adjust="None")

contrast(condif.emm, method=list(ptsd_fear_dACC - ptsd_fear_vmPFC), adjust="None")
contrast(condif.emm, method=list(ptsd_ext_vmPFC - ptsd_ext_dACC), adjust="None")

contrast(condif.emm, method=list(healthy_fear_dACC - ptsd_fear_dACC), adjust="None")
contrast(condif.emm, method=list(healthy_ext_vmPFC - ptsd_ext_vmPFC), adjust="None")


#workspace------------------------------------------------------------------------------------
basic <- lmer(hc_tail_ret_uni~phase*condition*group + (1|subject),data=emo,REML=FALSE)
anova(basic)
summary(basic)

mod <- lmer(pfc_diff_ers~amyg_cem_ret_uni*phase*condition*group + (1|subject),data=emoz,REML=FALSE)
anova(mod)
summary(mod)
modsum <- summary(mod)
#anova(basic,mod)

probe_interaction(mod,"hc_tail_ret_uni",modx="condition",mod2="phase",cond.int=TRUE,interval=TRUE,centered="none",colors=con_pal)
emmip(mod,condition~ev,at=list(ev=seq(0,1,by=.1),condition=c("1","0")),style="factor")
emmeans(basic,pairwise ~ phase*condition,adjust='None')


hc_ret_uni <- lmer(pfc_diff_ers~hc_tail_ret_uni*phase*condition*group + hc_body_ret_uni*phase*condition*group + hc_head_ret_uni*phase*condition*group + phase*condition*group + (1|subject),data=emoz,REML=FALSE)
hc_ers <- lmer(pfc_diff_ers~hc_tail_ers*phase*condition*group + hc_body_ers*phase*condition*group + hc_head_ers*phase*condition*group + phase*condition*group + (1|subject),data=emoz,REML=FALSE)
amyg_ret_uni <- lmer(pfc_diff_ers~amyg_bla_ret_uni*phase*condition*group + amyg_cem_ret_uni*phase*condition*group + phase*condition*group + (1|subject),data=emoz,REML=FALSE)
amyg_ers <- lmer(pfc_diff_ers~amyg_bla_ers*phase*condition*group + amyg_cem_ers*phase*condition*group + phase*condition*group + (1|subject),data=emoz,REML=FALSE)


# gPPI --------------------------------------------------------------------
df <- read.csv('cleaned_gPPI_lmm.csv')
anova(basic)
basic <- lmer(hc_head~phase*condition*group + (1|subject),data=df,REML=FALSE)
summary(basic)

mod <- lmer(pfc_diff_ers~hc_tail*phase*condition*group + (1|subject),data=df,REML=FALSE)
summary(mod)
probe_interaction(mod,"hc_head",modx="phase",cond.int=TRUE,interval=TRUE,centered="none",colors=phase_pal)
