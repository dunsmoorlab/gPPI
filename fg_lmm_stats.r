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
require(afex)
afex::set_sum_contrasts()

phases <- c('baseline','acquisition','extinction')
emo_phases <- c('acquisition','extinction')
cons <- c('CS+','CS-')
groups <- c('healthy','ptsd')
group_pal <- colors()[c(565,555)]
phase_pal <- colors()[c(454,614)]
con_pal <- colors()[c(17,92)]
#setwd('/Users/ach3377/Documents/gPPI')
setwd('C:\\Users\\ACH\\Documents\\gPPI')
#df <- read.csv('all_data_lmm.csv')
#ers <- df[which(df$phase %in% c('baseline','acquisition','extinction')),]
#emo <- df[which(df$phase %in% c('acquisition','extinction')),]
#emoz <- read.csv('emo_data_lmm_zscore.csv')

emm_options(lmer.df="satterthwaite")
emm_options(lmerTest.limit = 500)
#PFC encoding-retrieval similarity----------------------------------------------------------------------
pfc <- read.csv('pfc_ers_cleaned_lmm.csv')
#i do this so the phases are in order alphabetically (baseline, conditioning, extinction)
pfc <- pfc %>% mutate(phase = recode(phase, "acquisition" = "conditioning"))

pfc.mod <- mixed(ers ~ condition*phase*roi*group + (1|subject), data=pfc, REML=FALSE, method="LRT")#, test_intercept=TRUE)
summary(pfc.mod)
anova(pfc.mod)
  
  #planned comparisons of CS+ vs. CS-
pfc.csdif <- emmeans(pfc.mod, revpairwise ~ condition|phase*roi*group, lmer.df="asymptotic")
confint(pfc.csdif$contrasts)
summary(pfc.csdif$contrasts)
p.adjust(summary(pfc.csdif$contrasts)$p.value, method="fdr") <.05

#post-hoc between roi and between group comparisons
pfc.csmean <- emmeans(pfc.mod, ~condition*phase*roi*group, adjust='None', lmer.df="asymptotic")

#these are manual codes that set up the double subtraction we want
#position in the vectors corresponds to row in the pfc.csmean object
healthy_acq_dACC  <- c(0,0,-1,1,0,0 ,0,0,0,0,0,0 ,0,0,0,0,0,0, 0,0,0,0,0,0)
healthy_ext_dACC  <- c(0,0,0,0,-1,1 ,0,0,0,0,0,0 ,0,0,0,0,0,0, 0,0,0,0,0,0)
healthy_acq_vmPFC <- c(0,0,0,0,0,0 ,0,0,-1,1,0,0 ,0,0,0,0,0,0, 0,0,0,0,0,0)
healthy_ext_vmPFC <- c(0,0,0,0,0,0 ,0,0,0,0,-1,1 ,0,0,0,0,0,0, 0,0,0,0,0,0)

ptsd_acq_dACC     <- c(0,0,0,0,0,0 ,0,0,0,0,0,0 ,0,0,-1,1,0,0, 0,0,0,0,0,0)
ptsd_ext_dACC     <- c(0,0,0,0,0,0 ,0,0,0,0,0,0 ,0,0,0,0,-1,1, 0,0,0,0,0,0)
ptsd_acq_vmPFC    <- c(0,0,0,0,0,0 ,0,0,0,0,0,0 ,0,0,0,0,0,0, 0,0,-1,1,0,0)
ptsd_ext_vmPFC    <- c(0,0,0,0,0,0 ,0,0,0,0,0,0 ,0,0,0,0,0,0, 0,0,0,0,-1,1)

group.con.list <- list(healthy_acq_dACC - healthy_acq_vmPFC,
                       healthy_ext_vmPFC - healthy_ext_dACC,
                       ptsd_acq_dACC - ptsd_acq_vmPFC,
                       ptsd_ext_vmPFC - ptsd_ext_dACC,
                       healthy_acq_dACC - ptsd_acq_dACC,
                       healthy_ext_vmPFC - ptsd_ext_vmPFC)

group.cons <- contrast(pfc.csmean, method=group.con.list, adjust="None")
confint(group.cons)
summary(group.cons)
p.adjust(summary(group.cons)$p.value, method="fdr")


#need to fdr these somehow

subcort <- read.csv('subcort_ers_cleaned_lmm.csv')
#amygdala------------------------------------------------------------------------------------
amyg <- subcort[which(subcort$roi %in% c('amyg_cem_ers','amyg_bla_ers')),]

amyg.mod <- lmer(ers ~ condition*phase*roi*group + (1|subject), data=amyg, REML=FALSE)
summary(amyg.mod)

#analyze preconditioning (baseline) first
amyg.pre <- amyg[which(amyg$phase %in% c('baseline')),]
amyg.pre.mod <- lmer(ers ~ condition*roi*group + (1|subject), data=amyg.pre, REML=FALSE)
summary(amyg.pre.mod)
anova(amyg.pre.mod)
#planned comparisons of CS+ vs. CS-
amyg.pre.csdif <- emmeans(amyg.pre.mod, revpairwise ~ condition|roi*group)
confint(amyg.pre.csdif) #point-estimate and confidence intervals
summary(amyg.pre.csdif) #t values
p.adjust(summary(amyg.pre.csdif$contrasts)$p.value, method="fdr") #fdr corrected p-values

#shows CeM > BLA for ptsd for pre-conditioning but not very interpretable
emmeans(amyg.pre.mod, pairwise ~ roi|group) 

#emotional reinstatement
amyg.emo <- amyg[which(amyg$phase %in% c('acquisition','extinction')),]
amyg.emo.mod <- lmer(ers ~ phase*condition*roi*group + (1|subject), data=amyg.emo, REML=FALSE)
summary(amyg.emo.mod)
anova(amyg.emo.mod)
#planned comparisons of CS+ vs. CS-
amyg.emo.csdif <- emmeans(amyg.emo.mod, revpairwise ~ condition|phase*roi*group, adjust=FALSE, lmer.df="satterthwaite")
confint(amyg.emo.csdif) #point-estimate and confidence intervals
summary(amyg.emo.csdif) #t values
p.adjust(summary(amyg.emo.csdif$contrasts)$p.value, method="fdr") #fdr corrected p-values

#no point doing post-hoc within the amygdala

#hippocampus------------------------------------------------------------------------------------
hpc <- subcort[which(subcort$roi %in% c('hc_head_ers','hc_body_ers','hc_tail_ers')),]
#i do this so the phases are in order alphabetically (baseline, conditioning, extinction)
hpc <- hpc %>% mutate(phase = recode(phase, "acquisition" = "conditioning"))

hpc.mod <- mixed(ers ~ condition*phase*roi*group + (1|subject), data=hpc, REML=FALSE, method="LRT")#, test_intercept=TRUE))
summary(hpc.mod)
anova(hpc.mod)

hpc.phasedif <- emmeans(hpc.mod, revpairwise ~ phase|group*roi, adjust="None", lmer.df="asymptotic")
summary(hpc.phasedif)
p.adjust(summary(hpc.phasedif$contrasts)$p.value, method="fdr") <.05#fdr corrected p-values



#context evidence prediction------------------------------------------------------------------
df <- read.csv('all_data_lmm.csv')
ctx.emo <- df[which(df$phase %in% c('acquisition','extinction')),]
ctx.emo.mod <- lmer(ev ~ phase*condition*group + (1|subject), data=ctx.emo, REML=FALSE)
summary(ctx.emo.mod)
anova(ctx.emo.mod)


#doesn't survive
ctx.emo.phasedif <- emmeans(ctx.mod, pairwise ~ phase|group, adjust=FALSE)
p.adjust(summary(ctx.emo.phasedif$contrasts)$p.value, method="fdr") #fdr corrected p-values

#using to predict the pfc difference split metric
ctx.pfc.mod <- lmer(pfc_diff_ers ~ ev*phase*condition*group + (1|subject), data = ctx.emo, REML=FALSE)
summary(ctx.pfc.mod)
anova(ctx.pfc.mod)

ctx.cstrend <- emtrends(ctx.pfc.mod, ~ condition, var="ev")
summary(ctx.cstrend)
test(ctx.cstrend)

ctx.evcon <- emmip(ctx.pfc.mod,condition~ev, at=list(ev=seq(0,1,by=0.001)), CIs=TRUE, plotit=FALSE)
ggplot(data=ctx.evcon, aes(x=ev,y=yvar, color=condition)) + geom_line() + geom_ribbon(aes(ymax=UCL, ymin=LCL, fill=condition), alpha=0.4)

#workspace------------------------------------------------------------------------------------


mod <- lmer(pfc_diff_ers~ev*phase*condition*group + (1|subject),data=ctx.emo,REML=FALSE)
anova(mod)
summary(mod)
modsum <- summary(mod)
#anova(basic,mod)

a <- probe_interaction(mod,"ev",modx="condition",cond.int=TRUE,interval=TRUE,centered="none",colors=con_pal)
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
