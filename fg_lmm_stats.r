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
require(ez)
afex::set_sum_contrasts()

phases <- c('baseline','acquisition','extinction')
emo_phases <- c('acquisition','extinction')
cons <- c('CS+','CS-')
groups <- c('healthy','ptsd')
group_pal <- colors()[c(565,555)]
phase_pal <- colors()[c(454,614)]
con_pal <- colors()[c(17,92)]
emm_options(lmer.df="asymptotic")

#setwd('/Users/ach3377/Documents/gPPI')
setwd('C:\\Users\\ACH\\Documents\\gPPI')

pfc <- read.csv('pfc_ers_cleaned_lmm.csv')
subcort <- read.csv('subcort_ers_cleaned_lmm.csv')
cb <- read.csv('cb_response_rsa.csv')
cb <- cb %>% mutate(phase = recode(phase, "acquisition" = "conditioning"))
ins <- cb[which(cb$roi == 'ant_ins'),]
precun <- cb[which(cb$roi == 'precun'),]
#i do this so the phases are in order alphabetically (baseline, conditioning, extinction)
pfc <- pfc %>% mutate(phase = recode(phase, "acquisition" = "conditioning"))
subcort <- subcort %>% mutate(phase = recode(phase, "acquisition" = "conditioning"))


#PFC encoding-retrieval similarity----------------------------------------------------------------------
pfc.mod <- mixed(ers ~ condition*phase*roi*group + (1|subject), data=pfc, REML=FALSE, method="LRT")#, test_intercept=TRUE)
pfc.mod.int <- emmeans(pfc.mod, ~1)
summary(pfc.mod.int)
test(pfc.mod.int)$p.value
anova(pfc.mod)
  
  #planned comparisons of CS+ vs. CS-
pfc.csdif <- emmeans(pfc.mod, revpairwise ~ condition|phase*roi*group, adjust="None")
pfc.csdif.ci <- confint(pfc.csdif$contrasts)
summary(pfc.csdif$contrasts)
p.adjust(summary(pfc.csdif$contrasts)$p.value, method="fdr")

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

pfc.con.list <- list(healthy_acq_dACC - healthy_acq_vmPFC,
                     healthy_ext_vmPFC - healthy_ext_dACC,
                     ptsd_acq_dACC - ptsd_acq_vmPFC,
                     ptsd_ext_vmPFC - ptsd_ext_dACC)

pfc.cons <- contrast(pfc.csmean, method=pfc.con.list, adjust="None")
confint(pfc.cons)
summary(pfc.cons)
p.adjust(summary(pfc.cons)$p.value, method="fdr")

#these are the within ROI comparisons for between phase
pfc.con.list2 <- list(healthy_acq_dACC - healthy_ext_dACC,
                      healthy_ext_vmPFC - healthy_acq_vmPFC,
                      ptsd_acq_dACC - ptsd_ext_dACC,
                      ptsd_ext_vmPFC - ptsd_acq_vmPFC)
pfc.cons2 <- contrast(pfc.csmean, method=pfc.con.list2, adjust="None")
confint(pfc.cons2)
summary(pfc.cons2)

#just really focused interaction
hdf <- pfc[which(pfc$group == 'healthy' & pfc$phase %in% c('conditioning','extinction')),]
hdf.mod <- mixed(ers ~ condition*phase*roi + (1|subject), data=hdf, REML=FALSE, method="LRT")
anova(hdf.mod)

pdf <- pfc[which(pfc$group == 'ptsd' & pfc$phase %in% c('conditioning','extinction')),]
pdf.mod <- mixed(ers ~ condition*phase*roi + (1|subject), data=pdf, REML=FALSE, method="LRT")
anova(pdf.mod)

#######################################################################
#anterior insula
ins.mod <- mixed(ers ~ condition*phase*group + (1|subject), data=ins, REML=FALSE, method="LRT")
anova(ins.mod)
ins.csdif <- emmeans(ins.mod, revpairwise ~ condition|phase*group, adjust="None")
summary(ins.csdif)
p.adjust(summary(ins.csdif$contrasts)$p.value, method="fdr")


#precuneus
precun.mod <- mixed(ers ~ condition*phase*group + (1|subject), data=precun, REML=FALSE, method="LRT")
anova(precun.mod)
precun.csdif <- emmeans(precun.mod, revpairwise ~ condition|phase*group, adjust="None")
summary(precun.csdif)
p.adjust(summary(precun.csdif$contrasts)$p.value, method="fdr")



# pfc US reinforcement ----------------------------------------------------
acq <- pfc[which(pfc$phase %in% c('conditioning') & pfc$condition == 'CS+'),]

us.mod <- mixed(ers ~ shock*roi*group + (1|subject),data=acq, REML=FALSE, method="LRT")
anova(us.mod)

us.means <- emmeans(us.mod, ~ shock|roi*group)
emmip(us.mod, shock~roi*group, CIs=TRUE)

sub.acq <- subcort[which(subcort$phase == 'conditioning' & subcort$condition == 'CS+'),]
amyg.acq <- sub.acq[which(sub.acq$roi %in% c('amyg_bla','amyg_cem')),]
hpc.acq <- sub.acq[which(sub.acq$roi %in% c('hc_head','hc_body','hc_tail')),]

amyg.us.mod <- mixed(ers ~ shock*roi*group + (1|subject), data=amyg.acq, REML=FALSE, method="LRT")
anova(amyg.us.mod)
emmeans(amyg.us.mod, pairwise ~ shock|roi*group, adjust="None")
emmip(amyg.us.mod, shock~roi*group, CIs=TRUE)

hpc.us.mod <- mixed(ers ~ shock*roi*group + (1|subject), data=hpc.acq, REML = FALSE, method="LRT")
anova(hpc.us.mod)
emmip(hpc.us.mod, shock~roi*group, CIs=TRUE)

#amygdala------------------------------------------------------------------------------------
amyg <- subcort[which(subcort$roi %in% c('amyg_cem','amyg_bla')),]

amyg.mod <- mixed(ers ~ condition*phase*roi*group + (1|subject), data=amyg, REML=FALSE, method="LRT")
anova(amyg.mod)

amyg.mod.int <- emmeans(amyg.mod, ~1)
test(amyg.mod.int)
summary(amyg.mod.int)

#doing the cs comps for continuity
amyg.cond <- emmeans(amyg.mod, revpairwise ~ condition|phase*roi*group, adjust="None")
summary(amyg.cond)
amyg.p <- p.adjust(summary(amyg.cond$contrasts)$p.value, method="fdr")



#hippocampus------------------------------------------------------------------------------------
hpc <- subcort[which(subcort$roi %in% c('hc_head','hc_body','hc_tail')),]

hpc.mod <- mixed(ers ~ condition*phase*roi*group + (1|subject), data=hpc, REML=FALSE, method="LRT")#, test_intercept=TRUE))
summary(hpc.mod)
anova(hpc.mod)

hpc.mod.int <- emmeans(hpc.mod, ~1)
summary(hpc.mod.int)
test(hpc.mod.int)$p.value

hpc.phasedif <- emmeans(hpc.mod, revpairwise ~ phase|group*roi, adjust="None")
confint(hpc.phasedif$contrasts)[p.adjust(summary(hpc.phasedif$contrasts)$p.value, method="fdr") < .05,]
summary(hpc.phasedif$contrasts)[p.adjust(summary(hpc.phasedif$contrasts)$p.value, method="fdr") < .05,]
p.adjust(summary(hpc.phasedif$contrasts)$p.value, method="fdr")[p.adjust(summary(hpc.phasedif$contrasts)$p.value, method="fdr") < .05] #fdr corrected p-values

#testing double dissociation in hippocampus
hpc.phasemean <- emmeans(hpc.mod, ~phase*roi*group)

healthy_acq_head <- c(0,0,0, 0,1,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0) 
healthy_ext_head <- c(0,0,0, 0,0,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0) 
healthy_acq_tail <- c(0,0,0, 0,0,0, 0,1,0, 0,0,0, 0,0,0, 0,0,0) 
healthy_ext_tail <- c(0,0,0, 0,0,0, 0,0,1, 0,0,0, 0,0,0, 0,0,0) 
  
ptsd_acq_head    <- c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,1,0, 0,0,0) 
ptsd_ext_head    <- c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,1, 0,0,0) 
ptsd_acq_tail    <- c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,1,0) 
ptsd_ext_tail    <- c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,1)

hpc.con.list <- list(healthy_acq_tail - healthy_acq_head,
                     healthy_ext_head - healthy_ext_tail,
                     ptsd_acq_tail - ptsd_acq_head,
                     ptsd_ext_head - ptsd_ext_tail)

hpc.cons <- contrast(hpc.phasemean, method=hpc.con.list, adjust="None")
confint(hpc.cons)
summary(hpc.cons)
p.adjust(summary(hpc.cons)$p.value, method="fdr")

#doing all 3 bc maybe it will work
hpc.roidif <- emmeans(hpc.mod, pairwise ~ roi|phase*group, adjust="None")
summary(hpc.roidif$contrasts)[p.adjust(summary(hpc.roidif$contrasts)$p.value, method='fdr') < .05,]

#just showing that the condition comp is null in the hippocampus
hpc.cond <- emmeans(hpc.mod, revpairwise ~ condition|phase*roi*group, adjust="None")
summary(hpc.cond)
p.adjust(summary(hpc.cond$contrasts)$p.value, method="fdr")

#just the interaction I care about for realz
hdf <- hpc[which(hpc$group == 'healthy' & hpc$phase %in% c('conditioning','extinction') & hpc$roi %in% c('hc_tail','hc_head')),]
hdf.mod <- mixed(ers ~ condition*phase*roi + (1|subject), data=hdf, REML=FALSE, method="LRT")
anova(hdf.mod)

pdf <- hpc[which(hpc$group == 'ptsd' & hpc$phase %in% c('conditioning','extinction') & hpc$roi %in% c('hc_tail','hc_head')),]
pdf.mod <- mixed(ers ~ condition*phase*roi + (1|subject), data=pdf, REML=FALSE, method="LRT")
anova(pdf.mod)


#subcortical univariate predictions------------------------------------------------------------------
df <- read.csv('all_data_lmm.csv')
emo <- df[which(df$phase %in% c('acquisition','extinction')),] #emo = emotional contexts

tail.uni.mod <- mixed(pfc_diff_ers ~ hc_tail_ret_uni*condition*phase*group + (1|subject), data=emo, REML=FALSE, method="LRT")
body.uni.mod <- mixed(pfc_diff_ers ~ hc_body_ret_uni*condition*phase*group + (1|subject), data=emo, REML=FALSE, method="LRT")
head.uni.mod <- mixed(pfc_diff_ers ~ hc_head_ret_uni*condition*phase*group + (1|subject), data=emo, REML=FALSE, method="LRT")
bla.uni.mod <- mixed(pfc_diff_ers ~ amyg_bla_ret_uni*condition*phase*group + (1|subject), data=emo, REML=FALSE, method="LRT")
cem.uni.mod <- mixed(pfc_diff_ers ~ amyg_cem_ret_uni*condition*phase*group + (1|subject), data=emo, REML=FALSE, method="LRT")

pred.rows <- c(1,5,6,8,11,12,13,15) #we just care about main effect and interactions with the continuous predictor
anova(tail.uni.mod)[pred.rows,]
tail.uni.slope <- emtrends(tail.uni.mod, ~1, var="hc_tail_ret_uni")
tail.uni.csdif <- emtrends(tail.uni.mod, revpairwise~condition, var="hc_tail_ret_uni")
summary(tail.uni.csdif$contrasts)$p.value
confint(tail.uni.csdif)
tail.uni.out <- emtrends(tail.uni.mod, ~condition|phase, var="hc_tail_ret_uni")
write.csv(data.frame(tail.uni.out),'hc_tail_slopes_con_phase.csv')

anova(body.uni.mod)[pred.rows,]
body.uni.slope <- emtrends(body.uni.mod, ~1, var="hc_body_ret_uni")
body.uni.dif <- emtrends(body.uni.mod, revpairwise ~ condition|phase, var="hc_body_ret_uni")
p.adjust(summary(body.uni.dif$contrasts)$p.value, method="fdr")
confint(body.uni.dif$contrasts)
body.uni.out <- emtrends(body.uni.mod, ~condition|phase, var="hc_body_ret_uni")
write.csv(data.frame(body.uni.out),'hc_body_slopes_con_phase.csv')


anova(head.uni.mod)[pred.rows,]
head.uni.slope <- emtrends(head.uni.mod, ~1, var="hc_head_ret_uni")

anova(bla.uni.mod)[pred.rows,]
bla.uni.slope <- emtrends(bla.uni.mod, ~1, var="amyg_bla_ret_uni")

anova(cem.uni.mod)[pred.rows,]
cem.uni.slope <- emtrends(cem.uni.mod, ~1, var="amyg_cem_ret_uni")

uni.slopes.out <- list(
  as.numeric(as.data.frame(tail.uni.slope)[1,c(2,5,6)]),
  as.numeric(as.data.frame(body.uni.slope)[1,c(2,5,6)]),
  as.numeric(as.data.frame(head.uni.slope)[1,c(2,5,6)]),
  as.numeric(as.data.frame(bla.uni.slope)[1,c(2,5,6)]),
  as.numeric(as.data.frame(cem.uni.slope)[1,c(2,5,6)]))
names(uni.slopes.out) <- c("tail","body","head","bla","cem")
write.csv(as.data.frame(uni.slopes.out),"uni_slopes_pfc_pred.csv")

#subcortical ERS predictions------------------------------------------------------------------
df <- read.csv('all_data_lmm.csv')
emo <- df[which(df$phase %in% c('acquisition','extinction')),] #emo = emotional contexts

tail.ers.mod <- mixed(pfc_diff_ers ~ hc_tail_ers*condition*phase*group + (1|subject), data=emo, REML=FALSE, method="LRT")
body.ers.mod <- mixed(pfc_diff_ers ~ hc_body_ers*condition*phase*group + (1|subject), data=emo, REML=FALSE, method="LRT")
head.ers.mod <- mixed(pfc_diff_ers ~ hc_head_ers*condition*phase*group + (1|subject), data=emo, REML=FALSE, method="LRT")
bla.ers.mod <- mixed(pfc_diff_ers ~ amyg_bla_ers*condition*phase*group + (1|subject), data=emo, REML=FALSE, method="LRT")
cem.ers.mod <- mixed(pfc_diff_ers ~ amyg_cem_ers*condition*phase*group + (1|subject), data=emo, REML=FALSE, method="LRT")

pred.rows <- c(1,5,6,8,11,12,13,15) #we just care about main effect and interactions with the continuous predictor

anova(tail.ers.mod)[pred.rows,]
tail.ers.slope <- emtrends(tail.ers.mod, ~1, var="hc_tail_ers")
tail.ers.dif <- emtrends(tail.ers.mod, pairwise~group|phase, var="hc_tail_ers", adjust="None")
p.adjust(test(tail.ers.dif)$p.value)
tail.ers.vis <- emmip(tail.ers.mod, phase|group~hc_tail_ers, at=list(hc_tail_ers=seq(-.2,.3,by=0.01)), CIs=TRUE, plotit=F)
tail.vis.out <- ggplot(data=tail.ers.vis, aes(x=hc_tail_ers,y=yvar, color=phase, linetype=group)) + geom_line() + geom_ribbon(aes(ymax=UCL, ymin=LCL, fill=phase), alpha=0.4)
write.csv(tail.vis.out$data,'hc_tail_pfc_diff_lines.csv')

tail.ers.phasedif <- emtrends(tail.ers.mod, revpairwise~phase|condition*group, var="hc_tail_ers")


anova(body.ers.mod)[pred.rows,]
body.ers.slope <- emtrends(body.ers.mod, ~1, var="hc_body_ers")

anova(head.ers.mod)[pred.rows,]
head.ers.slope <- emtrends(head.ers.mod, ~1, var="hc_head_ers")
head.ers.dif <- emtrends(head.ers.mod, pairwise~phase|group, var="hc_head_ers")
confint(head.ers.dif)
test(head.ers.dif)
p.adjust(test(head.ers.dif)$p.value)
head.ers.vis <- emmip(head.ers.mod, phase|group~hc_head_ers, at=list(hc_head_ers=seq(-.3,.5,by=0.01)), CIs=TRUE, plotit=F)
head.vis.out <- ggplot(data=head.ers.vis, aes(x=hc_head_ers,y=yvar, color=phase, linetype=group)) + geom_line() + geom_ribbon(aes(ymax=UCL, ymin=LCL, fill=phase), alpha=0.4)
write.csv(head.vis.out$data,'hc_head_pfc_diff_lines.csv')

anova(bla.ers.mod)[pred.rows,]
bla.ers.slope <- emtrends(bla.ers.mod, ~1, var="amyg_bla_ers")
bla.ers.dif <- emtrends(bla.ers.mod, pairwise~condition|group, var="amyg_bla_ers")
test(bla.ers.dif)


anova(cem.ers.mod)[pred.rows,]
cem.ers.slope <- emtrends(cem.ers.mod, ~1, var="amyg_cem_ers")
cem.ers.dif <- emtrends(cem.ers.mod, ~condition|phase*group, var="amyg_cem_ers")
test(cem.ers.dif)

#hello.mod <- mixed(pfc_diff_ers ~ hc_head_ret_uni*hc_head_ers*condition*phase*group + (1|subject), data=emo, REML=FALSE, method="LRT")
ers.slopes.out <- list(
  as.numeric(as.data.frame(tail.ers.slope)[1,c(2,5,6)]),
  as.numeric(as.data.frame(body.ers.slope)[1,c(2,5,6)]),
  as.numeric(as.data.frame(head.ers.slope)[1,c(2,5,6)]),
  as.numeric(as.data.frame(bla.ers.slope)[1,c(2,5,6)]),
  as.numeric(as.data.frame(cem.ers.slope)[1,c(2,5,6)]))
names(ers.slopes.out) <- c("tail","body","head","bla","cem")
write.csv(as.data.frame(ers.slopes.out),"ers_slopes_pfc_pred.csv")

#anterior double mod-------------------------------------------------------------
head.combo.mod <- mixed(pfc_diff_ers ~ hc_head_ers*hc_head_ret_uni*condition*phase*group + (1|subject), data=emo, REML=FALSE, method="LRT")
anova(head.combo.mod)
head.combo.ers <- emtrends(head.combo.mod, ~1, var="hc_head_ers")
head.combo.uni <- emtrends(head.combo.mod, ~1, var="hc_head_ret_uni")


#context evidence prediction------------------------------------------------------------------
df <- read.csv('all_data_lmm.csv')
emo <- df[which(df$phase %in% c('acquisition','extinction')),] #emo = emotional contexts

#maybe cut
ctx.solo <- mixed(ev ~ condition*phase*group + (1|subject), data=emo, REML=FALSE, method ="LRT")
anova(ctx.solo)

ctx.mod <- mixed(pfc_diff_ers ~ ev*condition*phase*group + (1|subject), data=emo, REML=FALSE, method="LRT")
ctx.slope <- emtrends(ctx.mod, ~1, var="ev")
anova(ctx.mod)
ctx.cs <- emtrends(ctx.mod, ~condition, var="ev")
test(ctx.cs)
p.adjust(test(ctx.cs)$p.value)

ctx.slopes.out <-  as.data.frame(ctx.cs)[c(1,2),c(1,2,5,6)]
write.csv(ctx.slopes.out, 'ev_slopes_pfc_pred.csv')


#ctx.cs.vis <- emmip(effectsize::standardize(ctx.mod,select="ev",include_response=FALSE), condition~ev, at=list(ev=seq(-1,1,by=0.01)), CIs=TRUE, plotit=F)
#ctx.cs.out <- ggplot(data=ctx.cs.vis, aes(x=ev,y=yvar, color=condition)) + geom_line() + geom_ribbon(aes(ymax=UCL, ymin=LCL, fill=condition), alpha=0.4)
#write.csv(ctx.cs.out$data,'ev_pfc_diff_lines.csv')



#analysis of recognition memory------------------------------------------------------------------
pfc <- read.csv('pfc_ers_cleaned_lmm.csv')
#i do this so the phases are in order alphabetically (baseline, conditioning, extinction)
pfc <- pfc %>% mutate(phase = recode(phase, "acquisition" = "conditioning"))

pfc.mod <- mixed(ers ~ mem_acc*condition*phase*roi*group + (1|subject), data=pfc, REML=FALSE, method="LRT")#, test_intercept=TRUE)
pfc.mod.int <- emmeans(pfc.mod, ~1)
summary(pfc.mod.int)
test(pfc.mod.int)$p.value
anova(pfc.mod)


beh <- pfc %>% mutate(mem_acc = recode(mem_acc, "H" = 1, "M" = 0))
beh.mod <- mixed(mem_acc ~ condition*phase*group + (1|subject), data = beh, method="LRT", family="binomial")
anova(beh.mod)
high <- glmer(mem_acc ~ condition*phase*group + (1|subject), family="binomial", data=beh)

#amygdala memory------------------------------------------------------------------------------------
amyg <- subcort[which(subcort$roi %in% c('amyg_cem','amyg_bla')),]

amyg.mod <- mixed(ers ~ mem_acc*condition*phase*roi*group + (1|subject), data=amyg, REML=FALSE, method="LRT")
anova(amyg.mod)

amyg.mod.int <- emmeans(amyg.mod, ~1)
test(amyg.mod.int)
summary(amyg.mod.int)

dif <- emmeans(amyg.mod, pairwise ~ mem_acc|condition*phase*roi)
dif.ci <- confint(dif$contrasts)
summary(dif$contrasts)
p.adjust(summary(pfc.csdif$contrasts)$p.value, method="fdr")


# hippocampus memory ------------------------------------------------------
hpc <- subcort[which(subcort$roi %in% c('hc_head','hc_body','hc_tail')),]
#i do this so the phases are in order alphabetically (baseline, conditioning, extinction)
hpc <- hpc %>% mutate(phase = recode(phase, "acquisition" = "conditioning"))

hpc.mod <- mixed(ers ~ mem_acc*condition*phase*roi*group + (1|subject), data=hpc, REML=FALSE, method="LRT")#, test_intercept=TRUE))
summary(hpc.mod)
anova(hpc.mod)

dif <- emmeans(hpc.mod, revpairwise ~ mem_acc|phase*group*condition)
summary(dif$contrasts)

# pfc_diff conn --------------------------------------------------------------
df <- read.csv('conn_pfc_diff_lmm.csv')
df <- df[which(df$phase %in% c('acquisition','extinction')),]

hc_head.mod <- mixed(conn ~ condition*phase*group + (1|subject), data=df[which(df$seed %in% c('hc_head')),], REML=FALSE, method="LRT")
anova(hc_head.mod)
test(emmeans(hc_head.mod, ~ condition*phase*group,adjust="None"))

hc_tail.mod <- mixed(conn ~ condition*phase*group + (1|subject), data=df[which(df$seed %in% c('hc_tail')),], REML=FALSE, method="LRT")
anova(hc_tail.mod)
test(emmeans(hc_tail.mod, ~ condition*phase*group,adjust="None"))

amyg_bla.mod <- mixed(conn ~ condition*phase*group + (1|subject), data=df[which(df$seed %in% c('amyg_bla')),], REML=FALSE, method="LRT")
anova(amyg_bla.mod)
emmeans(amyg_bla.mod, ~ condition*phase*group,adjust="None")

amyg_cem.mod <- mixed(conn ~ condition*phase*group + (1|subject), data=df[which(df$seed %in% c('amyg_cem')),], REML=FALSE, method="LRT")
anova(amyg_cem.mod)
emmeans(amyg_cem.mod, revpairwise ~ phase|group,adjust="None")


#conn stats--------------------------------------------------------------
df <- read.csv('conn_stats_lmm.csv')
df <- df[which(df$phase %in% c('acquisition','extinction')),]
V <- df[which(df$target == 'vmPFC'),]
D <- df[which(df$target == 'dACC'),]
BLA <- df[which(df$target == 'amyg_bla'),]
CEM <- df[which(df$target == 'amyg_cem'),]
aHPC <- df[which(df$target == 'hc_head'),]

head_to_pfc <- df[which(df$seed == 'hc_head' & df$target %in% c('vmPFC','dACC')),]
bla_to_pfc <- df[which(df$seed == 'amyg_bla' & df$target %in% c('vmPFC','dACC')),]
cem_to_pfc <- df[which(df$seed == 'amyg_cem' & df$target %in% c('vmPFC','dACC')),]
tail_to_pfc <- df[which(df$seed == 'hc_tail' & df$target %in% c('vmPFC','dACC')),]
head_to_ext <- df[which(df$seed == 'hc_head' & df$target %in% c('vmPFC','dACC','amyg_bla')),]



v_to_sub <- df[which(df$seed == 'vmPFC' & df$target %in% c('amyg_cem','amyg_bla','hc_head','hc_tail')),]
d_to_sub <- df[which(df$seed == 'dACC' & df$target %in% c('amyg_cem','amyg_bla','hc_head','hc_tail')),]

hc_head.D <- mixed(conn ~ condition*phase*group + (1|subject), data=D[which(D$seed == 'hc_head'),], REML=FALSE, method="LRT")
anova(hc_head.D)
hc_head.D.em <- emmeans(hc_head.D, revpairwise ~ phase|group,adjust="None")
summary(hc_head.D.em)
p.adjust(summary(hc_head.D.em$contrasts)$p.value, method="fdr")


hc_head.V <- mixed(conn ~ condition*phase*group + (1|subject), data=V[which(V$seed == 'hc_head'),], REML=FALSE, method="LRT")
anova(hc_head.V)
hc_head.V.em <- emmeans(hc_head.V, revpairwise ~ phase|group,adjust="None")
summary(hc_head.V.em)
p.adjust(summary(hc_head.V.em$contrasts)$p.value, method="fdr")


hc_tail.aHPC <- mixed(conn ~ condition*phase*group + (1|subject), data=aHPC[which(aHPC$seed == 'hc_tail'),], REML=FALSE, method="LRT")
anova(hc_tail.aHPC)
hc_tail.aHPC.em <- emmeans(hc_tail.aHPC, revpairwise ~ phase|group,adjust="None")
summary(hc_tail.aHPC.em)
p.adjust(summary(hc_head.V.em$contrasts)$p.value, method="fdr")



hc_head.BLA <- mixed(conn ~ condition*phase*group + (1|subject), data=BLA[which(BLA$seed == 'hc_head'),], REML=FALSE, method="LRT")
anova(hc_head.BLA)
test(emmeans(hc_head.BLA, revpairwise ~ condition|phase*group,adjust="None"))


vmPFC.BLA <- mixed(conn ~ condition*phase*group + (1|subject), data=BLA[which(BLA$seed == 'vmPFC'),], REML=FALSE, method="LRT")
anova(vmPFC.BLA)
vmPFC.BLA.em <- emmeans(vmPFC.BLA, revpairwise ~ condition|phase*group, adjust="None")
summary(vmPFC.BLA.em)
p.adjust(summary(vmPFC.BLA.em$contrasts)$p.value, method="fdr")

dACC.BLA <- mixed(conn ~ condition*phase*group + (1|subject), data=BLA[which(BLA$seed == 'dACC'),], REML=FALSE, method="LRT")
anova(dACC.BLA)
emmeans(dACC.BLA, revpairwise ~ condition*phase*group, adjust="None")

BLA.CEM <- mixed(conn ~ condition*phase*group + (1|subject), data=CEM[which(CEM$seed == 'amyg_bla'),], REML=FALSE, method="LRT")
anova(BLA.CEM)
BLA.CEM.em <- emmeans(BLA.CEM, revpairwise ~ condition|phase*group, adjust="None")
summary(BLA.CEM.em)
p.adjust(summary(BLA.CEM.em$contrasts)$p.value, method="fdr")


head.pfc.mod <- mixed(conn ~ condition*phase*target*group + (1|subject), data=head_to_pfc, REML=FALSE, method="LRT")
anova(head.pfc.mod)
head.pfc.em <- emmeans(head.pfc.mod, revpairwise ~ target|phase*group, adjust="None")
summary(head.pfc.em)
p.adjust(summary(head.pfc.em$contrasts)$p.value, method="fdr")

head.ext.mod <- mixed(conn ~ condition*phase*target*group + (1|subject), data=head_to_ext, REML=FALSE, method="LRT")
anova(head.ext.mod)
emmeans(head.ext.mod, revpairwise ~ target|group, adjust="None")

bla.pfc.mod <- mixed(conn ~ condition*phase*target*group + (1|subject), data=bla_to_pfc, REML=FALSE, method="LRT")
anova(bla.pfc.mod)
emmeans(bla.pfc.mod, revpairwise ~ phase, adjust="None")

cem.pfc.mod <- mixed(conn ~ condition*phase*target*group + (1|subject), data=cem_to_pfc, REML=FALSE, method="LRT")
anova(cem.pfc.mod)
emmeans(cem.pfc.mod, revpairwise ~ target|group, adjust="None")

tail.pfc.mod <- mixed(conn ~ condition*phase*target*group + (1|subject), data=tail_to_pfc, REML=FALSE, method="LRT")
anova(tail.pfc.mod)
tail.pfc.em <- emmeans(tail.pfc.mod, revpairwise ~ target|group, adjust="None")
summary(tail.pfc.em)
p.adjust(summary(tail.pfc.em$contrasts)$p.value, method="fdr")


a_to_pfc <- df[which(df$seed == 'animal' & df$target %in% c('dACC') & df$csp_cat == 'animal' & df$condition == 'CS+'),]
t_to_pfc <- df[which(df$seed == 'tool' & df$target %in% c('dACC') & df$csp_cat == 'tool' & df$condition == 'CS+'),]
cat_to_pfc <- bind_rows(a_to_pfc,t_to_pfc)

a.pfc.mod <- mixed(conn ~ phase*group + (1|subject), data=a_to_pfc, REML=FALSE, method="LRT")
anova(a.pfc.mod)
a.pfc.em <- emmeans(a.pfc.mod, revpairwise ~ phase|group, adjust="None")
summary(a.pfc.em)
p.adjust(summary(a.pfc.em$contrasts)$p.value, method="fdr")     

t.pfc.mod <- mixed(conn ~ phase*group + (1|subject), data=t_to_pfc, REML=FALSE, method="LRT")
anova(t.pfc.mod)
t.pfc.em <- emmeans(t.pfc.mod, revpairwise ~ phase|group, adjust="None")
summary(t.pfc.em)
p.adjust(summary(t.pfc.em$contrasts)$p.value, method="fdr")

cat.pfc.mod <- mixed(conn ~ phase*group + (1|subject), data=cat_to_pfc, REML=FALSE, method="LRT")
anova(cat.pfc.mod)
cat.pfc.em <- emmeans(cat.pfc.mod, revpairwise ~ phase|group, adjust="None")
summary(cat.pfc.em)
p.adjust(summary(cat.pfc.em$contrasts)$p.value, method="fdr")

v.sub.mod <- mixed(conn ~ condition*phase*target*group + (1|subject), data=v_to_sub, REML=FALSE, method="LRT")
anova(v.sub.mod)
emmeans(v.sub.mod, revpairwise ~ target|phase*group, adjust="None")

d.sub.mod <- mixed(conn ~ condition*phase*target*group + (1|subject), data=d_to_sub, REML=FALSE, method="LRT")
anova(d.sub.mod)
emmeans(d.sub.mod, revpairwise ~ condition|group, adjust="None")

# memory behavioral data --------------------------------------------------
df <- read.csv('mem_hit_rate_sub_means.csv')
mem.res <- ezANOVA(df,dv=.(mem_acc),within=.(phase,condition),between=.(group),wid=.(subject),type=3)
