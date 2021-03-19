
#analyze preconditioning (baseline) first
pfc.pre <- pfc[which(pfc$phase %in% c('baseline')),]
pfc.pre.mod <- lmer(ers ~ condition*roi*group + (1|subject), data=pfc.pre, REML=FALSE)
summary(pfc.pre.mod)
anova(pfc.pre.mod)
#planned comparisons of CS+ vs. CS-
pfc.pre.csdif <- emmeans(pfc.pre.mod, revpairwise ~ condition|roi*group)
confint(pfc.pre.csdif) #point-estimate and confidence intervals
summary(pfc.pre.csdif) #t values
p.adjust(summary(pfc.pre.csdif$contrasts)$p.value, method="fdr") #fdr corrected p-values


#emotional reinstatement
pfc.emo <- pfc[which(pfc$phase %in% c('acquisition','extinction')),]
pfc.emo.mod <- lmer(ers ~ phase*condition*roi*group + (1|subject), data=pfc.emo, REML=FALSE)
summary(pfc.emo.mod)
anova(pfc.emo.mod)
#planned comparisons of CS+ vs. CS-
pfc.emo.csdif <- emmeans(pfc.emo.mod, revpairwise ~ condition|phase*roi*group, adjust=FALSE, lmer.df="satterthwaite")
confint(pfc.emo.csdif) #point-estimate and confidence intervals
summary(pfc.emo.csdif) #t values
p.adjust(summary(pfc.emo.csdif$contrasts)$p.value, method="fdr") #fdr corrected p-values

#post-hoc between roi and between group comparisons

#these are manual codes that set up the double subtraction we want
#position in the vectors corresponds to row in the pfc.emo.csmean object
#analyze preconditioning (baseline) first
hpc.pre <- hpc[which(hpc$phase %in% c('baseline')),]
hpc.pre.mod <- lmer(ers ~ condition*roi*group + (1|subject), data=hpc.pre, REML=FALSE)
summary(hpc.pre.mod)
anova(hpc.pre.mod)
#planned comparisons of CS+ vs. CS-
hpc.pre.csdif <- emmeans(hpc.pre.mod, revpairwise ~ condition|roi*group)
confint(hpc.pre.csdif) #point-estimate and confidence intervals
summary(hpc.pre.csdif) #t values
p.adjust(summary(hpc.pre.csdif$contrasts)$p.value, method="fdr") #fdr corrected p-values

#emotional reinstatement
hpc.emo <- hpc[which(hpc$phase %in% c('acquisition','extinction')),]
hpc.emo.mod <- lmer(ers ~ phase*condition*roi*group + (1|subject), data=hpc.emo, REML=FALSE)
summary(hpc.emo.mod)
anova(hpc.emo.mod)
#planned comparisons of CS+ vs. CS-
hpc.emo.csdif <- emmeans(hpc.emo.mod, revpairwise ~ condition|phase*roi*group, adjust=FALSE, lmer.df="satterthwaite")
confint(hpc.emo.csdif) #point-estimate and confidence intervals
summary(hpc.emo.csdif) #t values
p.adjust(summary(hpc.emo.csdif$contrasts)$p.value, method="fdr") #fdr corrected p-values

#phase by ROI comparisons
hpc.emo.phasedif <- emmeans(hpc.emo.mod, pairwise ~ phase|roi, ajust=FALSE, lmer.df="satterthwaite")
confint(hpc.emo.phasedif)
summary(hpc.emo.phasedif)
p.adjust(summary(hpc.emo.phasedif$contrasts)$p.value, method="fdr") #fdr corrected p-values

#phase difference and group difference
hpc.emo.groupdif <- emmeans(hpc.emo.mod, pairwise ~ group|phase*roi, ajust=FALSE, lmer.df="satterthwaite")
summary(hpc.emo.groupdif)
p.adjust(summary(hpc.emo.groupdif$contrasts)$p.value, method="fdr") #fdr corrected p-values



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
