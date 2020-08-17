require(lme4)
require(ez)
require(afex)
require(lmerTest) #Package must be installed first
require(emmeans)
`%notin%` <- Negate(`%in%`)

setwd("C:/Users/ACH/Documents/gPPI")


#set up output structure
#rois <- list('A32sg','A32p','A24cd','A24rv','A14m','A11m','A13','A10m','A9m','A8m','A6m')
rois <- list('rACG','sgACC')
phases <- list('baseline','acquisition','early_extinction','extinction')
#phases <- list('baseline','acquisition','extinction')
groups <- list('healthy','ptsd')
conditions <- list('CS+','CS-')

lmm_res <- expand.grid(roi=rois,phase=phases,group=groups,condition=conditions)
lmm_res$beta <- 0
lmm_res$chisq <- 0
lmm_res$pval <- 0

#load the full dataset
df <- read.csv('full_ers_split_ext.csv')
df$subject <- as.factor(df$subject)
df$high_confidence_accuracy <- as.factor(df$high_confidence_accuracy)

for (GROUP in groups){
  for (PHASE in phases){
    for (ROI in rois){
      for (CON in conditions){
        #dat <- df[with(df,subject %notin% c(18,20,120) & group == GROUP & roi == ROI & encode_phase == PHASE & trial_type == CON),]
        dat <- df[with(df,group == GROUP & roi == ROI & encode_phase == PHASE & trial_type == CON),]
        wo_mem <- lmer(rsa ~ (1|subject), REML = FALSE, data = dat)
        w_mem <- lmer(rsa ~ (1|subject) + high_confidence_accuracy, REML = FALSE, data = dat)
        lrt_res <- anova(wo_mem,w_mem)
        lmm_summary <- summary(w_mem)
        
        lmm_res[with(lmm_res, group == GROUP & roi == ROI & phase == PHASE & condition == CON),'beta'] <- lmm_summary$coefficients[2]
        lmm_res[with(lmm_res, group == GROUP & roi == ROI & phase == PHASE & condition == CON),'chisq'] <- lrt_res$Chisq[2]
        lmm_res[with(lmm_res, group == GROUP & roi == ROI & phase == PHASE & condition == CON),'pval'] <- lrt_res$`Pr(>Chisq)`[2]
        
      }
    }
  }
}

lmm_res <- as.data.frame(lmm_res)
lmm_res <- apply(lmm_res,2,as.character)

write.csv(lmm_res,'rsa_lmm_results.csv')





##################
ROI <- 'A14m'
GROUP <- 'ptsd'
PHASE <- 'extinction'
CON <- 'CS+'



c <- df[with(df,group==GROUP & roi == ROI & encode_phase == PHASE & trial_type == CON),]
#c <- df[with(df,group==GROUP & roi == ROI & encode_phase == PHASE),]
#c$mem <- meanCenter(c$resp_num)

##################
#here down is actually correct
wo_mem <- lmer(rsa ~ (1|subject), REML = FALSE, data = c)
w_mem <- lmer(rsa ~ (1|subject) + high_confidence_accuracy, REML = FALSE, data = c)
anova(wo_mem,w_mem)

#maybe useful graphing stuff
#emmip(std,  high_confidence_accuracy)
#emmeans(std, pairwise ~ high_confidence_accuracy)

#c$preds = predict(std, c, type = "response")

#require(ggplot2)

#ggplot(c, aes(x=high_confidence_accuracy, y=preds)) +
#  geom_point(shape=16, cex=1) +
#  xlab("x1") + ylab("predicted y") +
#  theme(axis.title=element_text(size=26), axis.text=element_text(size=16))


require(permuco)

df <- read.csv('amyg_cem_conn.csv')
str(df)
#df$cope <- as.factor(df$cope)
df$group <- as.factor(df$group)
#df$seed <- as.factor(df$seed)
df$target <- as.factor(df$target)
df$subject <- as.factor(df$subject)

form <- conn ~ group* target + Error(subject / (target))
amyg.aov <- aovperm(form,df,np=10000)


df <- read.csv('split_out_conn.csv')
str(df)
df$condition <- as.factor(df$condition)
df$group <- as.factor(df$group)
df$seed <- as.factor(df$seed)
df$target <- as.factor(df$target)
df$subject <- as.factor(df$subject)
df$phase <- as.factor(df$phase)

sdf = df[with(df,seed == 'amyg_cem'),]
cres <- ezANOVA(sdf,dv=.(conn),wid=.(subject),within=.(condition,phase,target),between=.(group))


