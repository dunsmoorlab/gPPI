from fg_config import *
import nibabel as nib
from nilearn.image import get_data, concat_imgs
from roi_rsa import group_roi_rsa

def apply_mask(mask=None,target=None):

    coor = np.where(mask == 1)
    values = target[coor]
    if values.ndim > 1:
        values = np.transpose(values) #swap axes to get feature X sample
    return values

runs = [1,2,3]
rois = ['hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']
phase3 = ['baseline','acquisition','extinction']

def collect_univariate()
    dfs = {}
    df = pd.DataFrame({'beta':0.0},index=pd.MultiIndex.from_product(
                                [all_sub_args,phase3,cons],
                                names=['subject','phase','condition'])).sort_index()

    for sub in all_sub_args:
        print(sub)
        subj = bids_meta(sub)
        
        masks = {roi:get_data(f'{subj.masks}/{roi}.nii.gz') for roi in rois}
        betas = concat_imgs([nib.load(f'{subj.beta}/memory_run-0{run}_beta.nii.gz') for run in runs]).get_fdata()

        subdf = pd.read_csv(f'{subj.rsa}/fs_mask_roi_ER.csv')
        subdf = subdf[subdf.roi == 'sgACC'].reset_index(
            ).rename(columns={'index':'trial_num'}
            ).set_index(['encode_phase','trial_type']
            ).sort_index(
            ).dropna(subset=['response']#sets us up to use .loc for stability
            ).drop(columns=['roi','rsa','stimulus','memory_condition','low_confidence_accuracy','high_confidence_accuracy','phase','CSp_trial','CSm_trial','response'])

        for roi in rois:
            roi_data = apply_mask(mask=masks[roi],target=betas)
            subdf[roi] = roi_data[subdf.trial_num].mean(axis=1)
        dfs[sub] = subdf.reset_index()
            # for phase in phase3:
            #     for con in cons:
            #         idx = subdf.loc[(phase,con),'trial_num'].values
            #         df.loc[(sub,roi,phase,con),'beta'] = roi_data[idx].mean()
    # df.to_csv('subcortical_betas.csv')
        
    df = pd.concat(dfs.values())
def add_ers_to_full_df()
    df = pd.read_csv('subcort_betas_lmm.csv').set_index(['subject','encode_phase','trial_type','trial_num']).sort_index()
    dfs = {}
    for sub in all_sub_args:
        subj = bids_meta(sub)
        subdf = pd.read_csv(f'{subj.rsa}/fs_mask_roi_ER.csv')
        subdf = subdf[subdf.roi == 'sgACC'].reset_index(
            ).rename(columns={'index':'trial_num'}
            ).set_index('stimulus'
            ).sort_index().rename(columns={'rsa':'vmPFC_ers'})
        dacc = pd.read_csv(f'{subj.rsa}/fs_mask_roi_ER.csv')
        dacc = dacc[dacc.roi == 'rACC'].set_index('stimulus').sort_index().rename(columns={'rsa':'dACC_ers'})[['dACC_ers']]
        subdf = pd.concat((subdf,dacc),axis=1)

        subdf = subdf.reset_index(
            ).set_index(['subject','encode_phase','trial_type','trial_num']
            ).sort_index(
            ).dropna(subset=['response']#sets us up to use .loc for stability
            ).drop(columns=['roi','stimulus','memory_condition','low_confidence_accuracy','high_confidence_accuracy','phase','CSp_trial','CSm_trial','response'])
        dfs[sub] = subdf
    dfs = pd.concat(dfs.values())
    df = pd.concat((df,dfs),axis=1)
    df = df.reset_index().rename(columns={'encode_phase':'phase','trial_type':'condition'}).drop(columns='trial_num')
    df.to_csv('subcort_betas_lmm.csv',index=False)
    for sub in all_sub_args:
        df.loc[sub,('hc_tail','hc_body','hc_head','amyg_bla','amyg_cem','vmPFC_ers','dACC_ers')] = df.loc[sub,('hc_tail','hc_body','hc_head','amyg_bla','amyg_cem','vmPFC_ers','dACC_ers')].apply(zscore)
'''graphing of the univariate data'''
from robust_corr import *
from scipy.stats import zscore
df = pd.read_csv('ers_subcort_betas_diff.csv').set_index(['group','phase'])#,'subject'])
# for group in groups:
#     for phase in phase3:
#         df.loc[(group,phase),('dACC_ers','vmPFC_ers','hc_tail','hc_body','hc_head','amyg_bla','amyg_cem')] = df.loc[(group,phase),('dACC_ers','vmPFC_ers','hc_tail','hc_body','hc_head','amyg_bla','amyg_cem')].apply(zscore)
df = df.reset_index().set_index(['group','phase','subject'])
# betas = pd.read_csv('subcortical_betas.csv').set_index(['condition','group','roi','phase','subject']).sort_index()
# betas = betas.loc['CS+'] - betas.loc['CS-']
# sns.catplot(data=betas.reset_index(),x='roi',y='beta',hue='group',col='phase',kind='bar',col_order=phase3,order=rois)

# rsa = pd.read_csv('pfc_ers_cleaned.csv').set_index(['condition','group','roi','phase','subject']).sort_index()
# rsa = rsa.loc['CS+'] - rsa.loc['CS-']


for seed in rois:
    for target in ['vmPFC_ers','dACC_ers']:
        for phase in phase3:
            print(seed,target,phase)
            for group in ['healthy','ptsd']:
                # print(pg.corr(betas.loc[(group,seed,phase),'beta'], rsa.loc[(group,target,phase),'rsa'])[['r','p-val']])
                # skipped_corr(betas.loc[(group,seed,phase),'beta'], rsa.loc[(group,target,phase),'rsa'])
                c = skipped_corr(df.loc[('healthy',phase),seed],df.loc[('healthy',phase),target], return_dist = True)
                p = skipped_corr(df.loc[('ptsd',phase),seed],df.loc[('ptsd',phase),target], return_dist = True)
                diff_p = np.min(((1 - np.mean((c - p) > 0)) * 2,(1 - np.mean((p - c) > 0)) * 2))
                print(f'difference P = {diff_p}')
            print('\n\n')
    input()

plot_full_skipped_corr(df.loc[('healthy','extinction'),'hc_head'], df.loc[('healthy','extinction'),'vmPFC_ers'],'Healthy Hc Head predicts vmPFC Extinction ERS','Hc Head activity (CS+ > CS-)','vmPFC Extinction ERS (CS+ > CS-)')
plot_full_skipped_corr(df.loc[('ptsd','extinction'),'amyg_bla'], df.loc[('ptsd','extinction'),'dACC_ers'],'PTSD Amyg BLA predicts dACC Extinction ERS','Amyg BLA activity (CS+ > CS-)','dACC Extinction ERS (CS+ > CS-)')
plot_full_skipped_corr(df.loc[('ptsd','extinction'),'hc_tail'], df.loc[('ptsd','extinction'),'dACC_ers'],'PTSD Hc Tail predicts dACC Extinction ERS','Hc Tail activity (CS+ > CS-)','dACC Extinction ERS (CS+ > CS-)')


plot_full_skipped_corr(betas.loc[('ptsd','hc_tail','extinction'),'beta'], rsa.loc[('ptsd','dACC','extinction'),'rsa'],'ptsd hc_tail predict dACC extinction ERS')




plot_full_skipped_corr(betas.loc[('ptsd','hc_head','baseline'),'beta'], rsa.loc[('ptsd','dACC','baseline'),'rsa'],'ptsd hc_head predict dACC baseline ERS')
plot_full_skipped_corr(betas.loc[('ptsd','amyg_bla','extinction'),'beta'], rsa.loc[('ptsd','dACC','extinction'),'rsa'],'ptsd amyg_bla predict dACC extinction ERS')

#######looking at maybe CS+E vs. CS+A in the amygdala?
betas = pd.read_csv('subcortical_betas.csv').set_index(['condition','group','roi','phase','subject']).sort_index()
# betas = betas.loc['CS+'] - betas.loc['CS-']
# sns.catplot(data=betas.reset_index(),x='roi',y='beta',hue='group',col='phase',kind='bar',col_order=phase3,order=rois)

# rsa = pd.read_csv('pfc_ers_cleaned.csv').set_index(['condition','group','roi','phase','subject']).sort_index()
# rsa = rsa.loc['CS+'] - rsa.loc['CS-']
'''hairbrain idea to predict combo phase and condition using all these features'''
from sklearn.linear_model import LogisticRegression
df = pd.read_csv('subcort_betas_lmm.csv').set_index(['subject','condition','phase']).sort_index()
logreg = LogisticRegression(solver='lbfgs',multi_class='multinomial')
coefs = {}
# features = ['hc_tail','hc_body','hc_head','amyg_bla','amyg_cem','vmPFC_ers','dACC_ers']
features = ['amyg_bla','amyg_cem']
for sub in all_sub_args:
    logreg.fit(df.loc[(sub),(features)].values,df.loc[(sub),'phasecon'].values)
    coefs[sub] = pd.DataFrame({class_: logreg.coef_[c] for c, class_ in enumerate(logreg.classes_)},index=features)
    # coefs[sub] = pd.DataFrame(logreg.coef_.reshape(-1,1),columns=['beta'],index=features)
    
    coefs[sub]['subject']=sub
cdf = pd.concat(coefs.values()).reset_index().rename(columns={'index':'feature'}).melt(id_vars=['subject','feature'],var_name='phasecon',value_name='beta')
# cdf = pd.concat(coefs.values()).reset_index().rename(columns={'index':'feature'})

cdf['group'] = cdf.subject.apply(lgroup)

sns.catplot(data=cdf,x='feature',y='beta',hue='phasecon',row='group',kind='bar')
# sns.catplot(data=cdf,x='group',hue='feature',y='beta',kind='bar')


'''clean data for export to R'''
#first just the model of activity
betas = pd.read_csv('subcortical_betas.csv')
betas['group'] = betas.subject.apply(lgroup)
betas = betas.set_index(['condition','group','roi','phase','subject']).sort_index()
# betas.to_csv('subcortical_betas.csv')

#next ERS 
c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)

cdf = c.df.dropna(subset=['response'])
pdf = p.df.dropna(subset=['response'])
cdf = cdf.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = pdf.groupby(['trial_type','encode_phase','roi','subject']).mean()
rsa = pd.concat((cdf,pdf)).reset_index()
rsa = rsa[rsa.roi.isin(['sgACC','rACC'])]
rsa.roi = rsa.roi.apply(pfc_rename)
rsa['group'] = rsa.subject.apply(lgroup)
rsa = rsa.rename(columns={'trial_type':'condition','encode_phase':'phase'})
rsa = rsa.set_index(['condition','group','roi','phase','subject']).sort_index()[['rsa']]
# rsa.to_csv('pfc_ers_cleaned.csv')

#combine betas and ers
betas = betas.unstack(level='roi')
betas.columns = betas.columns.droplevel(0)
betas = betas.rename_axis(None, axis=1)

rsa = rsa.unstack(level='roi')
rsa.columns = rsa.columns.droplevel(0)
rsa = rsa.rename_axis(None, axis=1).rename(columns={'dACC':'dACC_ers','vmPFC':'vmPFC_ers'})

df = pd.concat((rsa,betas),axis=1)
df.to_csv('ers_subcort_betas_full.csv')