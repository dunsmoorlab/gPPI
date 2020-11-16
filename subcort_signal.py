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
    df = pd.DataFrame({'beta':0.0},index=pd.MultiIndex.from_product(
                                [all_sub_args,rois,phase3,cons],
                                names=['subject','roi','phase','condition'])).sort_index()

    for sub in all_sub_args:
        print(sub)
        subj = bids_meta(sub)
        
        masks = {roi:get_data(f'{subj.masks}/{roi}.nii.gz') for roi in rois}
        betas = concat_imgs([nib.load(f'{subj.beta}/memory_run-0{run}_beta.nii.gz') for run in runs]).get_fdata()

        subdf = pd.read_csv(f'{subj.rsa}/fs_mask_roi_ER.csv')
        subdf = subdf[subdf.roi == 'sgACC'].reset_index(
            ).rename(columns={'index':'trial_num'}
            ).drop(columns=['roi','rsa']
            ).set_index(['encode_phase','trial_type']
            ).sort_index(
            ).dropna(subset=['response'])#sets us up to use .loc for stability

        for roi in rois:
            roi_data = apply_mask(mask=masks[roi],target=betas)
            for phase in phase3:
                for con in cons:
                    idx = subdf.loc[(phase,con),'trial_num'].values
                    df.loc[(sub,roi,phase,con),'beta'] = roi_data[idx].mean()
    df.to_csv('hc_head_subcortical_betas.csv')


'''graphing of the univariate data'''
betas = pd.read_csv('subcortical_betas.csv').set_index(['condition','group','roi','phase','subject']).sort_index()
betas = betas.loc['CS+'] - betas.loc['CS-']

rsa = pd.read_csv('pfc_ers_cleaned.csv').set_index(['condition','group','roi','phase','subject']).sort_index()
rsa = rsa.loc['CS+'] - rsa.loc['CS-']

for seed in rois:
    for target in ['vmPFC','dACC']:
        for phase in phase3:
            print(seed,target,phase)
            for group in ['healthy','ptsd']:
                print(pg.corr(betas.loc[(group,seed,phase),'beta'], rsa.loc[(group,target,phase),'rsa'])[['r','p-val']])
            print('\n\n')
    input()

# sns.catplot(data=df,x='roi',y='beta',hue='group',row='phase',col='condition',kind='bar')

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