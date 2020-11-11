from fg_config import *
import nibabel as nib
from nilearn.image import get_data, concat_imgs

def apply_mask(mask=None,target=None):

    coor = np.where(mask == 1)
    values = target[coor]
    if values.ndim > 1:
        values = np.transpose(values) #swap axes to get feature X sample
    return values

runs = [1,2,3]
# rois = ['hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']
rois = ['hc_head']
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
df = pd.read_csv('subcortical_betas.csv')
df['group'] = df.subject.apply(lgroup)

sns.catplot(data=df,x='roi',y='beta',hue='group',row='phase',col='condition',kind='bar')


