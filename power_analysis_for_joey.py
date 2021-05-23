import os
import pandas as pd
import numpy as np

from fg_config import *

for sub in sub_args:
    subj = bids_meta(sub)

    print(sub)
    for mask in ['caudate']:
        
        in_mask = f'/mnt/c/Users/ACH/Desktop/standard/meta_{mask}.nii.gz'
        out_mask = f'{subj.masks}/meta_{mask}.nii.gz'
        
        cmd = f'flirt -ref {subj.refvol_brain} \
                      -in {in_mask} \
                      -out {out_mask} \
                      -applyxfm \
                      -init {subj.std2ref} \
                      -interp nearestneighbour'
        
        os.system(cmd)

def apply_mask(mask=None,target=None):

    coor = np.where(mask == 1)
    values = target[coor]
    if values.ndim > 1:
        values = np.transpose(values) #swap axes to get feature X sample
    return values

from nilearn.image import get_data
from fg_config import *

rois = ['rACC','insula','thalamus','caudate']

df = pd.DataFrame({'beta':0.0},index=pd.MultiIndex.from_product([['CS+','CS-'],rois,sub_args],names=['condition','roi','subject']))

for sub in sub_args:
    subj = bids_meta(sub)
    print(sub)
    for condition, con in zip(['CS+','CS-'],['acquisition_CSp','acquisition_CSm']):
        img = get_data(f'{subj.weights}/{con}.nii.gz')
        
        for roi in rois:
            try:
                mask = get_data(f'{subj.masks}/meta_{roi}.nii.gz')
            except:
                mask = get_data(f'{subj.masks}/{roi}_mask.nii.gz')

            df.loc[(condition,roi,sub),'beta'] = apply_mask(mask=mask,target=img).mean() 

df = df.reset_index()
roi_rename = {'rACC':'dACC'}
df.roi = df.roi.apply(lambda x: roi_rename[x] if x in roi_rename else x)
df = df.set_index(['condition','subject']).sort_index()
df = df.reset_index()

for i in df.index:
    if df.loc[i,'roi'] == 'dACC' and df.loc[i,'subject'] == 26:
        df = df.drop(i)

df = df.set_index(['condition','subject']).sort_index()
results = df.groupby('roi').apply(lambda x: pg.ttest(x.loc[('CS+'),'beta'], x.loc[('CS-'),'beta'],paired=True))

print(results[['cohen-d']])
print(results['cohen-d'].apply(lambda x: pg.power_ttest(d=x,alpha=.05,power=.80,n=None,contrast='paired')))


df = df.reset_index().set_index(['condition','roi','subject'])
comp = df.loc['CS+'] - df.loc['CS-']

sns.set_context('notebook',font_scale=1.2)
fig, ax = plt.subplots()
comp = comp.reset_index()
sns.boxplot(data=comp,x='roi',y='beta', ax=ax)