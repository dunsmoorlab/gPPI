import os
import pickle
import numpy as np
import pandas as pd
import nibabel as nib

from fg_config import *
from bids_model import bids_events

from nilearn.input_data import NiftiMasker as nm
from nilearn.image import new_img_like
from collections import OrderedDict
from scipy.stats import ttest_1samp, ttest_ind


conditions = {'CS+': 'CSp',
              'CS-': 'CSm'}
groups = ['healthy','ptsd']
phase3 = ['baseline','acquisition','extinction']
phases = ['baseline','acquisition','early_extinction','extinction']
masker = nm(mask_img=std_2009_brain_mask_3mm)
masker.fit()

mats = {}
mem_
for phase in phases: mats[phase] = np.zeros((len(all_sub_args),69880))

for s, sub in enumerate(all_sub_args):
    print(sub)
    subj = bids_meta(sub)

    with open(os.path.join(subj.rsa,'sl_er.p'),'rb') as file:
        mat = pickle.load(file)

    mat = new_img_like(std_2009_brain_3mm,mat)
    mat = masker.transform(mat)

    df = pd.read_csv(os.path.join(subj.rsa,'fs_mask_roi_ER.csv'))
    df = df[df.roi == 'mOFC']
    df = df.drop(columns=['roi','rsa'])
    for cs in conditions:
        con = '%s_trial'%(conditions[cs])
        for i in range(1,9): 
            df.loc[ df[ df.encode_phase == 'extinction' ][ df[con] == i ].index,'encode_phase' ] = 'early_extinction'

    for phase in phases:
        csp = mat[df[df.encode_phase == phase][df.trial_type == 'CS+'].index,:].mean(axis=0)
        csm = mat[df[df.encode_phase == phase][df.trial_type == 'CS-'].index,:].mean(axis=0)
        mats[phase][s,:] = csp - csm

def arr_ttest_1samp(a):
    t, p = ttest_1samp(a,0)
    return np.array([t, p])

def arr_ttest_ind(a):
    t, p = ttest_ind(a[:24], a[24:])
    return np.array([t, p])

res = {}
for group in groups: res[group] = {}
res['comp'] = {}

for phase in phases:
    print(phase)
    res['healthy'][phase] = np.apply_along_axis(arr_ttest_1samp,0,mats[phase][:24])
    print('healthy done')
    res['ptsd'][phase] = np.apply_along_axis(arr_ttest_1samp,0,mats[phase][24:])
    print('ptsd done')
    res['comp'][phase] = np.apply_along_axis(arr_ttest_ind,0,mats[phase])

sl_dir = os.path.join(SCRATCH,'searchlight')
with open(os.path.join(sl_dir,'group_results_ttests.p'),'wb') as file:
    pickle.dump(res,file)

for phase in phases:
    for group in ['healthy','ptsd','comp']:
        nib.save(masker.inverse_transform(res[group][phase][0]), os.path.join(sl_dir,'%s_%s_tmap.nii.gz'%(group,phase)))
        nib.save(masker.inverse_transform(1 - res[group][phase][1]), os.path.join(sl_dir,'%s_%s_pmap.nii.gz'%(group,phase)))