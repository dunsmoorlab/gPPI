import os
import pickle
import numpy as np
import pandas as pd
import nibabel as nib

from fg_config import *
from bids_model import bids_events

from nilearn.input_data import NiftiMasker
from nilearn.image import new_img_like
from collections import OrderedDict
from scipy.stats import ttest_1samp, ttest_ind, wilcoxon


conditions = {'CS+': 'CSp',
              'CS-': 'CSm'}
groups = ['healthy','ptsd']
phase3 = ['baseline','acquisition','extinction']
phases = ['baseline','acquisition','early_extinction','extinction']
masker = NiftiMasker(mask_img=std_2009_brain_mask_3mm)
masker.fit()

mats = {}
# mem_
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
    # for cs in conditions:
    #     con = '%s_trial'%(conditions[cs])
    #     for i in range(1,9): 
    #         df.loc[ df[ df.encode_phase == 'extinction' ][ df[con] == i ].index,'encode_phase' ] = 'early_extinction'

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


'''FEARMEM SEARCHLIGHT ANALYSES'''
'''ERS'''
sl_dir = os.path.join(SCRATCH,'searchlight')
conditions = {'CS+': 'CSp',
              'CS-': 'CSm'}
phases = ['baseline','acquisition','extinction']

def aggregate_sm_ers():
    masker = NiftiMasker(mask_img=std_2009_brain_mask_3mm)
    masker.fit()

    sm = pd.read_csv('sm_by_stim_ref.csv').set_index(['subject','stimulus'])

    mats = {}
    for encode in phases:
        mats[encode] = {}
        for con in conditions:
            mats[encode][con] = {}
            for response in phases:
                mats[encode][con][response] = {}

    for s, sub in enumerate(xcl_sub_args):
        print(sub)
        subj = bids_meta(sub)

        with open(os.path.join(subj.rsa,'sl_er.p'),'rb') as file:
            mat = pickle.load(file)

        mat = new_img_like(std_2009_brain_3mm,mat)
        mat = masker.transform(mat)

        sub_df = pd.read_csv(f'{subj.rsa}/fs_mask_roi_ER.csv')
        sub_df = sub_df[sub_df.roi == 'sgACC']
        sub_df = sub_df.drop(columns=['roi','rsa']).reset_index(drop=True)
        sub_df['source_memory'] = ''
        for i, stim in enumerate(sub_df.stimulus):
            sub_df.loc[i,'source_memory'] = sm.loc[(sub,stim),'source_memory']

        for encode in phases:
            for con in cons:
                for response in phases:
                     
                    _slice = sub_df[sub_df.encode_phase == encode][sub_df.trial_type == con][sub_df.source_memory == response].index
                    if _slice.size == 0:
                        pass
                    else:
                        mats[encode][con][response][sub] = mat[_slice,:].mean(axis=0)
    
    with open(os.path.join(sl_dir,'sm_ers_data.p'),'wb') as file:
        pickle.dump(mats,file)

def arr_wilcoxon(a,n,tail):
    w, p = wilcoxon(a[:n], a[n:],alternative=tail, zero_method='wilcox', correction=True)
    return np.array([w, p])

def ers_contrasts(con1=[],con2=[],tail='',mats={},masker=None):
    data1 = mats[con1[0]][con1[1]][con1[2]]
    data2 = mats[con2[0]][con2[1]][con2[2]]
    _subs = [s for s in xcl_sub_args if s in data1.keys() and s in data2.keys()]
    
    _size = len(_subs)
    pe = np.zeros((_size*2,69880))

    for s, sub in enumerate(_subs):
        pe[s,:] = data1[sub]
        pe[_size+s,:] = data2[sub]

    cope = np.apply_along_axis(arr_wilcoxon,0,pe,_size,tail)

    Cons = {'+':'p','-':'m'}
    name = f'CS{Cons[con1[1][-1]]}{con1[0][0].upper()}_{con1[2][0].upper()}__CS{Cons[con2[1][-1]]}{con2[0][0].upper()}_{con2[2][0].upper()}'

    nib.save(masker.inverse_transform(1 - cope[1]), f'{sl_dir}/{name}_1-pmap.nii.gz')


def afni_3dttest(con1=[],con2=[],tail='',mats={},masker=None,std=nib.load(std_2009_brain_mask_3mm)):
    data1 = mats[con1[0]][con1[1]][con1[2]]
    data2 = mats[con2[0]][con2[1]][con2[2]]
    _subs = [s for s in xcl_sub_args if s in data1.keys() and s in data2.keys()]

    Cons = {'+':'p','-':'m'}
    name1 = f'CS{Cons[con1[1][-1]]}{con1[0][0].upper()}_{con1[2][0].upper()}'
    name2 = f'CS{Cons[con2[1][-1]]}{con2[0][0].upper()}_{con2[2][0].upper()}'
    name = f'{name1}__{name2}'

    out_dir = f'{sl_dir}/{name}'
    mkdir(out_dir)

    setA = ''
    setB = ''

    for s, sub in enumerate(_subs):
        file1 = f'{out_dir}/{sub}_{name1}.nii.gz'
        file2 = f'{out_dir}/{sub}_{name1}.nii.gz'

        nib.save(new_img_like(std,masker.inverse_transform(data1[sub]).get_fdata(),copy_header=True), file1)
        nib.save(new_img_like(std,masker.inverse_transform(data2[sub]).get_fdata(),copy_header=True), file2)

        setA += f'{file1} '
        setB += f'{file1} '

    cmd = f'3dttest++ -setA {setA} -setB {setB} -paired -prefix {out_dir}/{name}_ttest'
    os.system(cmd)

with open(f'{sl_dir}/sm_ers_data.p','rb') as file:
    mats = pickle.load(file)
masker = NiftiMasker(mask_img=std_2009_brain_mask_3mm)
masker.fit()

afni_3dttest(con1=['baseline','CS+','baseline'],
             con2=['baseline','CS+','acquisition'],
             tail='greater',mats=mats,masker=masker)
afni_3dttest(con1=['baseline','CS-','acquisition'],
             con2=['baseline','CS+','acquisition'],
             tail='greater',mats=mats,masker=masker)


ers_contrasts(con1=['baseline','CS+','baseline'],
             con2=['baseline','CS+','acquisition'],
             tail='greater',mats=mats,masker=masker)

ers_contrasts(con1=['baseline','CS-','acquisition'],
             con2=['baseline','CS+','acquisition'],
             tail='greater',mats=mats,masker=masker)