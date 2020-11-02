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
phases = ['acquisition','extinction']

masker = NiftiMasker(mask_img=std_2009_brain_mask_3mm)
masker.fit()


def sub_imgs():
    mats = {}
    for phase in phases: mats[phase] = np.zeros((len(all_sub_args),69880))
    for s, sub in enumerate(all_sub_args):
        print(sub)
        subj = bids_meta(sub)
        
        with open(os.path.join(subj.rsa,'sl_er.p'),'rb') as file:
            mat = pickle.load(file)
        mat = new_img_like(std_2009_brain_3mm,mat)
        mat = masker.transform(mat)

        df = pd.read_csv(os.path.join(subj.rsa,'fs_mask_roi_ER.csv'))
        df = df[df.roi == 'sgACC'].reset_index(
            ).rename(columns={'index':'trial_num'}
            ).drop(columns=['roi','rsa']
            ).set_index(['encode_phase','trial_type']
            ).sort_index(
            ).dropna(subset=['response'])

        for phase in phases:
            for con in conditions:
                est = mat[df.loc[(phase,con),'trial_num'].values,:].mean(axis=0)
                
