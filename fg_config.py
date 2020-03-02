import os
import sys

import pandas as pd
import numpy as np
import seaborn as sns

sub_args = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,23,24,25,26]
p_sub_args = [101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117,118, 120, 121, 122, 123, 124, 125]
all_sub_args = sub_args + p_sub_args

subjects = {'control':sub_args,
            'ptsd':p_sub_args,
            'all':all_sub_args}
# gpal = list((wes_palettes['Zissou'][0],wes_palettes['Royal1'][1]))
cpal = ['darkorange','grey']
cpoint = sns.color_palette(cpal,n_colors=2,desat=.75)
WORK = '/work/05426/ach3377/lonestar/'
HOME = '/home1/05426/ach3377/'
SCRATCH = '/scratch/05426/ach3377/'
gPPI_codebase = HOME + 'gPPI/'

def mkdir(path,local=False):
    if not local and not os.path.exists(path):
        os.makedirs(path)
def lgroup(x):
    if x > 100: return 'ptsd'
    else: return 'control'

#these are BIDS-app made
bids_dir = os.path.join(SCRATCH,'fc-bids')
deriv    = os.path.join(bids_dir, 'derivatives')
prep_dir = os.path.join(deriv,'fmriprep')
fs_dir   = os.path.join(deriv,'freesurfer')

#these are user made
model   = os.path.join(deriv,'model');#mkdir(model)
preproc = os.path.join(deriv,'preproc');#mkdir(preproc)



std_1mm_brain = os.path.join(WORK,'standard','MNI152_T1_1mm_brain.nii.gz')
std_3mm_brain = os.path.join(WORK,'standard','MNI152_T1_3mm_brain.nii.gz')
std_3mm_brain_mask = os.path.join(WORK,'standard','MNI152_T1_3mm_brain_mask.nii.gz')

tasks = {'baseline':{'n_trials':48,'ses':1},
         'acquisition':{'n_trials':48,'ses':1},
         'extinction':{'n_trials':48,'ses':1},
         'renewal':{'n_trials':24,'ses':2},
         'memory_run-01':{'n_trials':80,'ses':2},
         'memory_run-02':{'n_trials':80,'ses':2},
         'memory_run-03':{'n_trials':80,'ses':2},
         'localizer_run-01':{'n_trials':24,'ses':2},
         'localizer_run-02':{'n_trials':24,'ses':2}
         }


class bids_meta(object):

    def __init__(self, sub):

        if sub == 'gusbrain':
            self.fs_id = 'gusbrain'
            self.fsub = 'gusbrain'
        elif sub == 'fs':
            self.fs_id = 'fsaverage'
            self.fsub = 'fsaverage'
        else:
            
            local = False
            if 'ach' in sys.base_exec_prefix:
                local = True
            
            self.num = int(sub)
            
            self.fsub = 'sub-FC{0:0=3d}'.format(self.num)

            self.subj_dir = os.path.join(bids_dir,self.fsub)
            self.prep_dir = os.path.join(prep_dir,self.fsub)
            self.fs_dir   = os.path.join(fs_dir,self.fsub)

            self.model_dir   = os.path.join(model,self.fsub);mkdir(self.model_dir,local)
            self.feat_dir    = os.path.join(self.model_dir,'feats');mkdir(self.feat_dir,local)
            self.preproc_dir = os.path.join(preproc,self.fsub);mkdir(self.preproc_dir,local)
            
            self.reference    = os.path.join(self.preproc_dir,'reference');mkdir(self.reference,local)
            self.t1           = os.path.join(self.reference,'T1w.nii.gz')
            self.t1_mask      = os.path.join(self.reference,'T1w_mask.nii.gz')
            self.t1_brain     = os.path.join(self.reference,'T1w_brain.nii.gz')
            self.refvol       = os.path.join(self.reference,'boldref.nii.gz')
            self.refvol_mask  = os.path.join(self.reference,'boldref_mask.nii.gz')
            self.refvol_brain = os.path.join(self.reference,'boldref_brain.nii.gz')

            self.func = os.path.join(self.preproc_dir,'func');mkdir(self.func,local)
            self.beta = os.path.join(self.preproc_dir,'lss_betas');mkdir(self.beta,local) 

            self.fs_regmat = os.path.join(self.reference,'RegMat.dat')
            self.faa       = os.path.join(self.reference,'aparc+aseg.nii.gz')
            
            self.masks  = os.path.join(self.preproc_dir,'masks');mkdir(self.masks,local)
            self.weights = os.path.join(self.preproc_dir,'rsa_weights');mkdir(self.weights,local)

            self.rsa = os.path.join(self.model_dir,'rsa_results');mkdir(self.rsa,local)
    def cs_lookup(self):    
        if self.meta['DataFile.Basename'][0][0] == 'A':
            self.csplus = 'animal'
            self.csminus = 'tool'
        elif self.meta['DataFile.Basename'][0][0] == 'T':
            self.csplus = 'tool'
            self.csminus = 'animal'
