import os
import sys

import pandas as pd
import numpy as np
import seaborn as sns

from wesanderson import wes_palettes
from matplotlib import rcParams
from wesanderson import wes_palettes
from numpy.random import RandomState
from sklearn.linear_model import LogisticRegression
from matplotlib.ticker import MultipleLocator
from scipy.special import expit
from matplotlib.patches import Patch
from nibabel.brikhead import parse_AFNI_header

sub_args = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,23,24,25,26]
p_sub_args = [101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117,118, 120, 121, 122, 123, 124, 125]
all_sub_args = sub_args + p_sub_args
smt_sub_args = [2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,101,102,103,104,105,106,107,108,109,110,111,112,113,114,116,117,118,120]
xcl_sub_args = [2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,19,21,101,102,103,104,105,106,107,108,109,110,111,112,113,114,116,117,118]

subjects = {'healthy':sub_args,
            'ptsd':p_sub_args,
            'all':all_sub_args}
groups = ['healthy','ptsd']
memcon = ['encoding','retrieval']
levels = ['item','set']
mems = ['hit','miss']
cons = ['CS+','CS-']
consp = ['CSp','CSm']
# rois = ['mOFC','dACC','amyg','hpc','ins','hc_head','hc_body','hc_tail','rh_hc_head','rh_hc_body','rh_hc_tail','lh_hc_head','lh_hc_body','lh_hc_tail','amyg_bla','amyg_cem','rh_amyg_bla','rh_amyg_cem','lh_amyg_bla','lh_amyg_cem']
# rois = ['mOFC','dACC','amyg','hpc','ins','lh_amyg','rh_amyg','lh_hpc','rh_hpc','sgACC','rACC','rSMA','rACG','hc_head','hc_body','hc_tail','rh_hc_head','rh_hc_body','rh_hc_tail','lh_hc_head','lh_hc_body','lh_hc_tail','amyg_bla','amyg_cem','rh_amyg_bla','rh_amyg_cem','lh_amyg_bla','lh_amyg_cem']  
phases = {'baseline':24,'acquisition':24,'early_extinction':8,'extinction':16}
# phases = {'baseline':24,'acquisition':24,'extinction':24}
subs = range(24)
conds = ['CSp','CSm']
bn_rois = ['A32sg','A32p','A24cd','A24rv','A14m','A11m','A13','A10m','A9m','A8m','A6m']
groups = ['healthy','ptsd']

gpal = list((wes_palettes['Zissou'][0],wes_palettes['Royal1'][1]))
cpal = ['darkorange','grey']
cpoint = sns.color_palette(cpal,n_colors=2,desat=.75)
spal = list((wes_palettes['Darjeeling1'][-1],wes_palettes['Darjeeling1'][0],wes_palettes['Darjeeling1'][1],))
phase_pal = sns.color_palette(['darkgrey','darkmagenta','seagreen'],desat=1)
sns.set_style('ticks', {'axes.spines.right':False, 'axes.spines.top':False})
# sns.set_style({'axes.facecolor':'.9','figure.facecolor':'.9'})
sns.set_context('notebook',font_scale=1.4)
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Arial'
rcParams['savefig.dpi'] = 300
rcParams['axes.titlepad'] = 30

WORK = '/work/05426/ach3377/lonestar/'
HOME = os.path.expanduser('~')
SCRATCH = '/scratch/05426/ach3377/'
gPPI_codebase = HOME + '/gPPI/'

def mkdir(path,local=False):
    if not local and not os.path.exists(path):
        os.makedirs(path)
def lgroup(x):
    if x > 100: return 'ptsd'
    else: return 'healthy'
def pfc_rename(x):
    if x == 'sgACC':
        return 'vmPFC'
    elif x == 'rACC':
        return 'dACC'
    else:
        return x
def amyg_rename(x):
    if x == 'amyg_bla':
        return 'Amyg. BLA'
    elif x == 'amyg_cem':
        return 'Amyg. CeM'
    else:
        return x
#these are BIDS-app made
bids_dir = os.path.join(SCRATCH,'fc-bids') if 'ach' not in sys.base_exec_prefix else os.path.join('/Volumes','DunsmoorRed','fc-bids')
deriv    = os.path.join(bids_dir, 'derivatives')
prep_dir = os.path.join(deriv,'fmriprep')
fs_dir   = os.path.join(deriv,'freesurfer')

#these are user made
model   = os.path.join(deriv,'model');#mkdir(model)
preproc = os.path.join(deriv,'preproc');#mkdir(preproc)
group_masks = os.path.join(deriv,'group_masks')



std_1mm_brain = os.path.join(WORK,'standard','MNI152_T1_1mm_brain.nii.gz')
std_3mm_brain = os.path.join(WORK,'standard','MNI152_T1_3mm_brain.nii.gz')
std_3mm_brain_mask = os.path.join(WORK,'standard','MNI152_T1_3mm_brain_mask.nii.gz')
std_2009_brain = os.path.join(SCRATCH,'standard','MNI152NLin2009cAsym_T1_1mm_brain.nii.gz')
std_2009_brain_mask = os.path.join(SCRATCH,'standard','MNI152NLin2009cAsym_T1_1mm_brain_mask.nii.gz')
std_2009_brain_3mm = os.path.join(SCRATCH,'standard','MNI152NLin2009cAsym_T1_3mm_brain.nii.gz')
std_2009_brain_mask_3mm = os.path.join(SCRATCH,'standard','MNI152NLin2009cAsym_T1_3mm_brain_mask.nii.gz')
gm_3mm_thr = os.path.join(SCRATCH,'standard','gm_3mm_thr.nii.gz')
gm_1mm_thr = os.path.join(SCRATCH,'standard','gm_1mm_thr.nii.gz')


tasks = {'baseline':{'n_trials':48,'ses':1,'n_tr':259},
         'acquisition':{'n_trials':48,'ses':1,'n_tr':259},
         'extinction':{'n_trials':48,'ses':1,'n_tr':259},
         'renewal':{'n_trials':24,'ses':2,'n_tr':135},
         'memory_run-01':{'n_trials':80,'ses':2,'n_tr':310},
         'memory_run-02':{'n_trials':80,'ses':2,'n_tr':310},
         'memory_run-03':{'n_trials':80,'ses':2,'n_tr':310},
         'localizer_run-01':{'n_trials':24,'ses':2,'n_tr':240},
         'localizer_run-02':{'n_trials':24,'ses':2,'n_tr':240}
         }

slices={'CS+':{
                 'baseline':{'encoding':slice(0,24),
                             'retrieval':slice(144,168)},
              
              'acquisition':{'encoding':slice(24,48),
                             'retrieval':slice(168,192)},

         'early_extinction':{'encoding':slice(48,56),
                            'retrieval':slice(192,200)},
               
               'extinction':{'encoding':slice(56,72),
                             'retrieval':slice(200,216)}},
        'CS-':{
                 'baseline':{'encoding':slice(72,96),
                            'retrieval':slice(216,240)},
              
              'acquisition':{'encoding':slice(96,120),
                            'retrieval':slice(240,264)},

         'early_extinction':{'encoding':slice(120,128),
                            'retrieval':slice(264,272)},
               
               'extinction':{'encoding':slice(128,144),
                            'retrieval':slice(272,288)}}}

slice3={'CS+':{
                 'baseline':{'encoding':slice(0,24),
                             'retrieval':slice(144,168)},
              
              'acquisition':{'encoding':slice(24,48),
                             'retrieval':slice(168,192)},
               
               'extinction':{'encoding':slice(48,72),
                             'retrieval':slice(192,216)}},
        'CS-':{
                 'baseline':{'encoding':slice(72,96),
                            'retrieval':slice(216,240)},
              
              'acquisition':{'encoding':slice(96,120),
                            'retrieval':slice(240,264)},
               
               'extinction':{'encoding':slice(120,144),
                            'retrieval':slice(264,288)}}}

split_slices={'CS+':{
                 'baseline':{'encoding':slice(0,24),
                             'retrieval':slice(144,168)},
              
        'early_acquisition':{'encoding':slice(24,36),
                            'retrieval':slice(168,180)},

         'late_acquisition':{'encoding':slice(36,48),
                            'retrieval':slice(180,192)},

         'early_extinction':{'encoding':slice(48,60),
                            'retrieval':slice(192,204)},
               
          'late_extinction':{'encoding':slice(60,72),
                            'retrieval':slice(204,216)}},
        
        'CS-':{
                 'baseline':{'encoding':slice(72,96),
                            'retrieval':slice(216,240)},
              
        'early_acquisition':{'encoding':slice(96,108),
                            'retrieval':slice(240,252)},

         'late_acquisition':{'encoding':slice(108,120),
                            'retrieval':slice(252,264)},

         'early_extinction':{'encoding':slice(120,132),
                            'retrieval':slice(264,276)},
               
          'late_extinction':{'encoding':slice(132,144),
                            'retrieval':slice(276,288)}}}


mem_slices = {'CS+':{
                 'baseline':slice(0,24),
              
              'acquisition':slice(24,48),

                'extinction':slice(48,72),

                     'foil':slice(72,120)},

        'CS-':{
                 'baseline':slice(120,144),
              
              'acquisition':slice(144,168),

               'extinction':slice(168,192),

                     'foil':slice(192,240)}}

def apply_mask(mask=None,target=None):

    coor = np.where(mask == 1)
    values = target[coor]
    if values.ndim > 1:
        values = np.transpose(values) #swap axes to get sample X feature
    return values

def p2z(p,tail):
    '''returns a z-score corresponding a to the pvalue, given the tail (1 or 2 sided)'''
    from scipy.stats import norm
    return norm.ppf(1-p*tail)


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
            
            self.ref2std      = os.path.join(self.reference,'ref2std.mat')
            self.std2ref      = os.path.join(self.reference,'std2ref.mat')

            self.ref2std3     = os.path.join(self.reference,'ref2std3.mat')
            self.std32ref     = os.path.join(self.reference,'std32ref.mat')


            self.ref2t1       = os.path.join(self.reference,'ref2t1.mat')
            self.t12std_nii   = os.path.join(self.reference,'t12std')
            self.t12std       = os.path.join(self.reference,'t12std.mat')
            self.t12std_warp  = os.path.join(self.reference,'t12std_warp')

            self.func = os.path.join(self.preproc_dir,'func');mkdir(self.func,local)
            self.beta = os.path.join(self.preproc_dir,'lss_betas');mkdir(self.beta,local) 

            self.fs_regmat = os.path.join(self.reference,'RegMat.dat')
            self.faa       = os.path.join(self.reference,'aparc+aseg.nii.gz')
            self.saa       = os.path.join(self.reference,'std_aparc+aseg.nii.gz')
            
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
