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

    std = nib.load(std_2009_brain_mask_3mm)
    for sub in all_sub_args:
        print(sub)
        subj = bids_meta(sub)
        out = f'{subj.rsa}/ers_sl_imgs'
        mkdir(out)

        with open(os.path.join(subj.rsa,'sl_er.p'),'rb') as file:
            mat = pickle.load(file)
        mat = new_img_like(std_2009_brain_3mm,mat)#need this for the inverse to work
        mat = masker.transform(mat)

        df = pd.read_csv(os.path.join(subj.rsa,'fs_mask_roi_ER.csv'))
        df = df[df.roi == 'sgACC'].reset_index(
            ).rename(columns={'index':'trial_num'}
            ).drop(columns=['roi','rsa']
            ).set_index(['encode_phase','trial_type']
            ).sort_index(
            ).dropna(subset=['response'])#sets us up to use .loc for stability

        for phase in phases:
            for con in conditions:
                est = mat[df.loc[(phase,con),'trial_num'].values,:].mean(axis=0)
                est = new_img_like(std,masker.inverse_transform(est).get_fdata(),copy_header=True)#need this for the header
                nib.save(est,f'{out}/{phase}_{conditions[con]}.nii.gz')

        #lets also do the subtraction i guess
        for phase in phases:
            csp = mat[df.loc[(phase,'CS+'),'trial_num'].values,:].mean(axis=0)
            csm = mat[df.loc[(phase,'CS-'),'trial_num'].values,:].mean(axis=0)
            est = csp - csm
            est = new_img_like(std,masker.inverse_transform(est).get_fdata(),copy_header=True)#need this for the header
            nib.save(est,f'{out}/{phase}_diff.nii.gz')

def one_samp_ttest(subs=None,phase=None,name=''):
    out_parent = '/scratch/05426/ach3377/searchlight/ers_comps'
    out_dir = f'/scratch/05426/ach3377/searchlight/ers_comps/{name}_{phase}'
    mkdir(out_dir)
    
    setA = ''
    for s, sub in enumerate(subs):
        setA += f'{bids_meta(sub).rsa}/ers_sl_imgs/{phase}_diff.nii.gz '
    
    n_cors = 'export OMP_NUM_THREADS=48'
    cd_cmd = f'cd {out_dir}'
    clustsim_cmd = f'3dttest++ -setA {setA} \
                               -Clustsim 48 \
                               -mask {gm_3mm_thr} \
                               -prefix {name}_{phase}_clst-ttest'
    
    script = f'{out_dir}/ttest_script.txt'
    os.system(f'rm {script}')
    
    for cmd in [n_cors, cd_cmd, clustsim_cmd]:
        os.system(f"echo {cmd} >> {script}")
    
    jobfile = f'/home1/05426/ach3377/gPPI/jobs/{name}_sl_ers_job.txt'
    os.system(f'rm {jobfile}')

    #not run here, just submiting a job
    os.system(f'echo singularity run --cleanenv \
                    /scratch/05426/ach3377/bids-apps/neurosft.simg \
                    bash -x {script} >> {jobfile}')

    os.system(f'launch -N 1 \
                       -n 1 \
                       -J 3dttest++ \
                       -s {jobfile} \
                       -m achennings@utexas.edu \
                       -p normal \
                       -r 0:05:00 \
                       -A LewPea_MRI_Analysis')
one_samp_ttest(subs=sub_args,phase='acquisition',name='healthy')
one_samp_ttest(subs=sub_args,phase='extinction',name='healthy')
one_samp_ttest(subs=p_sub_args,phase='acquisition',name='ptsd')
one_samp_ttest(subs=p_sub_args,phase='extinction',name='ptsd')

def paired_ttest(subs=None,name=''):
    out_parent = '/scratch/05426/ach3377/searchlight/ers_comps'
    out_dir = f'/scratch/05426/ach3377/searchlight/ers_comps/{name}'
    mkdir(out_dir)
    
    setA = ''
    setB = ''
    for s, sub in enumerate(subs):
        subj = bids_meta(sub)
        setA += f'{subj.rsa}/ers_sl_imgs/extinction_diff.nii.gz '
        setB += f'{subj.rsa}/ers_sl_imgs/acquisition_diff.nii.gz '

    n_cors = 'export OMP_NUM_THREADS=48'
    cd_cmd = f'cd {out_dir}'
    clustsim_cmd = f'3dttest++ -setA {setA} \
                               -setB {setB} \
                               -AminusB \
                               -paired \
                               -Clustsim 48 \
                               -mask {gm_3mm_thr} \
                               -prefix {name}_clst-ttest'
    
    script = f'{out_dir}/ttest_script.txt'
    os.system(f'rm {script}')
    
    for cmd in [n_cors, cd_cmd, clustsim_cmd]:
        os.system(f"echo {cmd} >> {script}")
    
    #run it
    os.system(f'singularity run --cleanenv \
                    /scratch/05426/ach3377/bids-apps/neurosft.simg \
                    bash -x {script}')
paired_ttest(subs=sub_args,name='healthy_phase_diff')
paired_ttest(subs=p_sub_args,name='ptsd_phase_diff')

def ind_ttest(phase=None,name=''):
    out_parent = '/scratch/05426/ach3377/searchlight/ers_comps'
    out_dir = f'/scratch/05426/ach3377/searchlight/ers_comps/{name}'
    mkdir(out_dir)
    
    setA = ''
    setB = ''
    for sub in sub_args:
        subj = bids_meta(sub)
        setA += f'{subj.rsa}/ers_sl_imgs/{phase}_diff.nii.gz '
    
    for sub in p_sub_args:
        subj = bids_meta(sub)
        setB += f'{subj.rsa}/ers_sl_imgs/{phase}_diff.nii.gz '
    
    n_cors = 'export OMP_NUM_THREADS=48'
    cd_cmd = f'cd {out_dir}'
    clustsim_cmd = f'3dttest++ -setA {setA} \
                               -setB {setB} \
                               -AminusB \
                               -Clustsim 48 \
                               -mask {gm_3mm_thr} \
                               -prefix {name}_clst-ttest'
    
    script = f'{out_dir}/ttest_script.txt'
    os.system(f'rm {script}')
    
    for cmd in [n_cors, cd_cmd, clustsim_cmd]:
        os.system(f"echo {cmd} >> {script}")
    
    #run it
    os.system(f'singularity run --cleanenv \
                    /scratch/05426/ach3377/bids-apps/neurosft.simg \
                    bash -x {script}')
ind_ttest(phase='acquisition',name='acq_group_diff')
ind_ttest(phase='extinction',name='ext_group_diff')

def ers_cluster(contrast=None,thr=0,nvox=0,mask='../../standard/gm_3mm_thr.nii.gz',tail=None):
        here = os.getcwd()
        folder = contrast
        name = contrast.split('/')[-1]
        os.chdir(folder)
        
        if tail == 'one-sided':
            side = '1sided RIGHT_TAIL'
        elif tail == 'two-sided':
            side = '2sided'

        cmap = f'{name}_ClusterMap.nii.gz';os.system(f'rm {cmap}')
        ceff = f'{name}_ClusterEffEst.nii.gz';os.system(f'rm {ceff}')

        ctxt = f'{name}_cluster.txt';os.system(f'rm {ctxt}')
        where = f'{name}_where.txt';os.system(f'rm {where}')
        
        cmd = f"3dClusterize -inset {name}_clst-ttest+orig \
                   -ithr 1 \
                   -idat 0 \
                   -mask {mask} \
                   -NN 3 \
                   -{side} p={thr} \
                   -clust_nvox {nvox} \
                   -pref_map {cmap} \
                   -pref_dat {ceff} > {ctxt}"

        
        os.system(cmd)
        
        if os.path.exists(f'{name}_ClusterMap.nii.gz'):
            w_cmd = f"whereami -coord_file {ctxt}'[1,2,3]' > {where}"
            os.system(w_cmd)
        

        os.chdir(here)

ers_cluster(contrast=f'{HOME}/Desktop/ers_comps/healthy_acquisition',thr=0.001,nvox=20,tail='one-sided')#20
ers_cluster(contrast=f'{HOME}/Desktop/ers_comps/healthy_extinction',thr=0.001,nvox=20,tail='one-sided')#20
ers_cluster(contrast=f'{HOME}/Desktop/ers_comps/ptsd_acquisition',thr=0.001,nvox=20,tail='one-sided')#20
ers_cluster(contrast=f'{HOME}/Desktop/ers_comps/ptsd_extinction',thr=0.001,nvox=21,tail='one-sided')#21

ers_cluster(contrast=f'{HOME}/Desktop/ers_comps/healthy_phase_diff',thr=0.001,nvox=16,tail='two-sided')
ers_cluster(contrast=f'{HOME}/Desktop/ers_comps/ptsd_phase_diff',thr=0.001,nvox=18,tail='two-sided')
ers_cluster(contrast=f'{HOME}/Desktop/ers_comps/acq_group_diff',thr=0.001,nvox=21,tail='two-sided')
ers_cluster(contrast=f'{HOME}/Desktop/ers_comps/ext_group_diff',thr=0.001,nvox=20,tail='two-sided')

# 3dClusterize -inset healthy_CSpE__CSpA_clst-ttest+tlrc -ithr 1 -idat 0 -mask /scratch/05426/ach3377/standard/gm_1mm_thr.nii.gz -NN 2 -2sided p=0.01 -clust_nvox 498 -pref_map healthy_CSpE__CSpA_ClusterMap.nii.gz -pref_dat healthy_CSpE__CSpA_ClusterEffEst.nii.gz > healthy_CSpE__CSpA_cluster.txt
