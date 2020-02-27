import os
import re

import pandas as pd
import nibabel as nib

from collections import OrderedDict
from nilearn.image import concat_imgs, get_data
from nilearn.signal import clean
from sklearn.linear_model import Ridge
from scipy.signal import convolve

from fg_config import *
from brainiak_utils import _double_gamma_hrf, convolve_hrf

class bids_events():
    
    def __init__(self,sub):

        self.subj = bids_meta(sub)
        # self.fsl_events()
        # self.confounds()
    
    #generate CS+, CS-, US timing files for use in FSL
    def fsl_events(self):

        outstr = {'CS+'      :'CSp',
                  'CS-'      :'CSm',
                  'animal'   :'Animal',
                  'tool'     :'Tool',
                  'indoor'   :'Indoor',
                  'outdoor'  :'Outdoor',
                  'scrambled':'Scrambled',
                  'rest'     :'Rest',
                  'US'       :'US'}

        #walk through every folder containing the raw data and find all the event files
        for folder in os.walk(self.subj.subj_dir):
            for file in folder[2]:
                if 'events' in file and '.tsv' in file:
                    events = pd.read_csv(os.path.join(self.subj.subj_dir,folder[0],file), sep='\t')
                    phase = re.search('task-(.*)_events',file)[1]
                    out = os.path.join(self.subj.model_dir,'%s'%(phase))
                    mkdir(out)

                    #for every trial type
                    for con in events.trial_type.unique():
                        con_timing = events[events.trial_type == con][['onset','duration']]
                        con_timing['PM'] = 1
                        con_timing.to_csv( os.path.join(out, '%s_all.txt'%(outstr[con])),
                            sep='\t', float_format='%.8e', index=False, header=False)

                    #need to generate the US timing file
                    if phase == 'acquisition':
                        con = 'US'
                        US = events[events.shock == 'CSUS'][['onset','duration']].copy()
                        US.onset += US.duration
                        US.duration = 0
                        US['PM'] = 1
                        US.to_csv( os.path.join(out, '%s_all.txt'%(outstr[con])),
                            sep='\t', float_format='%.8e', index=False, header=False)

                    #handle early/late for associative learning phases

    #collect confound regressors from fMRIprep
    def confounds(self):

        #confounds of interest
        COI = ['a_comp_cor_00','framewise_displacement','trans_x','trans_y','trans_z','rot_x','rot_y','rot_z']
        # COI = ['global_signal','csf','white_matter','framewise_displacement','trans_x','trans_y','trans_z','rot_x','rot_y','rot_z']
        
        #walk through every folder of fMRIprep output and find all the confound files
        for folder in os.walk(self.subj.prep_dir):
            for file in folder[2]:
                if 'confounds' in file and '.tsv' in file:
                    C = pd.read_csv(os.path.join(self.subj.prep_dir,folder[0],file), sep='\t')
                    run_COI = COI.copy()
                    for _c in C.columns:
                        if 'cosine' in _c or 'motion_outlier' in _c:
                                run_COI.append(_c)
                    C = C[run_COI]
                    C['constant'] = 1
                    C['framewise_displacement'][0] = 0
                    
                    phase = re.search('task-(.*)_desc',file)[1]
                    out = os.path.join(self.subj.model_dir,'%s'%(phase))
                    C.to_csv(os.path.join(out,'confounds.txt'),
                        sep='\t',float_format='%.8e', index=False, header=False)

class lss():
    def __init__(self,sub):
        self.subj = bids_meta(sub)
    
    def estimate(self):
        for folder in os.walk(self.subj.subj_dir):
            for file in folder[2]:
                if 'events' in file and '.tsv' in file:
                    events = pd.read_csv(os.path.join(self.subj.subj_dir,folder[0],file), sep='\t')
                    phase = re.search('task-(.*)_events',file)[1]
                    lss_dir = os.path.join(self.subj.model_dir,'%s'%(phase),'lss_betas')
                    mkdir(lss_dir)

                    trial_types = events.trial_type.unique()

                    for trial in range(events.shape[0]):
                        beta_folder = os.path.join(lss_dir,'trial_{0:0=2d}'.format(trial))
                        mkdir(beta_folder)
                        beta_trial = events.loc[trial,['onset','duration']]
                        beta_trial['PM'] = 1
                        beta_trial.to_csv(os.path.join(beta_folder,'beta.txt'),
                                    sep='\t', float_format='%.8e', index=False, header=False)

                        for i, condition in enumerate(trial_types):
                            #grab all trials of each condition
                            con_trials = np.where(events.trial_type == condition)[0]
                            #make sure we don't model the same trial twice
                            con_trials = [t for t in con_trials if t != trial]
                            #get the onsets/durations
                            con_events = events.loc[con_trials,['onset','duration']]
                            con_events['PM'] = 1
                            con_events.to_csv(os.path.join(beta_folder,'no_interest_%s.txt'%(i)),
                                            sep='\t', float_format='%.8e', index=False, header=False)

                    self._autofill_lss(lss_dir=lss_dir,phase=phase,n_trials=events.shape[0])

    def _autofill_lss(self,lss_dir,phase,n_trials):
            template = os.path.join(gPPI,'feats','template_lss_%s.fsf'%(phase))

            replacements = {'SUBID':self.subj.fsub}            
            #need to handle the special cases where the TR is longer
            if self.subj.num in [105,106]:
                replacements['TR_length'] = '2.23'
            else:
                replacements['TR_length'] = '2'

            for t in range(n_trials):
                trial = 'trial_{0:0=2d}'.format(t)
                replacements['TRIAL'] = trial
                
                beta_folder = os.path.join(lss_dir,trial)
                outfeat = os.path.join(beta_folder,'%s.fsf'%(trial))

                with open(template) as infile: 
                    with open(outfeat, 'w') as outfile:
                        for line in infile:
                            for src, target in replacements.items():
                                line = line.replace(src, target)
                            outfile.write(line)

                with open('jobs/lss_betas/%s_%s_job.txt'%(self.subj.fsub,phase),'w') as txt_file:
                    txt_file.write('feat %s'%(outfeat))
                #also go ahead and make the job script here
                os.system('echo "feat %s" >> jobs/lss_betas/%s_job.txt'%(outfeat,self.subj.fsub))

    def reconstruct(self):
        for task in tasks:
            lss_dir = os.path.join(self.subj.model_dir,task,'lss_betas')
            beta_fname = os.path.join(self.subj.beta,'%s_beta.nii.gz'%(task))
            
            beta_img = OrderedDict()
            concat_ = True
            for i in range(tasks[task]['n_trials']):
                trial = 'trial_{0:0=2d}'.format(i)
                try:
                    beta_img[i] = nib.load(os.path.join(lss_dir,trial,'%s.feat'%(trial),'stats','cope1.nii.gz'))
                except FileNotFoundError: 
                    if os.path.exists(os.path.join(lss_dir,trial+'+.feat')):
                        os.system('rm -r %s'%(os.path.join(lss_dir,trial,trial+'.feat')))
                        os.system('mv %s %s'%( os.path.join(lss_dir,trial+'+.feat'), os.path.join( lss_dir,trial,'%s.feat'%(trial) ) ) )    
                    else:
                        os.system('mv %s %s'%(os.path.join(lss_dir,trial+'.feat'),os.path.join(lss_dir,trial+'/') ))
                    try:
                        beta_img[i] = nib.load(os.path.join(lss_dir,trial,'%s.feat'%(trial),'stats','cope1.nii.gz'))
                    except: 
                        os.system('echo "%s" >> bad_lss.txt'%(os.path.join(self.subj.model_dir,task,'lss_betas',trial)))
                        concat_ = False

            if concat_:
                # concatenate them
                beta_img = concat_imgs(beta_img.values())
                nib.save(beta_img,beta_fname)

            #mask them too
            os.system('fslmaths %s -mas %s %s'%(beta_fname,self.subj.refvol_mask,beta_fname))

class gPPI():
    
    def __init__(self,sub,mask=None,phases=None):

        self.subj = bids_meta(sub)
        self.mask = self.load_mask(mask)
        self.data = self.load_clean_data(phases=phases)
        self.extract_timecourse(mask_name=mask)

    def load_mask(self,mask):
        
        if mask is not None:
            if 'mask' in mask or 'lh_' in mask or 'rh_' in mask: 
                fname = os.path.join(self.subj.masks,'%s.nii.gz'%(mask))
            else:
                fname = os.path.join(self.subj.masks,'%s_mask.nii.gz')

            return get_data(fname)
    
    def _apply_mask(self,mask=None,target=None):

        coor = np.where(mask == 1)
        values = target[coor]
        if values.ndim > 1:
            values = np.transpose(values) #swap axes to get feature X sample
        return values
    
    def load_clean_data(self,phases=None):
        if phases is 'all':
            phases = tasks
        elif type(phases) is str:
            phases = [phases]

        #load the data
        data = {phase:get_data(os.path.join(self.subj.func,file)) for phase in phases for file in os.listdir(self.subj.func) if phase in file}
        #mask the data and take mean timeseries
        data = {phase: self._apply_mask(mask=self.mask,target=data[phase]).mean(axis=1) for phase in data}
        #clean the data 
        data = {phase: clean(data[phase][:,np.newaxis], #need this here to be feature X sample after meaning
                        
                        confounds=pd.read_csv(os.path.join(self.subj.model_dir,phase,'confounds.txt'),
                                        sep='\t',header=None).values,
                
                        t_r=2,detrend=False,standardize='zscore')
                                                            for phase in data}

        return data

    #extract the givin timecourse for each run
    def extract_timecourse(self,mask_name=None): 
        #deconvolve
        self.neuronal = {phase: self._deconvolve(self.data[phase]) for phase in self.data}

        for phase in self.neuronal:
                df = pd.Series(self.neuronal[phase])
                df.to_csv(os.path.join(self.subj.model_dir,phase,'%s_neuronal_signal.txt'%(mask_name)),
                    sep='\t', float_format='%.8e', index=False, header=False)
        

    def _deconvolve(self,dat):
        #create an HRF in unit AU (0-1)
        hrf = np.array(_double_gamma_hrf(temporal_resolution=.5))
        hrf /= np.max(hrf)

        #how long the data is
        N = dat.shape[0] 
        
        #create basis set
        xb = self._dct_mat(N,N) 
        
        #convolve it with HRF
        Hxb = np.zeros((N,N)) 
        for i in range(N):
            Hx = convolve(xb[:,i],hrf)[:N]
            Hxb[:,i] = Hx
        
        #initialize the ridge regression
        reg = Ridge(alpha=1,solver='lsqr',fit_intercept=False,normalize=False,max_iter=1000)
        #fit it
        reg.fit(Hxb,dat)
        
        #transforms neuronal signal out of cosine space
        neuronal = np.matmul(xb,reg.coef_[0,:])
        # neuronal = [neuronal,np.newaxis]
        return neuronal

    #basis set function
    def _dct_mat(self,N_,K_):
        n = np.array((range(0,N_))).T
        C_ = np.zeros((n.shape[0],K_))
        C_[:,0] = np.ones(n.shape[0])/np.sqrt(N_)
        for q in range(1,K_):
            C_[:,q] = np.sqrt(2/N_)* np.cos( np.pi*(2*n)* (q) /(2*N_))
        return C_

def autofill_fsf(template='',ses=None):
    outstr = re.search('template_(.*)',template)[1]
    for sub in all_sub_args:
        subj = bids_meta(sub)
        replacements = {'SUBID':subj.fsub}
        
        #need to handle the special cases where the TR is longer
        if ses == 1 and sub in [105,106]:
            replacements['TR_length'] = '2.23'
        else:
            replacements['TR_length'] = '2'

        outfeat = os.path.join(subj.feat_dir,'%s_%s.fsf'%(subj.fsub,outstr))

        with open(os.path.join(gPPI,'feats','%s.fsf'%(template))) as infile: 
            with open(outfeat, 'w') as outfile:
                for line in infile:
                    for src, target in replacements.items():
                        line = line.replace(src, target)
                    outfile.write(line)

        #also go ahead and make the job script here
        os.system('echo "feat %s" >> jobs/%s_job.txt'%(outfeat,outstr))

def wrap_lss_jobs():
    for i, job in enumerate(os.listdir('jobs/lss_betas')):
        if '.txt' in job:
            os.system('launch -N 1 -n 12 -J lss_%s -s jobs/lss_betas/%s -m achennings@utexas.edu -p normal -r 09:00:00 -A fMRI-Fear-Conditioni'%(i,job))

    for job in ['acquisition','extinction','baseline']:
        os.system('launch -N 1 -n 12 -J %s -s jobs/%s_rsa_job.txt -m achennings@utexas.edu -p normal -r 12:00:00')

    #submit a bunch of jobs at once
    for job in range(28):
        os.system('launch -N 1 -n 12 -J flss_%s -s jobs/final_lss_job_%s.txt -m achennings@utexas.edu -p normal -r 12:00:00 -A LewPea_MRI_Analysis'%(job,job))

    #splitting up a bunch of jobs into different job scripts
    for i in range(bad.shape[0]):
        os.system('echo "%s" >> jobs/final_lss_job_%s.txt'%(bad[0][i],int(np.floor(i/12))))

def clean_bad_lss():
    bad = pd.read_csv('bad_lss.txt',header=None)
    for beta in bad[0]:
        trial = beta[-8:]
        for beta_path in [beta+'.feat', beta+'/'+trial+'.feat', beta+'+.feat', beta+'/'+trial+'+.feat']:
            if os.path.exists(beta_path):
                # print(os.path.exists(beta_path))
                os.system('rm -r %s'%(beta_path))
        os.system('echo "feat %s/%s.fsf" >> jobs/bad_lss_job.txt'%(beta,trial))
