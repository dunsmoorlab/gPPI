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
    
    def phase_events(self,phase):
        ses_ = tasks[phase]['ses']
        return pd.read_csv(os.path.join(self.subj.subj_dir,'ses-'+str(ses_),'func',self.subj.fsub+'_ses-'+str(ses_)+'_task-'+phase+'_events.tsv'),sep='\t')

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

    def mem_events(self):
        outstr = {'CS+'      :'CSp',
                  'CS-'      :'CSm'}

        #walk through every folder containing the raw data and find all the event files
        for folder in os.walk(self.subj.subj_dir):
            for file in folder[2]:
                if 'events' in file and '.tsv' in file and 'memory' in file:
                    events = pd.read_csv(os.path.join(self.subj.subj_dir,folder[0],file), sep='\t')
                    phase = re.search('task-(.*)_events',file)[1]
                    out = os.path.join(self.subj.model_dir,'%s'%(phase))

                    #for every trial type
                    for con in events.trial_type.unique():
                        for encode in events.encode_phase.unique():
                            _timing = events[events.trial_type == con][events.encode_phase == encode][['onset','duration']]
                            _timing['PM'] = 1
                            _timing.to_csv( os.path.join(out, '%s_%s_all.txt'%(outstr[con],encode)),
                            sep='\t', float_format='%.8e', index=False, header=False)

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
                        #if 'cosine' in _c or 'motion_outlier' in _c:
                        if 'cosine' in _c:
                                run_COI.append(_c)
                    C = C[run_COI]
                    #C['constant'] = 1
                    C['framewise_displacement'][0] = 0 #really I think I should actually set this value to be the mean of the column, so that in FSL it doesn't impact how the FD is demeaned? but the value is truly 0 here...
                    
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
            template = os.path.join(gPPI_codebase,'feats','lss','template_lss_%s.fsf'%(phase))

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

                #with open('jobs/lss_betas/%s_%s_job.txt'%(self.subj.fsub,phase),'w') as txt_file:
                #    txt_file.write('feat %s'%(outfeat))
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
    
    def __init__(self,sub,mask=None,phases='encode_mem'):

        self.subj = bids_meta(sub)
        if '_mask' in mask:
            self.mask_name = mask[:-5]
        else:
            self.mask_name = mask
        #self.mask = self.load_mask(mask)
        #self.data = self.load_clean_data(phases=phases)
        #self.extract_timecourse()
        # self.interact()
        self._autofill_fsf()

    def load_mask(self,mask):
        
        if mask is not None:
            if 'mask' in mask or 'lh_' in mask or 'rh_' in mask: 
                fname = os.path.join(self.subj.masks,'%s.nii.gz'%(mask))
            else:
                fname = os.path.join(self.subj.masks,'%s_mask.nii.gz')
        
        if not os.path.exists(fname):
            fname = os.path.join(self.subj.masks,'%s.nii.gz'%(mask))
        
        return get_data(fname)
    
    def _apply_mask(self,mask=None,target=None):

        coor = np.where(mask == 1)
        values = target[coor]
        if values.ndim > 1:
            values = np.transpose(values) #swap axes to get sample X feature
        return values
    
    def load_clean_data(self,phases=None):
        if phases == 'all':
            phases = tasks
        elif phases == 'encode_mem':
            phases = ['baseline','acquisition','extinction','memory_run-01','memory_run-02','memory_run-03']
        elif phases == 'mem':
            phases = ['memory_run-01','memory_run-02','memory_run-03']
        elif type(phases) is str:
            phases = [phases]
        
        #load the data
        # data = {phase: get_data(os.path.join(self.subj.model_dir,phase,'%s_%s_gPPI.feat'%(self.subj.fsub,phase),'filtered_func_data.nii.gz')) for phase in phases}
        data = {phase: get_data(os.path.join(self.subj.func,'%s_ses-%s_task-%s_space-MNI152NLin2009cAsym_desc-preproc_denoised_bold.nii.gz'%(self.subj.fsub,tasks[phase]['ses'],phase))) for phase in phases}
        # data = {phase: get_data(os.path.join(self.subj.func,file)) for phase in phases for file in os.listdir(self.subj.func) if phase in file}

        # mask the data and take mean timeseries
        data = {phase: self._apply_mask(mask=self.mask,target=data[phase]).mean(axis=1) for phase in data}
        
        # DO NOT clean the data 
        # 6/3 note - we are now using the denoised data, so this is unnecssary 
        # data = {phase: clean(data[phase][:,np.newaxis], #need this here to be feature X sample after meaning
                        
        #                 confounds=pd.read_csv(os.path.join(self.subj.model_dir,phase,'confounds.txt'),
        #                                 sep='\t',header=None).values,
                
        #                 t_r=2,detrend=False,standardize='zscore')
        #                                                     for phase in data}

        return data

    #extract the givin timecourse for each run
    def extract_timecourse(self): 
        #deconvolve
        # self.neuronal = {phase: self._deconvolve(self.data[phase]) for phase in self.data}


        for phase in self.data:
                # df = pd.Series(self.data[phase])
                out = os.path.join(self.subj.model_dir,phase,self.mask_name)
                mkdir(out)
                # df.to_csv(os.path.join(out,'%s_bold_signal.txt'%(self.mask_name)),
                    # sep='\t', float_format='%.8e', index=False, header=False)
                np.savetxt(os.path.join(out,'%s_bold_signal.txt'%(self.mask_name)),self.data[phase])
                # np.savetxt(os.path.join(out,'%s_neuronal_signal.txt'%(self.mask_name)),self.neuronal[phase])

    def interact(self):

        for phase in self.neuronal:
            out = os.path.join(self.subj.model_dir,phase,self.mask_name)
            for file in os.listdir(os.path.join(self.subj.model_dir,phase)):
                if 'confounds' not in file and '.txt' in file:
                    long_events = np.zeros(tasks[phase]['n_tr']*2)
                    events = pd.read_csv(os.path.join(self.subj.model_dir,phase,file),sep='\t',header=None)
                    if self.subj.num in [105,106] and tasks[phase]['ses'] == 1:
                        events *= (2/2.23)
                    for i in events.index:
                        onset = int(events.loc[i,0])
                        duration = int(events.loc[i,1])
                        long_events[onset:onset+duration] = 1
                    long_events = long_events[::2]
                    assert long_events.shape[0] == tasks[phase]['n_tr']

                    ppi = (self.neuronal[phase] - self.neuronal[phase].mean()) * long_events
                    ppi = pd.Series(ppi)
                    ppi.to_csv(os.path.join(out,'%s'%(file[:-4]+'_ppi'+file[-4:])),
                        sep='\t', float_format='%.8e', index=False, header=False)

    def _autofill_fsf(self):
        for phase in [task for task in tasks if 'mem' in task or tasks[task]['ses'] ==1]:
            if 'memory' in phase:       
                template = os.path.join(gPPI_codebase,'feats','gPPI','mem_encode_gPPI.fsf')

            elif phase == 'acquisition':
                template = os.path.join(gPPI_codebase,'feats','gPPI','acq_encode_gPPI.fsf')

            elif phase in ['baseline','extinction']:
                template = os.path.join(gPPI_codebase,'feats','gPPI','day1_encode_gPPI.fsf')

            replacements = {'SUBID':self.subj.fsub,
                            'RUNID':phase,
                            'ROI':self.mask_name}            
        
            #need to handle the special cases where the TR is longer
            #the reason we always do this is because in the templates I only put in the wildcard on day 1
            #i.e. TR_length is not in the day 2 templates, so it doesnt change anything
            if self.subj.num in [105,106]:
                replacements['TR_length'] = '2.23'
            else:
                replacements['TR_length'] = '2'
    
            outfeat = os.path.join(self.subj.feat_dir,'%s_%s_%s_mem_encode_gPPI.fsf'%(self.subj.fsub,phase,self.mask_name))

            with open(template) as infile: 
                with open(outfeat, 'w') as outfile:
                    for line in infile:
                        for src, target in replacements.items():
                            line = line.replace(src, target)
                        outfile.write(line)

            #also go ahead and make the job script here
            os.system('echo "feat %s" >> jobs/gPPI/%s_%s_gPPI_job.txt'%(outfeat,self.mask_name,phase)) 
        

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

# for roi in ['rh_hpc','hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']:
#     autofill_fsf(group=False,template='gPPI/mem_encode_lvl2_gPPI',ses=2,name='mem_encode_lvl2_gPPI',roi=roi)

def autofill_fsf(group=False,template='',ses=None,name=None,roi=None,run=None,subjects=all_sub_args):
    if 'template' in template: outstr = re.search('template_(.*)',template)[1]
    elif roi is not None:
        outstr = roi+'_'+name
    else:
        outstr = name
    if run is not None:
        outstr = name
    if group:
        if roi is not None: replacements = {'ROI':roi}
        if ses == 'mem':
            for cope in range(1,15):
                replacements['COPEID'] = 'cope%s'%(cope)
                outfeat = os.path.join(SCRATCH,'group_gPPI',roi,'%s_%s.fsf'%(outstr,cope))

                with open(os.path.join(gPPI_codebase,'feats','%s.fsf'%(template))) as infile: 
                    with open(outfeat, 'w') as outfile:
                        for line in infile:
                            for src, target in replacements.items():
                                line = line.replace(src, target)
                            outfile.write(line)

                #also go ahead and make the job script here
                os.system('echo "feat %s" >> jobs/%s_job.txt'%(outfeat,outstr))
        if ses == 1:
            for phase in ['baseline','acquisition','extinction']:
                outfeat = os.path.join(SCRATCH,'group_gPPI',roi,'%s_%s_%s.fsf'%(roi,phase,name))
                replacements['RUNID'] = phase
                with open(os.path.join(gPPI_codebase,'feats','%s.fsf'%(template))) as infile:
                    with open(outfeat, 'w') as outfile:
                        for line in infile:
                            for src, target in replacements.items():
                                line = line.replace(src, target)
                            outfile.write(line)
                os.system('echo "feat %s" >> jobs/%s_job.txt'%(outfeat,name))

    else:
        for sub in subjects:
            subj = bids_meta(sub)
            replacements = {'SUBID':subj.fsub}
            if roi is not None: replacements['ROI'] = roi
            if run is not None: replacements['RUNID'] = run
            #need to handle the special cases where the TR is longer
            if ses == 1 and sub in [105,106]:
                replacements['TR_length'] = '2.23'
            else:
                replacements['TR_length'] = '2'

            outfeat = os.path.join(subj.feat_dir,'%s_%s.fsf'%(subj.fsub,outstr))

            with open(os.path.join(gPPI_codebase,'feats','%s.fsf'%(template))) as infile: 
                with open(outfeat, 'w') as outfile:
                    for line in infile:
                        for src, target in replacements.items():
                            line = line.replace(src, target)
                        outfile.write(line)

            #also go ahead and make the job script here
            os.system('echo "feat %s" >> jobs/%s_job.txt'%(outfeat,outstr))

def wrap_lss_jobs():
    for sub in all_sub_args:
        subj = bids_meta(sub)
        os.system('launch -N 1 -n 24 -J lss_%s -s jobs/lss_betas/%s_job.txt -m achennings@utexas.edu -p normal -r 6:00:00 -A LewPea_MRI_Analysis'%(sub,subj.fsub))

    for i, job in enumerate(os.listdir('jobs/lss_betas')):
        if '.txt' in job:
            os.system('launch -N 1 -n 12 -J lss_%s -s jobs/lss_betas/%s -m achennings@utexas.edu -p normal -r 09:00:00 -A LewPea_MRI_Analysis'%(i,job))

    for job in ['acquisition','extinction','baseline']:
        os.system('launch -N 1 -n 12 -J %s -s jobs/%s_rsa_job.txt -m achennings@utexas.edu -p normal -r 12:00:00')

    for job in [1,2,3]:
        # os.system('launch -N 1 -n 24 -J %s -s jobs/memory_run-0%s_rsa_job.txt -m achennings@utexas.edu -p normal -r 2:00:00 -A LewPea_MRI_Analysis'%(job,job))
        os.system(f'launch -N 1 -n 24 -J m{job} -s jobs/gPPI/memory-{job}_gPPI_job.txt -m achennings@utexas.edu -p normal -r 10:00:00 -A LewPea_MRI_Analysis')

    #submit a bunch of jobs at once
    for job in range(28):
        os.system('launch -N 1 -n 12 -J flss_%s -s jobs/final_lss_job_%s.txt -m achennings@utexas.edu -p normal -r 12:00:00 -A LewPea_MRI_Analysis'%(job,job))

    #splitting up a bunch of jobs into different job scripts
    for i in range(bad.shape[0]):
        os.system('echo "%s" >> jobs/final_lss_job_%s.txt'%(bad[0][i],int(np.floor(i/12))))

    for run in [1,2,3]:
        for roi in ['dACC','mOFC','rh_hpc','lh_hpc','lh_amyg','rh_amyg']:
            os.system('launch -N 1 -n 24 -J %s_%s -s jobs/%s_memory_run-0%s_gPPI_job.txt -m achennings@utexas.edu -p normal -r 3:00:00 -A LewPea_MRI_Analysis'%(run,roi,roi,run))

    for run in [1,2]:
        os.system(f'launch -N 2 -n 24 -J localizer-{run} -s jobs/localizer-{run}_job.txt -m achennings@utexas.edu -p normal -r 1:00:00 -A LewPea_MRI_Analysis')
        os.system('launch -N 2 -n 24 -J localizer_lvl2 -s jobs/localizer_lvl2_job.txt -m achennings@utexas.edu -p normal -r 1:00:00 -A LewPea_MRI_Analysis')

    # for roi in ['rh_hpc','hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']:
    for roi in ['rh_hpc']:
        # os.system('launch -N 1 -n 24 -J %s_lvl2 -s jobs/%s_mem_encode_lvl2_gPPI_job.txt -m achennings@utexas.edu -p normal -r 00:45:00 -A LewPea_MRI_Analysis'%(roi,roi))        
        os.system('launch -N 1 -n 14 -J %s_lvl3 -s jobs/%s_group_cope_job.txt -m achennings@utexas.edu -p normal -r 2:00:00 -A LewPea_MRI_Analysis -d 2897228'%(roi,roi))
    os.system('launch -N 1 -n 18 -J day1_lvl3 -s jobs/group_day1_gPPI_job.txt -m achennings@utexas.edu -p normal -r 4:00:00 -A LewPea_MRI_Analysis'%(roi,roi))

    for phase in ['baseline','acquisition','extinction']:
        os.system('launch -N 1 -n 6 -J %s_lvl3 -s jobs/%s_group_gPPI_job.txt -m achennings@utexas.edu -p normal -r 3:00:00 -A LewPea_MRI_Analysis'%(phase,phase))

    #launch -N 2 -n 34 -J sm -s jobs/source_memory_job.txt -m achennings@utexas.edu -p normal -r 3:00:00 -A LewPea_MRI_Analysis
    #launch -N 2 -n 34 -J sm_reg -s jobs/source_memory_reg_job.txt -m achennings@utexas.edu -p normal -r 0:30:00 -A LewPea_MRI_Analysis -d 3051447
    #launch -N 3 -n 30 -J sm_lvl3 -s jobs/source_memory_lvl3_job.txt -m achennings@utexas.edu -p normal -r 3:00:00 -A LewPea_MRI_Analysis -d 3051447
    #launch -N 1 -n 64 -p largemem512GB -J 3dlme -s jobs/3dlme_job.txt -m achennings@utexas.edu -r 48:00:00 -A LewPea_MRI_Analysis
    #launch -N 48 -n 48 -J smooth -s jobs/afni_fwhm_job.txt -m achennings@utexas.edu -p normal -r 2:00:00 -A LewPea_MRI_Analysis
    #launch -N 1 -n 64 -p largemem512GB -J 3dmvm -s jobs/3dMVM_job.txt -m achennings@utexas.edu -r 48:00:00 -A LewPea_MRI_Analysis
    os.system(f'launch -N 2 -n 48 -J univariate -s jobs/memory_univariate_lvl2_job.txt -m achennings@utexas.edu -p normal -r 2:00:00 -A LewPea_MRI_Analysis')
    # for phase in ['baseline','acquisition','extinction','memory_run-01','memory_run-02','memory_run-03']:
    for phase in ['memory_run-01','memory_run-02','memory_run-03']:
        # for roi in ['rh_hpc']:
        # for roi in ['rACC','sgACC','rh_hpc','lh_hpc','lh_amyg','rh_amyg']:
        for roi in ['hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']:
        # for roi in ['rh_hc_tail','lh_hc_tail','rh_hc_body','lh_hc_body','rh_hc_head','lh_hc_head','rh_amyg_bla','lh_amyg_bla','rh_amyg_cem','lh_amyg_cem']:
            # os.system('launch -N 1 -n 24 -J %s_gPPI -s jobs/%s_gPPI_job.txt -m achennings@utexas.edu -p normal -r 2:00:00 -A LewPea_MRI_Analysis'%(phase,phase))
            os.system(f'launch -N 1 -n 24 -J {roi}_{phase} -s jobs/gPPI/{roi}_{phase}_gPPI_job.txt -m achennings@utexas.edu -p normal -r 01:30:00 -A LewPea_MRI_Analysis')

def clean_bad_lss():
    bad = pd.read_csv('bad_lss.txt',header=None)
    for beta in bad[0]:
        trial = beta[-8:]
        for beta_path in [beta+'.feat', beta+'/'+trial+'.feat', beta+'+.feat', beta+'/'+trial+'+.feat']:
            if os.path.exists(beta_path):
                # print(os.path.exists(beta_path))
                os.system('rm -r %s'%(beta_path))
        os.system('echo "feat %s/%s.fsf" >> jobs/bad_lss_job.txt'%(beta,trial))

def wrangle_first_level_rsa():
    mem_weights = { 1:'mem_baseline_CSp.nii.gz',
                    2:'mem_baseline_CSm.nii.gz',
                    3:'mem_acquisition_CSp.nii.gz',
                    4:'mem_acquisition_CSm.nii.gz',
                    5:'mem_extinction_CSp.nii.gz',
                    6:'mem_extinction_CSm.nii.gz',
                    7:'mem_foil_CSp.nii.gz',
                    8:'mem_foil_CSm.nii.gz'}

    for sub in all_sub_args:
        subj = bids_meta(sub)
        print(sub)
        # for phase in ['baseline','acquisition','extinction']:
        #     csp_cope = os.path.join(subj.model_dir,phase,'rsa_model.feat','stats','cope1.nii.gz')
        #     csm_cope = os.path.join(subj.model_dir,phase,'rsa_model.feat','stats','cope2.nii.gz')

        #     os.system('cp %s %s'%(csp_cope,os.path.join(subj.weights,'%s_CSp.nii.gz'%(phase))))
        #     os.system('cp %s %s'%(csm_cope,os.path.join(subj.weights,'%s_CSm.nii.gz'%(phase))))

        for cope in mem_weights:
            in_file = os.path.join(subj.model_dir,'feats','memory_rsa_lvl2.gfeat',f'cope{cope}.feat','stats','cope1.nii.gz')
            out_file = os.path.join(subj.weights,mem_weights[cope])

            os.system(f'cp {in_file} {out_file}')

def prep_lvl2():
    import os
    from fg_config import bids_meta, all_sub_args
    
    for sub in all_sub_args:
        print(sub)
        subj = bids_meta(sub)
        # for phase in ['memory_run-01','memory_run-02','memory_run-03']:
        for phase in ['baseline','acquisition','extinction']:
            regfolder = os.path.join(subj.model_dir,phase,'%s_%s_reg_gPPI.feat'%(subj.fsub,phase),'reg')

            for roi in ['rh_hpc','hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']:
                dest = os.path.join(subj.model_dir,phase,roi,'source.feat','reg')

                os.system('cp -r %s %s'%(regfolder,dest))
                

def motion_outlier_count():
    import os
    import re
    import matplotlib.pyplot as plt
    import seaborn as sns

    '''to build the output you need a list of all subjects "all_sub_args"
    and a list of all tasks "tasks"
    '''
    df = pd.DataFrame(index=pd.MultiIndex.from_product([all_sub_args,tasks]))
    df['censored'] = 0.
    for sub in all_sub_args:
        '''
        to make this work for CCX, comment out the next line
        and change the subj.prep_dir so it is the path to each subject's
        fMRIprep directory
        '''
        subj = bids_meta(sub)
        for folder in os.walk(subj.prep_dir): 
            for file in folder[2]:
                if 'confounds' in file and '.tsv' in file:
                    C = pd.read_csv(os.path.join(subj.prep_dir,folder[0],file), sep='\t')
                    mo = [c for c in C.columns if 'motion_outlier' in c]
                    task = re.search('task-(.*)_desc',file)[1]
                    df.loc[(sub,task),'censored'] = len(mo)/C.shape[0]

    sns.boxplot(data=df,y='censored')
    sns.swarmplot(data=df,y='censored',color='black')

from fg_config import *
def copypasta(subs=all_sub_args,feat=None):
    os.system('rm copypasta.txt')
    for sub in subs:
        subj = bids_meta(sub)
        os.system('echo "%s" >> copypasta.txt'%(os.path.join(subj.feat_dir,subj.fsub+'_'+feat+'.feat')))
