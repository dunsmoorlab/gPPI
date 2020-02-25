import os
import re

import pandas as pd
import nibabel as nib

from collections import OrderedDict
from nilearn.image import concat_imgs

from fg_config import *

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
    
    def lss_betas(self):
        for folder in os.walk(self.subj.subj_dir):
            for file in folder[2]:
                if 'events' in file and '.tsv' in file:
                    events = pd.read_csv(os.path.join(self.subj.subj_dir,folder[0],file), sep='\t')
                    phase = re.search('task-(.*)_events',file)[1]
                    lss_dir = os.path.join(self.subj.model_dir,'%s'%(phase),'lss_betas')
                    mkdir(lss_dir)

                    trial_types = events.trial_type.unique()

                    # for trial in range(events.shape[0]):
                    #     beta_folder = os.path.join(lss_dir,'trial_{0:0=2d}'.format(trial))
                    #     mkdir(beta_folder)
                    #     beta_trial = events.loc[trial,['onset','duration']]
                    #     beta_trial['PM'] = 1
                    #     beta_trial.to_csv(os.path.join(beta_folder,'beta.txt'),
                    #                 sep='\t', float_format='%.8e', index=False, header=False)

                    #     for i, condition in enumerate(trial_types):
                    #         #grab all trials of each condition
                    #         con_trials = np.where(events.trial_type == condition)[0]
                    #         #make sure we don't model the same trial twice
                    #         con_trials = [t for t in con_trials if t != trial]
                    #         #get the onsets/durations
                    #         con_events = events.loc[con_trials,['onset','duration']]
                    #         con_events['PM'] = 1
                    #         con_events.to_csv(os.path.join(beta_folder,'no_interest_%s.txt'%(i)),
                    #                         sep='\t', float_format='%.8e', index=False, header=False)

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

                # with open(template) as infile: 
                #     with open(outfeat, 'w') as outfile:
                #         for line in infile:
                #             for src, target in replacements.items():
                #                 line = line.replace(src, target)
                #             outfile.write(line)

                # with open('jobs/lss_betas/%s_%s_job.txt'%(self.subj.fsub,phase),'w') as txt_file:
                #     txt_file.write('feat %s'%(outfeat))
                #also go ahead and make the job script here
                os.system('echo "feat %s" >> jobs/lss_betas/%s_job.txt'%(outfeat,self.subj.fsub))

    def lss_reconstruct(self):
        for task in tasks:
            lss_dir = os.path.join(self.subj.model_dir,task,'lss_betas')
            beta_fname = os.path.join(self.subj.beta,'%s_beta.nii.gz'%(task))
            
            beta_img = OrderedDict()
            concat_ = True
            for i in range(tasks[task]['n_trials']):
                trial = 'trial_{0:0=2d}'.format(i)
                try:
                    beta_img[i] = nib.load(os.path.join(lss_dir,trial,'%s.feat'%(trial),'stats','cope1.nii.gz'))
                except:
                    print(self.subj.num,task,trial)
                    concat_ = False
    
            if concat_:
                # concatenate them
                beta_img = concat_imgs(beta_img.values())
                nib.save(beta_img,beta_fname)

            #mask them too
            os.system('fslmaths %s -mas %s %s'%(beta_fname,self.subj.refvol_mask,beta_fname))


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
        os.system('launch -N 1 -n 12 -J %s -s jobs/%s_rsa_job.txt -m achennings@utexas.edu -p normal -r ')
# q = [os.path.join(self.subj.prep_dir,folder[0],file) )
#                         for folder in os.walk(self.subj.prep_dir)
#                         for file in folder[2]
#                         if 'T1w_desc-brain_mask.nii.gz' in file]
