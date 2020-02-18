import os
import re

import pandas as pd

from fg_config import *

class bids_events():
    
    def __init__(self,sub):

        self.subj = bids_meta(sub)
        self.fsl_events()
        self.confounds()
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

        outfeat = os.path.join(subj.feat_dir,'%s_%s.feat'%(subj.fsub,outstr))

        with open(os.path.join(gPPI,'feats','%s.fsf'%(template))) as infile: 
            with open(outfeat, 'w') as outfile:
                for line in infile:
                    for src, target in replacements.items():
                        line = line.replace(src, target)
                    outfile.write(line)

        #also go ahead and make the job script here
        os.system('echo "feat %s" >> %s_job.txt'%(outfeat,outstr))


