import os
import re

import pandas as pd

from fc_config import *

class bids_events():
    
    def __init__(self,sub):

        self.subj = bids_meta(sub)

    def CS_fsl(self): #generate CS+, CS-, US timing files for use in FSL

        outstr = {'CS+':'CSp',
                  'CS-':'CSm'}

        for folder in os.walk(subj.subj_dir):
            for file in folder[2]:
                if 'events' in file and '.tsv' in file:
                    events = pd.read_csv(os.path.join(subj.subj_dir,folder[0],file), sep='\t')
                    phase = re.search('task-(.*)_events',file)[1]
                    out = os.path.join(subj.model_dir,'task-%s'%(phase))
                    mkdir(out)

                    if 'localizer' not in phase:
                        for con in events.trial_type.unique():
                            con_timing = events[events.trial_type == con][['onset','duration']]
                            con_timing['PM'] = 1
                            con_timing.to_csv( os.path.join(out, '%s_all.txt'%(outstr[con])),
                                sep='\t', float_format='%.8e', index=False, header=False)

                    #handle the localizer case
                    #handle early/late for associative learning phases
                    


