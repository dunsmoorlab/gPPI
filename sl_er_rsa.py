import os
import pickle
import numpy as np
import pandas as pd
import nibabel as nib

from nilearn.decoding import SearchLight
from sklearn.model_selection import GroupKFold
from nilearn.image import index_img, new_img_like, get_data
from sklearn.base import BaseEstimator, ClassifierMixin
from collections import OrderedDict

from fg_config import *
from bids_model import bids_events

class RSA_estimator(BaseEstimator, ClassifierMixin): #ClassifierMixin most likely interchangeable
    def __init__(self,subject=None):
        
        self.subject = subject #sklearn requires some arg

    def fit(self,X,y=None,*args,**kwargs):

        self.X_train = X #grab half the data for "training"
        self.y_train = y #print this for validation check
        return self

    def score(self,X,y=None):
                
        r = np.corrcoef(self.X_train,X)[0,1] #Compute RSA here
        z = np.arctanh(r) #Convert to fischer-z
        return z    

class ER_searchlight():
    def __init__(self,sub,process_mask=None,run=False,reg=False):
        self.verbose=False

        self.subj = bids_meta(sub)
        self.wb_mask = nib.load(std_2009_brain_mask_3mm)
        self.process_mask = nib.load(std_2009_brain_mask_3mm)
        self.refvol = nib.load(std_2009_brain_3mm)
    
        #some stuff to handle the beta weights
        weights = 'copes' 

        self.conditions = {'CS+': 'CSp',
                           'CS-': 'CSm'}
        
        self.w_slices = {'CS+':{
                             'baseline':slice(0,24),
                          'acquisition':slice(24,48),                     
                           'extinction':slice(48,72)},
                         'CS-':{
                             'baseline':slice(72,96),
                          'acquisition':slice(96,120),                     
                           'extinction':slice(120,144)}
                           }

        # self.w_phases = {'baseline':'baseline','fear_conditioning':'fear','extinction':'ext'}
    
        self.encoding_phases = ['baseline','acquisition','extinction']
        self.mem_phases = ['memory_run-01','memory_run-02','memory_run-03']
        self.phases = self.encoding_phases + self.mem_phases
        
        self.load_data()
        self.init_SL()
        if run: self.wrap_run_SL()
        if reg: self.reg_res()
    
    def load_data(self):
        
        self.encoding_data = OrderedDict()
        self.encoding_labels = OrderedDict()
        
        self.mem_data = OrderedDict()
        self.mem_labels = OrderedDict()

        self.W = OrderedDict()

        for phase in self.phases:
            
            beta_img = get_data(os.path.join(self.subj.beta, phase+'_beta_std.nii.gz'))
            
            if phase in self.encoding_phases:
                events = bids_events(self.subj.num).phase_events(phase=phase)
                for cs in self.conditions:
                    _con = self.conditions[cs]+'_'+'trial'
                    events[_con] = ''
                    events.loc[np.where(events.trial_type==cs)[0],_con] = [i for i in range(1,25)]

                events['phase'] = phase
                self.encoding_labels[phase] = events
                self.encoding_data[phase] = beta_img
        
                #load in the weights here
                self.W[phase] = OrderedDict()
                for con in self.conditions:
                    self.W[phase][con] = get_data(
                        os.path.join(
                        self.subj.weights,'%s_%s_std.nii.gz'%
                        (phase,self.conditions[con])))

            elif phase in self.mem_phases:
                events = bids_events(self.subj.num).phase_events(phase=phase)
                # foil_mask = events.memcond.isin(['Old'])
                events['phase'] = phase
                self.mem_labels[phase] = events          
                self.mem_data[phase] = beta_img

        self.encoding_labels = pd.concat(self.encoding_labels.values(),sort=False)
        self.encoding_labels.reset_index(inplace=True,drop=True) #need this bc phase events are all numbered the same
        # self.encoding_data = nib.concat_images(self.encoding_data.values(),axis=3)
        self.encoding_data = np.concatenate([self.encoding_data['baseline'],self.encoding_data['acquisition'],self.encoding_data['extinction']],axis=-1)

        self.mem_labels = pd.concat(self.mem_labels.values(),sort=False)
        foil_mask = self.mem_labels.memory_condition.isin(['Old'])
        self.mem_labels = self.mem_labels[foil_mask].reset_index(drop=True)
        # self.mem_data = image.index_img(nib.concat_images(self.mem_data.values(),axis=3), foil_mask)
        self.mem_data = np.concatenate([self.mem_data['memory_run-01'],self.mem_data['memory_run-02'],self.mem_data['memory_run-03']],axis=-1)
        foil_mask = foil_mask.values.ravel()
        self.mem_data = self.mem_data[:,:,:,foil_mask]

    def init_SL(self,process_mask=None):
        

        # if process_mask is not None:
        #     process_mask = nib.load(os.path.join(self.subj.roi,process_mask+'_mask.nii.gz'))
        
        self.sl = SearchLight(
                    mask_img=self.wb_mask, #whole brain mask
                    process_mask_img=self.process_mask, #where to run the searchlight (None = all)
                    radius=6, #radius in mm
                    estimator=RSA_estimator(self.subj.num), #just pass something for sklearn
                    scoring=None, #None defaults to our custom score() function in RSA_estimator
                    cv=GroupKFold(n_splits=2), #not used for RSA
                    verbose=0,
                    n_jobs=-1, #-1: all processors
                    )
    
    def run_SL(self,ER_item,labels,groups):

        self.sl.fit(ER_item,labels,groups)

        return self.sl.scores_

    def wrap_run_SL(self):
        labels = ['encode','mem']
        groups = [1,2]

        self.rsa = self.mem_labels.copy().drop(columns=['onset','duration','response_time'])
        for cs in self.conditions: self.rsa[self.conditions[cs]+'_'+'trial'] = ''

        self.ER_res = OrderedDict()
        
        for i, stim in enumerate(self.rsa.stimulus):
            print(i)
            mem_loc = i
            encoding_loc = np.where(self.encoding_labels.stimulus == stim)[0][0]
            _phase = self.rsa.loc[mem_loc,'encode_phase']
            _trial_type = self.rsa.loc[mem_loc,'trial_type']
            _con        = self.conditions[_trial_type]+'_'+'trial'
            self.rsa.loc[i,_con] = self.encoding_labels.loc[encoding_loc,_con] 
            
            encoding_trial = self.encoding_data[:,:,:,encoding_loc] * self.W[_phase][_trial_type]
            mem_trial      = self.mem_data[:,:,:,mem_loc] * self.W[_phase][_trial_type]

            ER_item = np.stack([encoding_trial,mem_trial],axis=-1)
            
            ER_item = new_img_like(self.refvol,ER_item,copy_header=True)

            self.ER_res[i] = self.run_SL(ER_item,labels,groups)
        
        self.ER_res = np.array([*self.ER_res.values()])
        self.ER_res = np.moveaxis(self.ER_res,0,-1)

        with open(os.path.join(self.subj.rsa,'sl_er.p'),'wb') as file:
            pickle.dump(self.ER_res,file)
    
    def reg_res(self):     
        self.save_str = os.path.join(self.subj.rsa,'item_ER.nii.gz')        

        ref2std = os.path.join(self.subj.subj_dir,'ref2std.mat')

        std = os.path.join(WORK,'standard','MNI152_T1_3mm_brain.nii.gz')

        os.system('flirt -in %s -ref %s -out %s -init %s -applyxfm'%(self.save_str,std,self.save_str,ref2std))
