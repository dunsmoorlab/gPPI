import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pingouin as pg
import nibabel as nib

from matplotlib.ticker import MultipleLocator, ScalarFormatter
from collections import OrderedDict
from nilearn import image
from nilearn.image import get_data
from nilearn.input_data import NiftiMasker
from scipy.stats import pearsonr, spearmanr, zscore
from sklearn.linear_model import LogisticRegression
from scipy.special import expit

from fg_config import *
from bids_model import bids_events
# from glm_timing import glm_timing

class roi_rsa():
    def __init__(self,sub=None,fs=True,hemi=False):

        self.subj = bids_meta(sub)
        self.verbose=False
        self.fs = fs
        self.hemi = hemi
        #set the phases
        self.encoding_phases = ['baseline','acquisition','extinction']
        self.mem_phases = ['memory_run-01','memory_run-02','memory_run-03']
        self.phases = self.encoding_phases + self.mem_phases
        
        #hardcode copes as weights (as opposed to zmaps)
        weights = 'copes'
        # self.weight_dir = os.path.join(data_dir,'rsa_%s'%(weights),self.subj.fsub)
        #need these to load the weights becauase I named them lazily
        self.conditions = {'CS+': 'CSp',
                           'CS-': 'CSm'}
        # self.w_phases = {'baseline'         : 'baseline',
        #                  'fear_conditioning': 'fear',
        #                  'extinction'       : 'ext'
        self.w_slices = {'CS+':{
                             'baseline':slice(0,24),
                          'acquisition':slice(24,48),                     
                           'extinction':slice(48,72)},
                         'CS-':{
                             'baseline':slice(72,96),
                          'acquisition':slice(96,120),                     
                           'extinction':slice(120,144)}
                           }

        
        #hardcode rois for now
        #if self.fs: self.rois = ['mOFC','dACC','amyg_cem','amyg_bla','hc_head','hc_body','hc_tail'] 
        # if self.fs: self.rois = ['mOFC','dACC','amyg','hpc','ins','hc_head','hc_body','hc_tail','rh_hc_head','rh_hc_body','rh_hc_tail','lh_hc_head','lh_hc_body','lh_hc_tail','amyg_bla','amyg_cem','rh_amyg_bla','rh_amyg_cem','lh_amyg_bla','lh_amyg_cem'] 
        if self.fs: self.rois = ['mOFC','dACC','amyg','hpc','ins','lh_amyg','rh_amyg','lh_hpc','rh_hpc','sgACC','rACC','rSMA','rACG','hc_head','hc_body','hc_tail','rh_hc_head','rh_hc_body','rh_hc_tail','lh_hc_head','lh_hc_body','lh_hc_tail','amyg_bla','amyg_cem','rh_amyg_bla','rh_amyg_cem','lh_amyg_bla','lh_amyg_cem']  
            # if hemi:
                # self.rois = ['rh_hc_head','rh_hc_body','rh_hc_tail','rh_amyg_bla','rh_amyg_cem',
                #              'lh_hc_head','lh_hc_body','lh_hc_tail','lh_amyg_bla','lh_amyg_cem']
            # else:
            #     self.rois = ['amyg_cem','amyg_bla','hc_head','hc_body','hc_tail']

        # else: self.rois = ['hippocampus','mOFC','dACC','insula','amygdala']
        
        #data needs to be loaded WITHOUT mask to facilitate more intricate analyses
        self.load_data() 
        # self.compute_item_rsa()
        # self.compute_cross_rsa()
        self.compute_mem_mats()

    def load_data(self):
            
        #encoding data & labels
        self.encoding_data = OrderedDict()
        self.encoding_labels = OrderedDict()

        #ls-u style beta weights
        self.W = OrderedDict()
        
        #retrieva data & labels
        self.mem_data = OrderedDict()
        self.mem_labels = OrderedDict()

        for phase in self.phases:

            beta_img = get_data(os.path.join(self.subj.beta, phase+'_beta.nii.gz'))

            if phase in self.encoding_phases:
                events = bids_events(self.subj.num).phase_events(phase=phase)
                for cs in self.conditions:
                    _con = self.conditions[cs]+'_'+'trial'
                    events[_con] = ''
                    events.loc[np.where(events.trial_type==cs)[0],_con] = [i for i in range(1,25)]
                    

                events['phase'] = phase #this isn't in the dataframe yet
                #save labels & data
                self.encoding_labels[phase] = events
                self.encoding_data[phase] = beta_img
                
                #load in the weights here
                self.W[phase] = OrderedDict()
                for con in self.conditions:
                    self.W[phase][con] = get_data(
                            os.path.join(
                            self.subj.weights,'%s_%s.nii.gz'%
                            (phase,self.conditions[con])))
        
            elif phase in self.mem_phases:
                events = bids_events(self.subj.num).phase_events(phase=phase)
                events['phase'] = phase #add phase to df
                #save labels & data
                self.mem_labels[phase] = events
                self.mem_data[phase] = beta_img
        
        #get 1 dataframe & image for encoding
        self.encoding_labels = pd.concat(self.encoding_labels.values(),sort=False)
        self.encoding_labels.reset_index(inplace=True,drop=True) #need this bc phase events are all numbered the same
        self.encoding_data = np.concatenate([self.encoding_data['baseline'],self.encoding_data['acquisition'],self.encoding_data['extinction']],axis=-1)
        print('ENCODING DATA SHAPE =',self.encoding_data.shape)
        #same for retrieval, except we have to remove foils
        self.mem_labels = pd.concat(self.mem_labels.values(),sort=False)
        self.all_mem_labels = self.mem_labels.copy().reset_index(drop=True)

        foil_mask = self.mem_labels.memory_condition.isin(['Old'])

        self.mem_labels = self.mem_labels[foil_mask].reset_index(drop=True)
        
        self.all_mem_data = np.concatenate([self.mem_data['memory_run-01'],self.mem_data['memory_run-02'],self.mem_data['memory_run-03']],axis=-1)
        self.mem_data = np.concatenate([self.mem_data['memory_run-01'],self.mem_data['memory_run-02'],self.mem_data['memory_run-03']],axis=-1)

        print('MEM_DATA SHAPE =',self.mem_data.shape)
        #self.mem_data = np.array([self.mem_data[:,:,:,i] for i in foil_mask if i is True])
        foil_mask = foil_mask.values.ravel()
        self.mem_data = self.mem_data[:,:,:,foil_mask]
        print('MEM_DATA SHAPE =',self.mem_data.shape)


        #we need this for the roi bootstrapping, probably best to do it here and broadcast
        self.dACC_nvox = np.where( get_data(os.path.join(self.subj.masks,'dACC_mask.nii.gz')) == 1 )[0].shape[0]

    def apply_mask(self,roi=None,target=None):
        #pass in an roi and a nifti image to return the data located in the roi mask
        try:
            mask = get_data(os.path.join(self.subj.masks,'%s.nii.gz'%(roi)))
        except:
            mask = get_data(os.path.join(self.subj.masks,'%s_mask.nii.gz'%(roi)))
        # else: mask = nib.load(os.path.join(self.subj.roi,'%s_mask.nii.gz'%(roi)))
        coor = np.where(mask == 1)
        values = target[coor]
        if values.ndim > 1:
            values = np.transpose(values) #swap axes to get feature X sample
        return values

    def boot_rsa(self,encoding_trial,mem_trial):
        assert encoding_trial.shape == mem_trial.shape
        n_boot = 1000
        boot_res = np.zeros(n_boot)
        for i in range(n_boot):
            _samp       = np.random.randint(low=0, high=encoding_trial.shape[0], size=self.dACC_nvox)
            boot_res[i] = np.arctanh( pearsonr( encoding_trial[_samp], mem_trial[_samp] )[0] )
        return boot_res.mean()

    def compute_item_rsa(self):
    
        self.rsa = self.mem_labels.copy().drop(columns=['onset','duration'])
        for cs in self.conditions: self.rsa[self.conditions[cs]+'_'+'trial'] = ''
        for roi in self.rois:
            print(roi)
            self.rsa[roi] = 0 #going by columns and then melting at the end
            
            #apply the mask to everything for this roi
            encoding_data = self.apply_mask(roi,self.encoding_data)
            mem_data      = self.apply_mask(roi,self.mem_data)
            W = {}
            for phase in self.encoding_phases:
                W[phase] = {}
                for con in self.conditions:
                    W[phase][con] = self.apply_mask(roi,self.W[phase][con])
            
            #great, you've got everything in the mask shape, now run the item rsa
            for i, stim in enumerate(self.rsa.stimulus):
                mem_loc = i
                encoding_loc = np.where(self.encoding_labels.stimulus == stim)[0][0]
                
                _phase = self.rsa.loc[mem_loc,'encode_phase']
                _trial_type = self.rsa.loc[mem_loc,'trial_type']
                _con        = self.conditions[_trial_type]+'_'+'trial'
                self.rsa.loc[i,_con] = self.encoding_labels.loc[encoding_loc,_con] 
                
                encoding_trial = encoding_data[encoding_loc] * W[_phase][_trial_type]
                mem_trial      = mem_data[mem_loc] * W[_phase][_trial_type]
                
                #if roi == 'mOFC':
                #    self.rsa.loc[i,roi] = self.boot_rsa(encoding_trial,mem_trial)
                #else:
                p = pearsonr(encoding_trial,mem_trial)[0]
                z = np.arctanh(p)
                self.rsa.loc[i,roi] = z
                    #np.corrcoef()
                

        self.rsa['subject'] = self.subj.num
        self.rsa = self.rsa.melt(id_vars=['subject','trial_type','stimulus','memory_condition','encode_phase','response',
                               'low_confidence_accuracy','high_confidence_accuracy','phase','CSp_trial','CSm_trial'],
                      value_vars=self.rois
                  ).rename(columns={'variable':'roi', 'value':'rsa'})

        if self.hemi:
                self.rsa.to_csv(os.path.join(self.subj.rsa,'roi_ER_HCA_hemi.csv'),index=False)
        
        elif self.fs:
            self.rsa.to_csv(os.path.join(self.subj.rsa,'fs_mask_roi_ER.csv'), index=False)

    def compute_cross_rsa(self):
        self.cross_mats = {}
        self.encoding_labels.phase = pd.Categorical(self.encoding_labels.phase,self.encoding_phases,ordered=True)
        encode_reorder = list(self.encoding_labels.sort_values(by=['trial_type','phase']).index)
        mem_reorder = [np.where(self.mem_labels.stimulus == self.encoding_labels.loc[i,'stimulus'])[0][0] for i in encode_reorder]
        for roi in self.rois:
            self.cross_mats[roi] = {}
            print(roi)
            
            #apply the mask to everything for this roi
            encoding_data = self.apply_mask(roi,self.encoding_data)[encode_reorder]
            mem_data      = self.apply_mask(roi,self.mem_data)[mem_reorder]
            W = {}
            for phase in self.encoding_phases:
                W[phase] = {}
                for con in self.conditions:
                    W[phase][con] = self.apply_mask(roi,self.W[phase][con])
                    encoding_data[self.w_slices[con][phase]] = encoding_data[self.w_slices[con][phase]] * W[phase][con]
                    mem_data[self.w_slices[con][phase]] = mem_data[self.w_slices[con][phase]] * W[phase][con]
       

            encode_mat = np.arctanh(np.corrcoef(encoding_data))
            encode_mat[np.eye(encode_mat.shape[0],dtype=bool)] = 0
       
            mem_mat = np.arctanh(np.corrcoef(mem_data))
            mem_mat[np.eye(mem_mat.shape[0],dtype=bool)] = 0

            ers_mat = np.arctanh(np.corrcoef(encoding_data,mem_data))
            ers_mat[np.eye(ers_mat.shape[0],dtype=bool)] = 0
            self.cross_mats[roi]['encoding'] = encode_mat
            self.cross_mats[roi]['retrieval'] = mem_mat
            self.cross_mats[roi]['ers'] = ers_mat
        with open(os.path.join(self.subj.rsa,'cross_mats.p'),'wb') as file:
            pickle.dump(self.cross_mats,file)

    def compute_mem_mats(self):
        self.mem_mats = {}
        mem_phase4 = ['baseline','acquisition','extinction','foil']
        self.all_mem_labels.encode_phase = pd.Categorical(self.all_mem_labels.encode_phase,mem_phase4,ordered=True)
        mem_reorder = list(self.all_mem_labels.sort_values(by=['trial_type','encode_phase']).index)
        for roi in self.rois:
            print(roi)

            mem_data = self.apply_mask(roi,self.all_mem_data)[mem_reorder]

            mem_mat = np.arctanh(np.corrcoef(mem_data))
            mem_mat[np.eye(mem_mat.shape[0],dtype=bool)] = 0
            self.mem_mats[roi] = mem_mat

        with open(os.path.join(self.subj.rsa,'mem_mats.p'),'wb') as file:
                    pickle.dump(self.mem_mats,file)

class group_roi_rsa():

    def __init__(self,group='control',ext_split=True,fs=True,hemi=False):
        
        self.group = group
        self.ext_split = ext_split
        self.fs = fs
        self.hemi = hemi

        if   group == 'control': self.subs = sub_args;
        elif group == 'ptsd':    self.subs = p_sub_args
        elif group == 'between': self.subs = all_sub_args

        

        if self.hemi:
            self.rois = ['rh_amyg_cem','lh_amyg_cem',
                      'rh_amyg_bla','lh_amyg_bla',
                      'rh_hc_head' ,'lh_hc_head',
                      'rh_hc_body' ,'lh_hc_body',
                      'rh_hc_tail' ,'lh_hc_tail']
        elif self.fs:
                # self.rois = ['vmPFC','dACC','amyg_cem','amyg_bla','hc_head','hc_body','hc_tail']
                # self.rois = ['fvmPFC','fdACC','mOFC','dACC','amyg','hpc','ins']
                #self.rois = ['mOFC','dACC','amyg','hpc','ins','hc_head','hc_body','hc_tail','rh_hc_head','rh_hc_body','rh_hc_tail','lh_hc_head','lh_hc_body','lh_hc_tail','amyg_bla','amyg_cem','rh_amyg_bla','rh_amyg_cem','lh_amyg_bla','lh_amyg_cem']
                self.rois = ['mOFC','dACC','amyg','hpc','ins','lh_amyg','rh_amyg','lh_hpc','rh_hpc','sgACC','rACC','rSMA','rACG','hc_head','hc_body','hc_tail','rh_hc_head','rh_hc_body','rh_hc_tail','lh_hc_head','lh_hc_body','lh_hc_tail','amyg_bla','amyg_cem','rh_amyg_bla','rh_amyg_cem','lh_amyg_bla','lh_amyg_cem']    


        self.conditions = {'CS+': 'CSp',
                           'CS-': 'CSm'}
        
        if self.ext_split: self.encoding_phases = ['baseline','acquisition','early_extinction','extinction']
        else:              self.encoding_phases = ['baseline','acquisition','extinction']

        self.load_rois()
        self.load_cross_mats()
        self.load_mem_mats()
        
    def load_rois(self):

        df = {} #ultimate output

        for sub in self.subs:
            subj = bids_meta(sub)
            # self.hemi: df[sub] = pd.read_csv(os.path.join(subj.rsa,'roi_ER_HCA_hemi.csv'))
            df[sub] = pd.read_csv(os.path.join('rsa_results',subj.fsub,'fs_mask_roi_ER.csv'))

        self.df = pd.concat(df.values()).reset_index(drop=True)
        #lets label the first 8 trials of extinction as "early_extinction"
        if self.ext_split:
            for cs in self.conditions:
                con = '%s_trial'%(self.conditions[cs])
                for i in range(1,9): 
                    self.df.loc[ self.df[ self.df.encode_phase == 'extinction' ][ self.df[con] == i ].index,'encode_phase' ] = 'early_extinction'
        self.df.encode_phase = pd.Categorical(self.df.encode_phase,self.encoding_phases,ordered=True)
        self.df.roi = pd.Categorical(self.df.roi,self.rois,ordered=True)
        self.df.high_confidence_accuracy = (self.df.high_confidence_accuracy == 'H').astype(int)
        self.df.low_confidence_accuracy = (self.df.low_confidence_accuracy == 'H').astype(int)
        if self.hemi: self.df['hemi'] = self.df.roi.apply(lambda x: x[0])

    def load_cross_mats(self):
        in_mats = {}
        for sub in self.subs:
            subj = bids_meta(sub)
            with open(os.path.join('rsa_results',subj.fsub,'cross_mats.p'),'rb') as file:
                in_mats[sub] = pickle.load(file)
        mats = {}
        for roi in self.rois:
            print(roi)
            mats[roi] = np.stack([in_mats[sub][roi]['ers'] for sub in in_mats])
        self.mats = mats

    def load_mem_mats(self):
        in_mats = {}
        for sub in self.subs:
            subj = bids_meta(sub)
            with open(os.path.join('rsa_results',subj.fsub,'mem_mats.p'),'rb') as file:
                in_mats[sub] = pickle.load(file)
        self.mem_mats = in_mats

    def vis_cross_mat(self,rois):
        matlines = {
        'CS+':{
                        'linestyle':'-',                
                         'baseline':{'color':'black',
                                     'trial':[0,24]},
                'fear_conditioning':{'color':'purple',
                                     'trial':[24,48]},                     
                       'extinction':{'color':'green',
                                     'trial':[48,72]}
                 },
         'CS-':{
                        'linestyle':':',
                         'baseline':{'color':'black',
                                     'trial':[72,96]},
                'fear_conditioning':{'color':'purple',
                                     'trial':[96,120]},                     
                       'extinction':{'color':'green',
                                     'trial':[120,144]}
                 }
         }

        if type(rois) is str: rois = [rois]
        # for roi in self.rois:
        for roi in rois:
            mat = self.mats[roi].mean(axis=0)
            mat = np.tanh(mat)
            mask = np.zeros_like(mat)
            mask[np.triu_indices_from(mask)] = True
            
            # _cmap = sns.diverging_palette(220, 20, sep=20, as_cmap=True)
            _cmap = sns.diverging_palette(33,282,center='dark',as_cmap=True)
            fig, ax = plt.subplots()
            sns.heatmap(mat,cmap='twilight_shifted',ax=ax,mask=mask,square=True,center=0)

            ax.set_title(self.group+' - '+roi+'_%s__%s'%(np.round(mat.max(),3),np.round(mat.min(),3)))

            ax.xaxis.set_major_locator(MultipleLocator(24))
            ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.yaxis.set_major_locator(MultipleLocator(24))
            ax.yaxis.set_major_formatter(ScalarFormatter())

            #lines to block out ERS from encoding and retrieval
            # ax.hlines(144,0,144,'grey')
            # ax.vlines(144,144,288,'grey')
           
            # for phase in self.encoding_phases:
            #     for con in self.conditions:
            #         lower = matlines[con][phase]['trial'][0]
            #         upper = matlines[con][phase]['trial'][1]
            #         color = matlines[con][phase]['color']
            #         ls    = matlines[con]['linestyle']
                    
            #         ax.hlines([upper,upper+144,upper+144,upper+120],
            #                   [lower,lower+144,lower,lower],
            #                   [upper,upper+144,upper,upper],
            #                  color,linestyle=ls)
            #         ax.vlines([lower,lower+144,lower,upper],
            #                   [upper,upper+144,upper+144,lower+144],
            #                   [lower,lower+144,lower+144,upper+144],
            #                  color,linestyle=ls)
    def mat_stats(self):
        stats = {}
        bc_overlap = {}
        for roi in self.rois:
            stats[roi] = {}
            bc_overlap[roi] = {}
            mat = self.mats[roi] #put everything back into pearson's r space before computing spearmans corr
            stats[roi]['bc_CS+'] = np.zeros(mat.shape[0])
            stats[roi]['bc_CS-'] = np.zeros(mat.shape[0])
            bc_overlap[roi]['CS+'] = {}
            bc_overlap[roi]['CS-'] = {}
            for i in range(mat.shape[0]):
                m = mat[i,:,:]
                # stats[roi]['bc_CS+'][i] = np.arctanh(spearmanr(m[0:24,24:48],m[144:168,168:192],axis=None)[0])#encoding baseline-to-conditioning
                # stats[roi]['bc_CS-'][i] = np.arctanh(spearmanr(m[72:96,96:120],m[216:240,240:264],axis=None)[0])#memory baseline-to-conditioning
                stats[roi]['bc_CS+'][i] = (m[0:24,24:48] - m[144:168,168:192]).mean()#encoding baseline-to-conditioning
                stats[roi]['bc_CS-'][i] = (m[72:96,96:120] - m[216:240,240:264]).mean()#memory baseline-to-conditioning
                bc_overlap[roi]['CS+'][self.subs[i]] = (m[0:24,24:48] - m[144:168,168:192]).mean(axis=0)
                bc_overlap[roi]['CS-'][self.subs[i]] = (m[72:96,96:120] - m[216:240,240:264]).mean(axis=0)

        self.stats = stats
        self.bc_overlap = bc_overlap

    def inspect(self):
        df = self.df.groupby(['subject','roi','encode_phase','trial_type']
                    ).mean(
                    ).reset_index()
        
        if self.hemi:
            df['hemi'] = df.roi.apply(lambda x: x[0])
            for hemi in ['r','l']:
                dat = df[df.hemi == hemi]
                dat.roi = pd.Categorical(dat.roi,dat.roi.unique(),ordered=True)
                g = sns.catplot(data=dat,y='rsa',x='trial_type',
                row='roi',col='encode_phase',
                kind='point',palette=cpal)
                for ax in g.axes.ravel(): ax.hlines(0,ax.get_xlim()[0],ax.get_xlim()[1])
        else:

            g = sns.catplot(data=df,y='rsa',x='trial_type',
                            row='roi',col='encode_phase',
                            kind='point',palette=cpal)
            for ax in g.axes.ravel(): ax.hlines(0,ax.get_xlim()[0],ax.get_xlim()[1])

    def graph(self):
        # if self.group == 'between': df = self.df #this doesn't really make sense with how i plan to graph the results
        # else: df = self.df.query('group == @self.group')
        df = self.df.groupby(['trial_type','encode','roi','subject']).mean()
        paired = df.loc['CS+'] - df.loc['CS-']
        paired.reset_index(inplace=True);df.reset_index(inplace=True)
        
        # emo_paired = self.df.query('group == @self.group')
        emo_paired = self.df.groupby(['encode','trial_type','roi','subject']).mean()
        emo_paired = emo_paired.loc['fear_conditioning'] - emo_paired.loc['extinction']
        emo_paired.reset_index(inplace=True)

        # sns.set_context('talk');sns.set_style('ticks')
        # fig, ax = plt.subplots(3,1,sharey=True)
        # for i, phase in enumerate(df.encode.unique().categories):
            
            #CS+ and CS- pointplot
            # sns.pointplot(data=df.query('encode == @phase'),x='roi',y='rsa',hue='trial_type',
            #               palette=cpal,capsize=.08,join=False,dodge=True,ax=ax[i])

            #CS+ and CS- swarmplot
            # sns.swarmplot(data=df.query('encode == @phase'),x='roi',y='rsa',hue='trial_type',
            #               palette=cpoint,linewidth=2,edgecolor='black',size=8,ax=ax[i])
            
            #CS+ and CS- boxplot
            # sns.boxplot(data=df.query('encode == @phase'),x='roi',y='rsa',hue='trial_type',
            #             palette=cpal,dodge=True,ax=ax[i])
            
            # ax[i].set_title('Phase = %s'%(phase))
            # ax[i].legend_.remove()

        #Paired, CS+ - CS- pointplot
        if self.ext_split: phase_pal = sns.color_palette(['black','darkmagenta','lightgreen','seagreen'],desat=.75)
        else:              phase_pal = sns.color_palette(['black','darkmagenta','lightgreen','seagreen'],desat=.75)
        
        fig, ax = plt.subplots()
        sns.pointplot(data=paired,x='roi',y='rsa',hue='encode',palette=phase_pal
                      ,capsize=.08,join=False,dodge=True,ax=ax)
        ax.hlines(y=0,xmin=ax.get_xlim()[0],xmax=ax.get_xlim()[1])
        plt.title('Relative RSA\n\n[CS+ - CS-]')
        sns.despine(ax=ax)

        fig, ax = plt.subplots()
        sns.pointplot(data=emo_paired,x='roi',y='rsa',hue='trial_type',
                      palette=cpal,dodge=True,capsize=.08,join=False,ax=ax)
        ax.hlines(y=0,xmin=ax.get_xlim()[0],xmax=ax.get_xlim()[1])
        plt.title('Relative RSA\n\n[Conditioning - Extinction]')
        sns.despine(ax=ax)

        for v in [0,1]:
            df = self.df[self.df.acc == v]
            paired = df.groupby(['trial_type','encode','roi','subject']).mean()
            paired = paired.loc['CS+'] - paired.loc['CS-']
            paired.reset_index(inplace=True);paired.reset_index(inplace=True)
            fig, ax = plt.subplots()
            sns.pointplot(data=paired,x='roi',y='rsa',hue='encode',palette=phase_pal
                      ,capsize=.08,join=False,dodge=True,ax=ax)
            ax.hlines(y=0,xmin=ax.get_xlim()[0],xmax=ax.get_xlim()[1])
            plt.title('Relative RSA\n\n[CS+ - CS-] Memory = %s'%(v))
            sns.despine(ax=ax)

            emo_paired = df.groupby(['encode','trial_type','roi','subject']).mean()
            emo_paired = emo_paired.loc['fear_conditioning'] - emo_paired.loc['extinction']
            emo_paired.reset_index(inplace=True)
            fig, ax = plt.subplots()
            sns.pointplot(data=emo_paired,x='roi',y='rsa',hue='trial_type',
                          palette=cpal,dodge=True,capsize=.08,join=False,ax=ax)
            ax.hlines(y=0,xmin=ax.get_xlim()[0],xmax=ax.get_xlim()[1])
            plt.title('Relative RSA\n\n[Conditioning - Extinction] Memory = %s'%(v))
            sns.despine(ax=ax)


    def cstrial(self):
        idx = pd.IndexSlice
        if self.hemi:
            for hemi in self.df.hemi.unique():
                df = self.df.set_index(['hemi','encode_phase','trial_type','roi']).sort_index()
                fig, ax = plt.subplots(1,len(self.encoding_phases),sharey=True)
                for i, phase in enumerate(self.encoding_phases):
                    dat = df.loc[(hemi,phase,'CS+')].reset_index()
                    dat.roi = pd.Categorical(dat.roi,dat.roi.unique())
                    sns.lineplot(data=dat,x='csp_trial',y='rsa',ax=ax[i],hue='roi',
                        palette='husl')
                    if i !=0: ax[i].legend_.remove()
                    ax[i].hlines(y=0,xmin=ax[i].get_xlim()[0],xmax=ax[i].get_xlim()[1],color='grey',linestyle='--')

                plt.sup_title(hemi)
        else:

            df = self.df.set_index(['encode_phase','trial_type','roi','subject']).sort_index()
            df['block'] = 0
            for con in self.conditions:
                df = df.sort_values(by='%s_trial'%(self.conditions[con]))
                for phase in self.encoding_phases:
                    for roi in self.rois:
                        for sub in self.subs:
                            df.loc[(phase,con,roi,sub),'block'] = np.repeat(range(1,7),4)
            df = df.reset_index()
            df = df.groupby(['trial_type','block','encode_phase','roi','subject']).mean()
            df = (df.loc['CS+','rsa'] - df.loc['CS-','rsa']).reset_index()
            df = df[df.roi.isin(['sgACC','rACC'])]
            df.roi = pd.Categorical(df.roi,df.roi.unique())
            df = df.set_index('encode_phase').sort_index()
            self.csdf = df

            fig, ax = plt.subplots(1,len(self.encoding_phases),sharey=True)
            for i, phase in enumerate(self.encoding_phases):
                dat = df.loc[phase].reset_index()
                
                sns.lineplot(data=dat,x='block',y='rsa',hue='roi',ax=ax[i],
                    palette='husl')
                if i!=0:ax[i].legend_.remove()
                xl=ax[i].get_xlim();ax[i].hlines(0,xl[0],xl[1],color='grey',linestyle='--')
            # if not self.HCA:
            # comp = df.reset_index().set_index(['roi','block','encode'])
            # comp = (comp.loc['dACC','rsa'] - comp.loc['mOFC','rsa']).reset_index().set_index('encode')
            # fig, ax = plt.subplots(1,len(self.encoding_phases),sharey=True)
            # for i, phase in enumerate(self.encoding_phases):
            #     sns.lineplot(data=comp.loc[phase],x='block',y='rsa',ax=ax[i])
            #     xl=ax[i].get_xlim();ax[i].hlines(0,xl[0],xl[1],color='grey',linestyle='--')


    def roi_logreg(self,acc='high_confidence_accuracy'):
        logreg = LogisticRegression(solver='lbfgs')
   
        def boot_roi_logreg(bdf,n_boot=1000):
            bdf = bdf.set_index('subject')
            boot_res = np.zeros(n_boot)
            for i in range(n_boot):
                _samp = np.random.choice(self.subs,24)    
                _y = bdf.loc[_samp,acc].values
                _X = bdf.loc[_samp,'rsa'].values.reshape(-1,1) #get it ready for logreg
                
                logreg.fit(_X,_y)
                boot_res[i] = logreg.coef_[0][0]
            return boot_res

        df = self.df.set_index(['encode_phase','roi','trial_type']).sort_index()
        betas = {}
        for phase in self.encoding_phases:
            betas[phase] = {};print(phase)
            for roi in ['sgACC','rSMA']:
                betas[phase][roi] = {};print(roi)
                for con in self.conditions: 
                    print(con);betas[phase][roi][con] = boot_roi_logreg(df.loc[(phase,roi,con)])

        self.betas = betas
    def bc_logreg(self,acc='high_confidence_accuracy'):
        logreg = LogisticRegression(solver='lbfgs')
        
        #need to match up the values to each baseline trial
        df = self.df.set_index(['encode','roi','trial_type','subject']).sort_index()
        df = df.loc['baseline']
        df['bc_overlap'] = 0
        for con in self.conditions:
            df = df.sort_values(by='%s_trial'%(self.conditions[con]))
            for roi in self.rois:
                for sub in self.subs:
                    df.loc[(roi,con,sub),'bc_overlap'] = self.bc_overlap[roi][con][sub]
        df = df.sort_index()
        def boot_bc_logreg(bdf,n_boot=1000):
            boot_res = np.zeros(n_boot)
            for i in range(n_boot):
                _samp = np.random.choice(self.subs,24)    
                _y = bdf.loc[_samp,acc].values
                _X = bdf.loc[_samp,'bc_overlap'].values.reshape(-1,1) #get it ready for logreg
                
                logreg.fit(_X,_y)
                boot_res[i] = logreg.coef_[0][0]
            return boot_res
        
        bc_betas = {}
        for roi in self.rois:
            bc_betas[roi] = {};print(roi)
            for con in self.conditions: 
                print(con);bc_betas[roi][con] = boot_bc_logreg(df.loc[(roi,con)])
        self.bc_betas = bc_betas


    def cross_logreg(self,rois=None,acc='hc_acc'):
        logreg = LogisticRegression(solver='lbfgs')
        if rois is None: rois = self.rois
        def boot_cross_logreg(bdf,n_boot=1000):
        
            boot_res = np.zeros( (n_boot, len(rois)) )
            for i in range(n_boot):
                _samp = np.random.choice(self.subs,len(self.subs))
                _y = bdf.loc[_samp,acc].values
                _X = bdf.loc[_samp,rois].values
                logreg.fit(_X,_y)
                boot_res[i,] = logreg.coef_[0]
            return boot_res
        
        #this chunk of unintelligible pandas code takes in the long form data frame,
        #groups everything by the unique trials and puts the RSA values in wide format
        df = self.df
        df = df.drop(columns=['memcond','phase','response','phase','csp_trial','csm_trial']
            ).set_index(['subject','stim','encode','trial_type','roi']
            ).unstack(level=-1)['rsa'].rename(columns=str
            ).reset_index(
            ).set_index(['subject','stim'])
        #we don't need repeated values for the memory, just grab the first one for each trial    
        df[acc] = self.df.groupby(['subject','stim']).first()[acc]
        #here we set up for the bootstrapping logreg
        df = df.reset_index(
            ).set_index(['encode','trial_type']
            ).sort_index(
            ).drop(columns='stim')

        betas = {}
        for phase in self.encoding_phases:
            betas[phase] = {};print(phase)
            for con in self.conditions:
                print(con);betas[phase][con] = boot_cross_logreg(df.loc[(phase,con)].set_index('subject'))

        self.betas = betas

    def vis_logreg(self):

        fig, ax = plt.subplots(1,4,sharey=True)

        for i, phase in enumerate(self.encoding_phases):
            g = self.betas[phase]
            gdf = pd.DataFrame.from_dict(g,orient='index').reset_index(
                ).melt(id_vars='index'
                ).rename(columns={'index':'roi','variable':'trial_type','value':'beta'})

            gdf = gdf.set_index(['roi', 'trial_type'])['beta'].apply(pd.Series).stack(
                ).reset_index().drop(columns='level_2').rename(columns={0:'beta'})
            # if not self.hemi: gdf.roi = pd.Categorical(gdf.roi,self.rois,ordered=True)
            gstats = gdf.groupby(['trial_type','roi']).mean()
            gstats['conf']  = gdf.groupby(['trial_type','roi']).apply(lambda x: np.percentile(x,[5,95]))
            
            sns.pointplot(data=gdf,x='roi',y='beta',hue='trial_type',
                      palette=cpal,dodge=True,capsize=.08,join=False,ax=ax[i],ci=None)
            ax[i].set_title(phase)
            ax[i].hlines(y=0,xmin=ax[i].get_xlim()[0],xmax=ax[i].get_xlim()[1])
            ax[i].set_xticklabels(ax[i].get_xticklabels(),rotation=90)


            for j, roi in enumerate(gdf.roi.unique()):
                x = ax[i].get_xticks()[j]
                ax[i].vlines(x-.05,ymin=gstats.loc[('CS+',roi),'conf'][0], ymax=gstats.loc[('CS+',roi),'conf'][1], color=cpal[0], linewidth=3)
                ax[i].vlines(x+.05,ymin=gstats.loc[('CS-',roi),'conf'][0], ymax=gstats.loc[('CS-',roi),'conf'][1], color=cpal[1], linewidth=3)
                # ax[i].errorbar(x+.05,gstats.loc[('CS-'),'beta'], yerr=np.concatenate(gstats.loc[('CS-'),'conf'].values).reshape(5,2).transpose(), color='black', linewidth=3,capsize=5,capthick=3)
            if i != 0: ax[i].legend_.remove()
            # ax[i] = sns.violinplot(data=gdf,x='roi',y='beta',hue='trial_type',scale='count',
            #           palette=cpal,split=True,join=False,ax=ax[i],cut=True)

    def vis_cross_logreg(self,rois=None):


        if rois is None: rois = self.rois
        fig, ax = plt.subplots(1,4,sharey=True)

        for i, phase in enumerate(self.encoding_phases):
            g = self.betas[phase]
            gdf = {}
            for con in self.conditions:
                gdf[con] = pd.DataFrame(g[con])
                gdf[con].columns = rois
                gdf[con]['trial_type'] = con
            gdf = pd.concat(gdf.values())

            gdf = gdf.melt(id_vars='trial_type'
                ).rename(columns={'variable':'roi','value':'beta'})

            #calculate the stats
            gstats = gdf.groupby(['trial_type','roi']).mean()
            gstats['conf']  = gdf.groupby(['trial_type','roi']).apply(lambda x: np.percentile(x,[5,95]))
            
            #graphing starts here
            sns.pointplot(data=gdf,x='roi',y='beta',hue='trial_type',
                      palette=cpal,dodge=True,capsize=.08,join=False,ax=ax[i],ci=None)
            ax[i].set_title(phase)
            ax[i].hlines(y=0,xmin=ax[i].get_xlim()[0],xmax=ax[i].get_xlim()[1])
            ax[i].set_xticklabels(ax[i].get_xticklabels(),rotation=45)

            for j, roi in enumerate(gdf.roi.unique()):
                x = ax[i].get_xticks()[j]
                ax[i].vlines(x-.05,ymin=gstats.loc[('CS+',roi),'conf'][0], ymax=gstats.loc[('CS+',roi),'conf'][1], color=cpal[0], linewidth=3)
                ax[i].vlines(x+.05,ymin=gstats.loc[('CS-',roi),'conf'][0], ymax=gstats.loc[('CS-',roi),'conf'][1], color=cpal[1], linewidth=3)
            if i !=0: ax[i].legend_.remove()

def vox_count():
    vox = {}
    for sub in all_sub_args:
        subj = bids_meta(sub)
        vox[sub] = {}
        for roi in ['hpc','mOFC','dACC','ins','amyg']:
            mask = get_data(os.path.join(subj.masks,'%s_mask.nii.gz'%(roi)))
            vox[sub][roi] = np.where(mask == 1)[0].shape[0]
    df = pd.DataFrame.from_dict(vox,orient='index'
                    ).reset_index(
                    ).melt(id_vars='index'
                    ).rename(columns={'index':'subject', 'variable':'roi', 'value':'nvox'})
    fig, ax = plt.subplots()
    sns.boxplot(data=df,x='roi',y='nvox',ax=ax)
    sns.despine(ax=ax,trim=True)

def resp_count():
    df = pd.read_csv(os.path.join(data_dir,'group_ER','all_template_df.csv'),index_col=0)
    df = df.groupby(['subject','encode','trial_type']).sum().reset_index()
    df['group'] = df.subject.apply(lgroup)
    df.encode = pd.Categorical(df.encode,categories=['baseline','fear_conditioning','extinction'],ordered=True)
    df[['hc_acc','acc']] /= 24
    for acc in ['hc_acc','acc']:
        fig, ax = plt.subplots(1,2,sharey=True)
        for i, group in enumerate(['control','ptsd']):
            dat = df.query('group == @group')
            sns.boxplot(data=dat,x='encode',y=acc,hue='trial_type',palette=cpal,ax=ax[i])
            # sns.swarmplot(data=dat,x='encode',y=acc,hue='trial_type',palette=cpal,ax=ax[i],
            #             linewidth=2,edgecolor='black',dodge=True,size=5)
        # ax.set_yticks(range(0,25,4))
            ax[i].yaxis.set_major_locator(MultipleLocator(.2))
            ax[i].yaxis.set_minor_locator(MultipleLocator(.1))
            ax[i].set_ylim(0,1.01)
            ax[i].legend_.remove()
            ax[i].set_title(group)

def copy_out():
    # from fg_config import *
    import os
    from glob import glob
    out = os.path.join(SCRATCH,'rsa_results');mkdir(out)
    for sub in all_sub_args:
        subj = bids_meta(sub)
        sub_out = os.path.join(out,subj.fsub);mkdir(out)
        os.system('cp -R %s %s'%(subj.rsa,sub_out))
    q = glob('/scratch/05426/ach3377/rsa_results/sub-FC***/sl_er.p')