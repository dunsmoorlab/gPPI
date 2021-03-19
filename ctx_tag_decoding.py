# import nibabel as nib
# import matplotlib.pyplot as plt
# import seaborn as sns

# from fg_config import *

# from nilearn.image import concat_imgs, get_data
# from sklearn.feature_selection import SelectKBest, f_classif
# from sklearn.linear_model import LogisticRegression
# from sklearn.pipeline import Pipeline
# from sklearn.model_selection import cross_val_score
# from sklearn.metrics import confusion_matrix, classification_report, precision_recall_fscore_support, roc_auc_score
# from sklearn.preprocessing import label_binarize

def relative_evidence(x): #takes in the decision function from a classifier and converts to non-multinomial evidence values
    return np.reciprocal(np.exp(np.copy(x)*-1)+1)

class sub_xval():

    def __init__(self,sub,masks=None):

        print(sub)
        self.subj = bids_meta(sub)
        self.data, self.labels = self.load_loc()
        self.decode(masks)
        self.prep_output()
        

    def load_loc(self):
        if self.subj.num == 107: runs = [1]
        else: runs = [1,2]

        data = concat_imgs([nib.load(f'{self.subj.beta}/localizer_run-0{run}_beta.nii.gz') for run in runs]).get_fdata()
        labels = pd.concat([pd.read_csv(f'{self.subj.subj_dir}/ses-2/func/{self.subj.fsub}_ses-2_task-localizer_run-0{run}_events.tsv',sep='\t') for run in runs]).reset_index(drop=True)
        self.orig_labels = labels.trial_type
        #this is what we did for FearCon so we are going to start with just outdoor vs. scrambled
        # trial_mask = labels.trial_type.isin(['outdoor','indoor','scrambled','rest','animal','tool'])
        # trial_mask = labels.trial_type.isin(['outdoor','indoor','scrambled'])
        
        #downsample and combine indoor/outdoor scenes
        indoor = labels.index[labels.trial_type == 'indoor'][::2]
        outdoor = labels.index[labels.trial_type == 'outdoor'][::2]
        scrambled = labels.index[labels.trial_type == 'scrambled']
        # animal = labels.index[labels.trial_type == 'animal'][::2]
        # tool = labels.index[labels.trial_type == 'tool'][::2]
        
        scenes = np.concatenate((indoor,outdoor))
        # stims = np.concatenate((animal,tool))

        labels.loc[scenes,'trial_type'] = 'scene'
        # labels.loc[stims,'trial_type'] = 'stim'
        trial_mask = np.concatenate((scenes,scrambled))#add stims here if necessary
        trial_mask.sort()


        data = data[:,:,:,trial_mask]
        labels = labels.loc[trial_mask,'trial_type']
        return data, labels
 
    def decode(self,masks):
        self.clf = Pipeline([ ('anova', SelectKBest(f_classif)), ('logreg', LogisticRegression(multi_class='multinomial')) ])

        self.scores = {}
        self.cmats = {}
        self.ev = {}
        self.xval_groups = np.repeat([1,2],self.labels.shape[0]/2)
        

        for roi in masks:
            roi_data = apply_mask(mask=get_data(f'{self.subj.masks}/{roi}.nii.gz'), target=self.data)
            self.scores[roi] = self.run_xval(data=roi_data,labels=self.labels)
            self.cmats[roi] = self.confusion_matrix(data=roi_data,labels=self.labels)
            self.ev[roi] = self.collect_ev(data=roi_data,labels=self.labels) 

    def run_xval(self, data=None, labels=None):
        if data.shape[1] < 500: k = 'all'
        else: k = 500
        self.clf.set_params(anova__k=k)
        
        score = 'roc_auc' if len(np.unique(labels)) == 2 else 'roc_auc_ovr'
        xval_score = cross_val_score(self.clf, data, labels, groups=self.xval_groups, cv=2, scoring=score)
        # proba_res = { phase: self.clf.predict_proba(test_data[phase]) for phase in test_data }

        return xval_score.mean()

    def confusion_matrix(self,data=None,labels=None):

        cmats = {}
        for group in np.unique(self.xval_groups):
            train_data = data[self.xval_groups == group]
            test_data = data[self.xval_groups != group]
            train_labels = labels[self.xval_groups == group]

            prediction = self.clf.fit(train_data,train_labels).predict(test_data)
            cmats[group] = confusion_matrix(train_labels,prediction,labels=self.clf.classes_,normalize='true')

        return np.mean([cmats[1],cmats[2]],axis=0)

    def collect_ev(self, data=None, labels=None):
        
        evs = {}
        for group in np.unique(self.xval_groups):
            train_data = data[self.xval_groups == group]
            test_data = data[self.xval_groups != group]
            train_labels = labels[self.xval_groups == group]

            prediction = self.clf.fit(train_data,train_labels).decision_function(test_data)
            evs[group] = relative_evidence(prediction)

        ev = np.concatenate((evs[1],evs[2]),axis=0)
        if len(self.clf.classes_) == 2: ev = np.vstack((ev,1-ev)).T
        return ev

    def prep_output(self):
        self.df = pd.DataFrame.from_dict(self.scores,orient='index'
                             ).reset_index(
                             ).rename(columns={0:'auc','index':'roi'})
        self.df['subject'] = self.subj.num

        self.ev_df = pd.DataFrame(self.ev['ppa'],columns=self.clf.classes_)
        self.ev_df['true_label'] = self.orig_labels
        self.ev_df = self.ev_df.melt(id_vars='true_label',value_vars=self.clf.classes_,var_name='category',value_name='proba')
        self.ev_df['subject'] = self.subj.num

class group_xval():

    def __init__(self,masks=None):

        self.collect_sub_dat(masks=masks)
        self.vis_results(df=self.df)
        # self.vis_cmats(masks=masks)

    def collect_sub_dat(self,masks=None):

        dfs = {}
        cmats = {}
        ev_dfs = {}
        for sub in all_sub_args:
            subxval = sub_xval(sub,masks)
            dfs[sub] = subxval.df
            cmats[sub] = subxval.cmats
            ev_dfs[sub] = subxval.ev_df
        self.classes = subxval.clf.classes_ #should be stable across everyone
        self.df = pd.concat(dfs.values())
        self.cmats = cmats
        self.ev_df = pd.concat(ev_dfs.values())
        
        self.df['group'] = self.df.subject.apply(lgroup)
        self.ev_df['group'] = self.ev_df.subject.apply(lgroup)

    def vis_results(self,df):
        fig, ax = plt.subplots()
        sns.pointplot(data=df,x='roi',y='auc',ax=ax,color='black')
        sns.swarmplot(data=df,x='roi',y='auc',ax=ax)
        ax.hlines(.5,*ax.get_xlim(),color='black')

        fig, ax = plt.subplots()
        sns.pointplot(data=df,x='roi',y='auc',ax=ax,hue='group',dodge=True,palette=gpal)
        sns.swarmplot(data=df,x='roi',y='auc',ax=ax,hue='group',palette=gpal)
        ax.hlines(.5,*ax.get_xlim(),color='black')
        ax.legend_.remove()

    def vis_cmats(self,masks=None):
        for roi in masks:
            avg = np.array([self.cmats[sub][roi] for sub in self.cmats])
            mat = avg.mean(axis=0)
            fig, ax = plt.subplots()
            sns.heatmap(mat,xticklabels=self.classes,yticklabels=self.classes,cmap='Purples',annot=True,cbar=False,ax=ax)
            ax.set_xlabel('predicted')
            ax.set_ylabel('true')
            plt.tight_layout()

    def vis_ev_df(self,df):
        df.true_label = df.true_label.apply(lambda x: 'stim' if x in ['animal','tool'] else x)
        df.true_label = df.true_label.apply(lambda x: 'scene' if x in ['indoor','outdoor'] else x)
        df = df.groupby(['true_label','category','group','subject']).mean().reset_index()
        fig, ax = plt.subplots()
        sns.barplot(data=df,x='true_label',y='proba',hue='category',ax=ax)
        ax.set_ylabel('Non-multinomial relative evidence')
        ax.legend(loc='center left',bbox_to_anchor=(1,.5))
        plt.tight_layout()

class sub_mem_decode():

    def __init__(self, sub, masks=None):

        print(sub)
        self.subj = bids_meta(sub)
        self.loc_data, self.loc_labels = sub_xval.load_loc(self)
        self.mem_data, self.mem_labels = self.load_mem_dat()
        self.wrap_decode(masks=masks)
        self.prep_output()

    def load_mem_dat(self):
        runs = [1,2,3]

        data = concat_imgs([nib.load(f'{self.subj.beta}/memory_run-0{run}_beta.nii.gz') for run in runs]).get_fdata()
        labels = pd.concat([pd.read_csv(f'{self.subj.subj_dir}/ses-2/func/{self.subj.fsub}_ses-2_task-memory_run-0{run}_events.tsv',sep='\t') for run in runs]).reset_index(drop=True)

        return data, labels

    def wrap_decode(self,masks):
        self.clf = Pipeline([ ('anova', SelectKBest(f_classif)), ('logreg', LogisticRegression(multi_class='multinomial')) ])

        self.proba_res = {}
        for roi in masks:
            mask_data = get_data(f'{self.subj.masks}/{roi}.nii.gz')
            train_data = apply_mask(mask=mask_data, target=self.loc_data)
            test_data = apply_mask(mask=mask_data, target=self.mem_data)
            self.proba_res[roi] = self.decode(roi=roi,train_data=train_data, train_labels=self.loc_labels, test_data=test_data, test_labels=self.mem_labels)

    def decode(self,roi,train_data,train_labels,test_data,test_labels):
        if train_data.shape[1] < 500: k = 'all'
        else: k = 500
        self.clf.set_params(anova__k=k)
        self.clf.fit(train_data,train_labels)
        # proba_res = self.clf.predict_proba(test_data)
        proba_res = relative_evidence(self.clf.decision_function(test_data))

        if len(self.clf.classes_) == 2:
            test_labels[self.clf.classes_[0]] = proba_res
            test_labels[self.clf.classes_[1]] = 1 - proba_res
        else:
            for class_ in self.clf.classes_:
                test_labels[class_] = proba_res[:,np.where(self.clf.classes_ == class_)[0][0]]
        
        test_labels['roi'] = roi
        
        return test_labels

    def prep_output(self):
        self.df = pd.concat(self.proba_res.values()).reset_index(drop=True)
        self.df['subject'] = self.subj.num

class group_mem_decode():

    def __init__(self,masks=None):

        self.collect_sub_dfs(masks=masks)
        # self.vis_results(df=self.df)

    def collect_sub_dfs(self,masks=None):

        dfs = {}
        for sub in all_sub_args:
            subxval = sub_mem_decode(sub,masks)
            dfs[sub] = subxval.df
        self.classes_ = subxval.clf.classes_
        self.df = pd.concat(dfs.values())

        self.df['group'] = self.df.subject.apply(lgroup)

    def vis_results(self,df):
        # class_ = 'scene'
        # q = d.df.groupby(['roi','encode_phase','trial_type','subject'])[class_].mean().reset_index()
        # q['group'] = q.subject.apply(lgroup)

        # sns.catplot(data=q,col='encode_phase',x='trial_type',hue='group',y=class_,palette=gpal,kind='swarm')
        # fig, ax = plt.subplots()
        # sns.pointplot(data=q,x='encode_phase',y=class_,ax=ax,hue='group',dodge=True,palette=gpal)
        # sns.swarmplot(data=q,x='encode_phase',y=class_,ax=ax,hue='group',palette=gpal,dodge=True)
        # ax.legend_.remove()


        # sns.displot(y=d.df.proba,kind='ecdf',rug=True)

        df = df.copy()
        df = df.dropna(subset=['response'])
        df = df.drop(columns=['roi','onset','duration','memory_condition','response','low_confidence_accuracy','high_confidence_accuracy','response_time'])
        df = df.rename(columns={'encode_phase':'phase','trial_type':'condition'})
        df = df.melt(id_vars=['group','phase','condition','subject','stimulus'],value_vars=self.classes_,var_name='category',value_name='proba')
        # sns.displot(data=df,y='proba',hue='category',kind='ecdf')

        sns.catplot(data=df.groupby(['group','subject','category','condition','phase']).mean().reset_index(),x='condition',y='proba',hue='category',row='group',col='phase',kind='bar')

        df = df.groupby(['category','group','phase','condition','subject']).mean().sort_index().loc['scene']
        # df.to_csv('ctx_ev.csv')
        df = df.reset_index()
        sns.catplot(data=df,x='condition',y='proba',hue='phase',col='group',palette=phase_pal,kind='swarm',dodge=True,hue_order=phase3+['foil'])
'''copy localizer events to tacc again'''
def copy_out():
    dest = '/Users/ach3377/Desktop/loc_data';mkdir(dest)
    source = f'{SCRATCH}/loc_data'
    # for sub in all_sub_args:
    for sub in [21,116]:
        subj = bids_meta(sub)
        # for run in [1,2]:
        for run in [2]:
            if sub == 107 and run == 2:
                pass
            else:
                # events = f'{subj.fsub}_ses-2_task-localizer_run-0{run}_events.tsv'
                # os.system(f'cp {subj.subj_dir}/ses-2/func/{events} {dest}/{events}')
                # os.system(f'cp {dest}/{events} {subj.subj_dir}/ses-2/func/{events}')

                bold = f'{subj.fsub}_ses-2_task-localizer_run-0{run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'
                # os.system(f'cp {subj.func}/{bold} {dest}/{bold}')
                os.system(f'cp {source}/{bold} {subj.func}/{bold}')

    _in = f'{subj.prep_dir}/ses-2/func/{subj.fsub}_ses-2_task-localizer_run-0{run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'
    _out = f'{subj.func}/{subj.fsub}_ses-2_task-localizer_run-0{run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'
    os.system(f'fslmaths {_in} -mas {subj.refvol_mask} {_out}')
def maskthing(subs):
    import os
    for sub in subs:
        print(sub)
        fsub = 'sub-FC{0:0=3d}'.format(sub)

        # os.system(f'flirt -in /Users/ach3377/Desktop/ppa_masks/group_ppa_mask.nii.gz -ref {subj.refvol_brain} -applyxfm -init {subj.std2ref} -interp nearestneighbour -out {subj.masks}/ppa.nii.gz')

        for category in ['animal','tool']:
            os.system(f'flirt -in /mnt/c/Users/ACH/Desktop/standard/{category}_fg.nii.gz \
                              -ref /mnt/d/fc-bids/derivatives/preproc/{fsub}/reference/boldref_brain.nii.gz \
                              -applyxfm \
                              -init /mnt/d/fc-bids/derivatives/preproc/{fsub}/reference/std2ref.mat -interp nearestneighbour \
                              -out /mnt/d/fc-bids/derivatives/preproc/{fsub}/masks/{category}_fg_mask.nii.gz')

    # for i in ev.index:
    #     g, ph, c, sub, stim, val = ev.loc[i,['group','phase','condition','subject','stimulus','proba']]
    #     p.loc[(g,ph,c,sub,stim),'ev'] = val