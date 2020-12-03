import nibabel as nib
import matplotlib.pyplot as plt
import seaborn as sns

from fg_config import *

from nilearn.image import concat_imgs, get_data
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score
from sklearn.metrics import confusion_matrix, classification_report, precision_recall_fscore_support, roc_auc_score
from sklearn.preprocessing import label_binarize


class sub_xval():

    def __init__(self,sub,masks=None):

        self.subj = bids_meta(sub)
        self.data, self.labels = self.load_loc()
        self.decode(masks)
        self.prep_output()
        

    def load_loc(self):
        if self.subj.num == 107: runs = [1]
        else: runs = [1,2]

        data = concat_imgs([nib.load(f'{self.subj.beta}/localizer_run-0{run}_beta.nii.gz') for run in runs]).get_fdata()
        labels = pd.concat([pd.read_csv(f'{self.subj.subj_dir}/ses-2/func/{self.subj.fsub}_ses-2_task-localizer_run-0{run}_events.tsv',sep='\t') for run in runs]).reset_index()

        #this is what we did for FearCon so we are going to start with just outdoor vs. scrambled
        trial_mask = labels.trial_type.isin(['outdoor','indoor','scrambled','rest','animal','tool'])

        data = data[:,:,:,trial_mask]
        labels = labels.trial_type[trial_mask]
        return data, labels

    def decode(self,masks):
        self.clf = Pipeline([ ('anova', SelectKBest(f_classif)), ('logreg', LogisticRegression(multi_class='ovr')) ])

        self.scores = {}
        self.cmats = {}
        self.xval_groups = np.repeat([1,2],self.labels.shape[0]/2)
        
        for roi in masks:
            roi_data = apply_mask(mask=get_data(f'{self.subj.masks}/{roi}.nii.gz'), target=self.data)
            self.scores[roi] = self.run_xval(data=roi_data,labels=self.labels)
            self.cmats[roi] = self.confusion_matrix(data=roi_data,labels=self.labels)

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

    def prep_output(self):
        self.df = pd.DataFrame.from_dict(self.scores,orient='index'
                             ).reset_index(
                             ).rename(columns={0:'auc','index':'roi'})
        self.df['subject'] = self.subj.num

class group_xval():

    def __init__(self,masks=None):

        self.collect_sub_dat(masks=masks)
        self.vis_results(df=self.df)
        self.vis_cmats(masks=masks)

    def collect_sub_dat(self,masks=None):

        dfs = {}
        cmats = {}
        for sub in all_sub_args:
            subxval = sub_xval(sub,masks)
            dfs[sub] = subxval.df
            cmats[sub] = subxval.cmats
        self.classes = subxval.clf.classes_ #should be stable across everyone
        self.df = pd.concat(dfs.values())
        self.cmats = cmats
        self.df['group'] = self.df.subject.apply(lgroup)

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

class sub_mem_decode():

    def __init__(self, sub, masks=None):

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
        self.clf = Pipeline([ ('anova', SelectKBest(f_classif)), ('logreg', LogisticRegression()) ])

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
        proba_res = self.clf.predict_proba(test_data)

        test_labels['proba'] = proba_res[:,np.where(self.clf.classes_ == 'outdoor')[0][0]]
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

        self.df = pd.concat([sub_mem_decode(sub,masks).df for sub in all_sub_args])
        self.df['group'] = self.df.subject.apply(lgroup)

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

        q = d.df.groupby(['roi','encode_phase','trial_type','subject'])['proba'].mean().reset_index()
        q['group'] = q.subject.apply(lgroup)

        fig, ax = plt.subplots()
        sns.catplot(data=q,col='encode_phase',x='trial_type',hue='group',y='proba',palette=gpal,kind='bar')
        sns.pointplot(data=q,x='roi',y='proba',ax=ax,hue='group',dodge=True,palette=gpal)
        sns.swarmplot(data=q,x='roi',y='proba',ax=ax,hue='group',palette=gpal,dodge=True)
        ax.legend_.remove()


        sns.displot(y=d.df.proba,kind='ecdf',rug=True)





'''copy localizer events to tacc again'''
dest = '/Users/ach3377/Desktop/loc_data';mkdir(dest)
# source = f'{SCRATCH}/loc_data'
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
            os.system(f'cp {subj.func}/{bold} {dest}/{bold}')
            # os.system(f'cp {source}/{bold} {subj.func}/{bold}')

_in = f'{subj.prep_dir}/ses-2/func/{subj.fsub}_ses-2_task-localizer_run-0{run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'
_out = f'{subj.func}/{subj.fsub}_ses-2_task-localizer_run-0{run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'
os.system(f'fslmaths {_in} -mas {subj.refvol_mask} {_out}')
