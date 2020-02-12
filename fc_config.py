import os
import pandas as pd
import numpy as np

sub_args = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,23,24,25,26]
p_sub_args = [101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117,118, 120, 121, 122, 123, 124, 125]
all_sub_args = sub_args + p_sub_args

subjects = {'control':sub_args,
            'ptsd':p_sub_args,
            'all':all_sub_args}

WORK = '/work/05426/ach3377/lonestar/'
HOME = '/home1/05426/ach3377/'
CODEBASE = HOME + 'CodeBase/'

data_dir = os.path.join(WORK,'fc')
group_glm = os.path.join(WORK,'group_glm')


std_1mm_brain = os.path.join(WORK,'standard','MNI152_T1_1mm_brain.nii.gz')
std_3mm_brain = os.path.join(WORK,'standard','MNI152_T1_3mm_brain.nii.gz')
std_3mm_brain_mask = os.path.join(WORK,'standard','MNI152_T1_3mm_brain_mask.nii.gz')


fsl_betas = {
    'baseline': '/run001_beta.nii.gz',
    'fear_conditioning': '/run002_beta.nii.gz',
    'extinction': '/run003_beta.nii.gz',
    'extinction_recall': '/run004_beta.nii.gz' ,
    'memory_run_1': '/run005_beta.nii.gz',
    'memory_run_2': '/run006_beta.nii.gz',
    'memory_run_3': '/run007_beta.nii.gz',
    'localizer_1': '/run008_beta.nii.gz',
    'localizer_2': '/run009_beta.nii.gz',
}

class meta(object):

	def __init__(self, sub):

		if sub == 'gusbrain':
			self.fs_id = 'gusbrain'
			self.fsub = 'gusbrain'
		elif sub == 'fs':
			self.fs_id = 'fsaverage'
			self.fsub = 'fsaverage'
		else:
			
			self.num = int(sub)
			
			self.fsub = 'Sub{0:0=3d}'.format(self.num)

			self.subj_dir = os.path.join(data_dir,self.fsub)
                        
			self.ref2std = os.path.join(self.subj_dir,'ref2std.mat')

			self.bold_dir = os.path.join(self.subj_dir,'bold')

			self.beta_dir = os.path.join(self.bold_dir,'beta')

			self.anatomy = os.path.join(self.subj_dir, 'anatomy')

			self.orig_brain = os.path.join(self.anatomy, 'orig_brain.nii.gz')

			self.mask = os.path.join(self.subj_dir, 'mask')

			self.roi = os.path.join(self.mask,'roi')

			self.refvol_be = os.path.join(WORK, 'FearCon', self.fsub, 'tacc_temp', 'be_refvol.nii.gz')

			self.reg_dir = os.path.join(self.subj_dir,'reg')

			self.model_dir = os.path.join(self.subj_dir, 'model')

			self.behavior = os.path.join(self.subj_dir,'behavior')
			
			self.sl_dir = os.path.join(self.model_dir,'searchlight')
			if not os.path.exists(self.sl_dir): os.mkdir(self.sl_dir)
			
			self.rsa = os.path.join(self.model_dir,'rsa')
			if not os.path.exists(self.rsa): os.mkdir(self.rsa)

			meta_file = os.path.join(self.subj_dir,'behavior','%s_elog.csv'%(self.fsub))
			if os.path.exists(meta_file):
				self.meta = pd.read_csv(meta_file)

				self.cs_lookup()

	def cs_lookup(self):	
		if self.meta['DataFile.Basename'][0][0] == 'A':
			self.csplus = 'animal'
			self.csminus = 'tool'
		elif self.meta['DataFile.Basename'][0][0] == 'T':
			self.csplus = 'tool'
			self.csminus = 'animal'
