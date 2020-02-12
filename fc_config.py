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
SCRATCH = '/SCRATCH/05426/ach3377/'
gPPI = HOME + 'gPPI/'

bids_dir = os.path.join(SCRATCH,'fc-bids')
deriv    = os.path.join(bids_dir, 'derivatives')
prep_dir = os.path.join(deriv,'fmriprep')
fs_dir   = os.path.join(deriv,'freesurfer')

# group_glm = os.path.join(WORK,'group_glm')


std_1mm_brain = os.path.join(WORK,'standard','MNI152_T1_1mm_brain.nii.gz')
std_3mm_brain = os.path.join(WORK,'standard','MNI152_T1_3mm_brain.nii.gz')
std_3mm_brain_mask = os.path.join(WORK,'standard','MNI152_T1_3mm_brain_mask.nii.gz')

class bids_meta(object):

	def __init__(self, sub):

		if sub == 'gusbrain':
			self.fs_id = 'gusbrain'
			self.fsub = 'gusbrain'
		elif sub == 'fs':
			self.fs_id = 'fsaverage'
			self.fsub = 'fsaverage'
		else:
			
			self.num = int(sub)
			
			self.fsub = 'sub-FC{0:0=3d}'.format(self.num)

			self.subj_dir = os.path.join(bids_dir,self.fsub)
            self.prep_dir = os.path.join(prep_dir,self.fsub)
            self.fs_dir   = os.path.join(fs_dir,self.fsub)

# def cs_lookup(self):	
# 	if self.meta['DataFile.Basename'][0][0] == 'A':
# 		self.csplus = 'animal'
# 		self.csminus = 'tool'
# 	elif self.meta['DataFile.Basename'][0][0] == 'T':
# 		self.csplus = 'tool'
# 		self.csminus = 'animal'
