from fg_config import *
import nibabel as nib
from nilearn.image import mean_img, threshold_img

class preproc():

    def __init__(self,sub):

        self.subj = bids_meta(sub)

    def be_t1(self):

        t1 = os.path.join(self.subj.prep_dir,'anat','%s_desc-prepoc_T1w.nii.gz'%(self.subj.fsub))
        t1_mask = os.path.join(self.subj.prep_dir,'anat','%s_desc-brain_mask.nii.gz')

        os.system('cp %s %s'%(t1,self.subj.t1))
        os.system('cp %s %s'%(t1_mask,self.subj.t1_mask))
        os.system('fslmaths %s -mas %s %s'%(self.subj.t1,self.subj.t1_mask,self.subj.t1_brain))

    def refvol(self):

        refs = [nib.load( os.path.join(prep_dir,folder[0],file) )
                        for folder in os.walk(prep_dir)
                        for file in folder[2]
                        if 'T1w_boldref.nii.gz' in file]
        ref = mean_img(refs)
        nib.save(ref,self.subj.refvol)
        
        masks = [nib.load( os.path.join(prep_dir,folder[0],file) )
                        for folder in os.walk(prep_dir)
                        for file in folder[2]
                        if 'T1w_desc-brain_mask.nii.gz' in file]
        mask = mean_img(masks)
        mask = threshold_img(mask,1)
        nib.save(mask,self.subj.refvol_mask)

        os.system('fslmaths %s -mas %s %s'%(self.subj.refvol,self.subj.refvol_mask,self.subj.refvol_brain))
        