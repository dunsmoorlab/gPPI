from fg_config import *
import nibabel as nib
from nilearn.image import mean_img, get_data, threshold_img, new_img_like
from subprocess import Popen

class fmriprep_preproc():

    def __init__(self,sub):

        self.subj = bids_meta(sub)
        self.space = 'MNI152NLin2009cAsym'
        # self.space = 'T1w'
    
    def be_t1(self): #brain extract T1w

        # t1 = os.path.join(self.subj.prep_dir,'anat','%s_desc-preproc_T1w.nii.gz'%(self.subj.fsub))
        # t1_mask = os.path.join(self.subj.prep_dir,'anat','%s_desc-brain_mask.nii.gz'%(self.subj.fsub))
        t1 = os.path.join(self.subj.prep_dir,'anat','%s_space-%s_desc-preproc_T1w.nii.gz'%(self.subj.fsub,self.space))
        t1_mask = os.path.join(self.subj.prep_dir,'anat','%s_space-%s_desc-brain_mask.nii.gz'%(self.subj.fsub,self.space))

        os.system('cp %s %s'%(t1,self.subj.t1))
        os.system('cp %s %s'%(t1_mask,self.subj.t1_mask))
        os.system('fslmaths %s -mas %s %s'%(self.subj.t1,self.subj.t1_mask,self.subj.t1_brain))

    def refvol(self): #create reference bold, bold mask, and mask the ref

        refs = [nib.load( os.path.join(self.subj.prep_dir,folder[0],file) )
                        for folder in os.walk(self.subj.prep_dir)
                        for file in folder[2]
                        if '%s_boldref.nii.gz'%(self.space) in file]
        ref = mean_img(refs)
        nib.save(ref,self.subj.refvol)
        
        masks = [nib.load( os.path.join(self.subj.prep_dir,folder[0],file) )
                        for folder in os.walk(self.subj.prep_dir)
                        for file in folder[2]
                        if '%s_desc-brain_mask.nii.gz'%(self.space) in file]
        mask = mean_img(masks)
        mask = threshold_img(mask,.1)
        nib.save(mask,self.subj.refvol_mask)

        os.system('fslmaths %s -bin %s'%(self.subj.refvol_mask,self.subj.refvol_mask))

        os.system('fslmaths %s -mas %s %s'%(self.subj.refvol,self.subj.refvol_mask,self.subj.refvol_brain))
    
    def be_bold(self): #mask all of the preproc bold files with the new mask

        for folder in os.walk(self.subj.prep_dir):
            for file in folder[2]:
                if '%s_desc-preproc_bold.nii.gz'%(self.space) in file:
                    os.system('fslmaths %s -mas %s %s'%(os.path.join(self.subj.prep_dir,folder[0],file), self.subj.refvol_mask, os.path.join(self.subj.func,file)))

    def bbreg(self):
        # reg = 'export SUBJECTS_DIR=%s; bbregister --s %s --mov %s --bold --reg %s'%(fs_dir,self.subj.fsub,self.subj.refvol_brain,self.subj.fs_regmat)
        # v2v = 'export SUBJECTS_DIR=%s; mri_vol2vol --subject %s --targ %s --mov %s --reg %s --nearest --inv --o %s'%(fs_dir,self.subj.fsub,os.path.join(self.subj.fs_dir,'mri','aparc+aseg.mgz'),self.subj.refvol_brain,self.subj.fs_regmat,self.subj.faa)
        # mask = 'fslmaths %s -mas %s %s'%(self.subj.faa,self.subj.refvol_mask,self.subj.faa)
        # for cmd in [reg,v2v,mask]: os.system(cmd)

        os.system('fslmaths %s/ses-1/func/%s_ses-1_task-extinction_space-%s_desc-aparcaseg_dseg.nii.gz -mas %s %s'%(self.subj.prep_dir,self.subj.fsub,self.space,self.subj.refvol_mask,self.subj.faa))

    def fs_mask(self):

        rois = {
                 'mOFC':{'l':1014,
                         'r':2014},

                 'amyg':{'l':18,
                         'r':54},

                 'hpc':{'l':17,
                        'r':53},

                 'dACC':{'l':1002,
                         'r':2002},

                 'ins':{'l':1035,
                        'r':2035}
                 }

        for roi in rois:
            for hemi in rois[roi]:
                #make the hemisphere specific masks
                out = os.path.join(self.subj.masks,'%sh_%s.nii.gz'%(hemi,roi))
                thr = rois[roi][hemi]
                cmd = 'fslmaths %s -thr %s -uthr %s -bin %s'%(self.subj.faa,thr,thr,out)
                os.system(cmd)
            #make the joint mask
            l = os.path.join(self.subj.masks,'lh_%s.nii.gz'%(roi))
            r = os.path.join(self.subj.masks,'rh_%s.nii.gz'%(roi))
            out = os.path.join(self.subj.masks,'%s_mask.nii.gz'%(roi))
            cmd = 'fslmaths %s -add %s -bin %s'%(l,r,out)
            os.system(cmd)

    def group_func_mask(self):
        
        rois = {'fvmPFC':os.path.join(SCRATCH,'roi_masks','thr_CSm_CSp.nii.gz'),
                'fdACC':os.path.join(SCRATCH,'roi_masks','thr_CSp_CSm.nii.gz')}
        vals = {'fvmPFC':[1014,2014,1026,2026],
                'fdACC':[1002,2002,1023,2023,1028,2028]}

        for roi in rois:
            out_mask = os.path.join(self.subj.masks,'%s_mask.nii.gz'%(roi))
            
            os.system('flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour'%(
                    rois[roi], self.subj.refvol_brain, self.subj.std2ref, out_mask))
            
            os.system('fslmaths %s -mul %s %s'%(out_mask,self.subj.faa,out_mask)) 

            mask_img = nib.load(out_mask)
            mask_dat = get_data(mask_img)

            coor = [np.where(mask_dat == val) for val in vals[roi]]
            #set all values to 0 except the ones we want
            mask_dat[:,:,:] = 0
            for val in coor: mask_dat[val] = 1

            #make it a nifti
            mask_img = new_img_like(mask_img,mask_dat,affine=mask_img.affine,copy_header=True)
            #and save
            nib.save(mask_img,os.path.join(self.subj.masks,'%s_mask.nii.gz'%(roi)))


    def fsl_reg(self):

        os.system('flirt -in %s -ref %s -dof 12 -omat %s'%(self.subj.refvol_brain, std_2009_brain, self.subj.ref2std))
        os.system('convert_xfm -omat %s -inverse %s'%(self.subj.std2ref,self.subj.ref2std))


        # standard      = '/work/IRC/ls5/opt/apps/fsl-5.0.10/data/standard/MNI152_T1_1mm_brain'
        # standard_head = '/work/IRC/ls5/opt/apps/fsl-5.0.10/data/standard/MNI152_T1_1mm'
        # standard_mask = '/work/IRC/ls5/opt/apps/fsl-5.0.10/data/standard/MNI152_T1_1mm_brain_mask_dil'
        
        # os.system('epi_reg --epi=%s --t1=%s --t1brain=%s --out=%s'%(self.subj.refvol_brain,self.subj.t1,self.subj.t1_brain,self.subj.ref2t1[:-4]))

        # os.system('flirt -in %s -ref %s -omat %s \
        #                 -cost corratio -dof 12 -searchrx -90 90 \
        #                 -searchry -90 90 -searchrz -90 90 -interp trilinear'%(self.subj.t1_brain,standard,self.subj.t12std))

        # os.system('fnirt --in=%s --aff=%s --cout=%s \
        #                  --config=T1_2_MNI152_2mm \
        #                  --ref=%s --refmask=%s --warpres=10,10,10'%(self.subj.t1,self.subj.t12std,self.subj.t12std_warp,standard_head,standard_mask))

        # os.system('applywarp --ref=%s --in=%s \
        #                      --out=%s --warp=%s \
        #                      --premat=%s --interp=nn'%(standard,self.subj.faa,self.subj.saa,self.subj.t12std_warp,self.subj.ref2t1))

def std_space_check():
    from nilearn.image import get_data
    pdir = 'D:\\fc-bids\\derivatives\\fmriprep'
    for sub in all_sub_args:
        print(sub)
        subj = bids_meta(sub)
        sub_dir = os.path.join(pdir,subj.fsub)
        aparcs = glob('%s/ses-*/func/%s_ses-*_task-**_space-MNI152NLin2009cAsym_desc-aparcaseg_dseg.nii.gz'%(sub_dir,subj.fsub))
        for a in aparcs:
            assert np.array_equal(get_data(aparcs[0]),get_data(a))

def group_std_masks():

    saa = {}
    
    thr = {'mOFC':     [1014,2014]
            'dACC':     [1002,2002],
            'lh_amyg':  [18],
            'rh_amyg':  [54],
            'lh_hpc':   [17],
            'rh_hpc':   [53]}

    masks = {}

    for sub in all_sub_args:
        subj = bids_meta(sub)

        os.system('flirt -in %s -out %s -ref %s -applyxfm -init %s -interp=nn'%(subj.faa,subj.saa,std_2009_brain,subj.ref2std))

        saa[sub] = get_data(subj.saa)

        


