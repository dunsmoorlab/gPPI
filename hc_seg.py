import os
import nibabel as nib
from fg_config import *
from subprocess import Popen

# os.system('source fsdev.sh')
#Popen(['export SUBJECTS_DIR=/scratch/05426/ach3377/fc-bids/derivatives/freesurfer'])

for sub in all_sub_args:
  print(sub)
  subj = bids_meta(sub)

  for hemi in ['r','l']:
    vol = os.path.join(subj.masks,hemi+'h_HCA.nii.gz')
    l2v =['mri_label2vol',
            '--subject', '%s'%(subj.fsub) , 
            '--seg', os.path.join(subj.fs_dir,'mri',hemi+'h.hippoAmygLabels-T1.v21.HBT.FSvoxelSpace.mgz'),
            '--temp', subj.refvol_brain,
            '--reg', subj.fs_regmat,
            '--fillthresh', '.3' ,
            '--o', vol]

    Popen(l2v).wait()

    hc_rois = {'hc_head':232,
               'hc_body':231,
               'hc_tail':226}
        
    for roi in hc_rois:
      thr = str(hc_rois[roi])
      out = os.path.join(subj.masks,hemi+'h_%s.nii.gz'%(roi))
      hc_cmd = ['fslmaths', vol,
                '-thr', thr,
                '-uthr', thr,
                '-bin', out]
      Popen(hc_cmd).wait()

    amyg_rois = {'amyg_bla':[7001,7003],
                 'amyg_cem':[7004,7010]}
    for roi in amyg_rois:
      lthr = str(amyg_rois[roi][0])
      uthr = str(amyg_rois[roi][1])
      out = os.path.join(subj.masks,hemi+'h_%s.nii.gz'%(roi))
      amyg_cmd = ['fslmaths', vol,
                            '-thr', lthr,
                           '-uthr', uthr,
                            '-bin', out ]
      Popen(amyg_cmd).wait()



  for roi in ['hc_head','hc_body','hc_tail','amyg_bla','amyg_cem']:
    cmd = ['fslmaths', os.path.join(subj.masks,'rh_'+roi+'.nii.gz'),
                   '-add', os.path.join(subj.masks,'lh_'+roi+'.nii.gz'),
                   '-bin', os.path.join(subj.masks,roi+'.nii.gz')]
    Popen(cmd).wait()


    #make the full circuit mask
  cmd = ['fslmaths', os.path.join(subj.masks,'hc_head.nii.gz'),
           '-add', os.path.join(subj.masks,'hc_body.nii.gz'),
           '-add', os.path.join(subj.masks,'hc_tail.nii.gz'),
           '-add', os.path.join(subj.masks,'amyg_bla.nii.gz'),
           '-add', os.path.join(subj.masks,'amyg_cem.nii.gz'),
           '-add', os.path.join(subj.masks,'mOFC_mask.nii.gz'),
           '-add', os.path.join(subj.masks,'dACC_mask.nii.gz'),
           '-bin', os.path.join(subj.masks,'circuit_mask.nii.gz')]
  Popen(cmd).wait()
    



