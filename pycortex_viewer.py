'''import MNI brain into pycortex'''
import cortex
import numpy as np
from nilearn import image
import matplotlib
from matplotlib import cm
import os
import nibabel as nib

for group in ['healthy','ptsd']:
    for phase in ['acquisition','extinction']:
        q = image.get_data(f'/mnt/c/Users/ACH/Desktop/ers_comps/{group}_{phase}/{group}_{phase}_ClusterEffEst.nii.gz')
        print(np.min(q[np.where(q!=0)]),np.max(q))
for group in ['healthy','ptsd']:
    for phase in ['acquisition','extinction']:
        for hemi in ['l','r']:
            os.system(f'mri_vol2surf --src ers_comps/{group}_{phase}/{group}_{phase}_ClusterEffEst.nii.gz \
             --srcreg ers_comps/surf/anat_3mm_to_1mm.dat --out ers_comps/surf/{hemi}h_{group}_{phase}.mgh \
             --hemi {hemi}h --projfrac 1 --fwhm 0')

acq_norm = matplotlib.colors.Normalize(vmin=.087,vmax=.59)
ext_norm = matplotlib.colors.Normalize(vmin=.1,vmax=.3)

R = cm.get_cmap('Reds')
R.set_bad(alpha=0)
R.set_under(alpha=0)

B = cm.get_cmap('Blues')
B.set_under(alpha=0)
B.set_bad(alpha=0)



data = {group:{phase:np.concatenate((np.array(image.get_data(f'ers_comps/surf/lh_{group}_{phase}.mgh').ravel()),np.array(image.get_data(f'ers_comps/surf/rh_{group}_{phase}.mgh').ravel()))) for phase in ['acquisition','extinction']} for group in ['healthy','ptsd']}

data['healthy']['acquisition'] = R(acq_norm(data['healthy']['acquisition']))
data['healthy']['extinction'] = B(ext_norm(data['healthy']['extinction']))

data['ptsd']['acquisition'] = R(acq_norm(data['ptsd']['acquisition']))
data['ptsd']['extinction'] = B(ext_norm(data['ptsd']['extinction']))

verts = []
for group in ['healthy','ptsd']:
    for phase in ['acquisition','extinction']:
        verts.append(cortex.VertexRGB(data[group][phase][:,0],data[group][phase][:,1],data[group][phase][:,2],subject='MNI2009c',alpha=data[group][phase][:,3],description=f'{group}_{phase}')        )

verts = {str(i):verts[i] for i in range(4)}

cortex.webshow(verts)






cortex.freesurfer.import_subj('MNI2009c',freesurfer_subject_dir='/mnt/c/Users/ACH/Desktop',)

cortex.xfm.Transform.from_freesurfer('/mnt/c/Users/ACH/Desktop/MNI2009c/mri/transforms/talairach.xfm',
                                    '/mnt/c/Users/ACH/Desktop/standard/MNI152NLin2009cAsym_T1_1mm_brain.nii.gz',
                                    'MNI2009c',
                                    freesurfer_subject_dir='/mnt/c/Users/ACH/Desktop',
        ).save(subject='MNI2009c',
               name='tal', # places into subj/transforms/xfm_name
               xfmtype='coord')



ref = nib.load('standard/MNI152NLin2009cAsym_T1_1mm_brain.nii.gz')
mask = 'standard/combo_mask.nii.gz'
new_mask = image.new_img_like(ref,image.get_data(mask),copy_header=True)
nib.save(new_mask,'standard/good_combo_mask.nii.gz')

for hemi in ['l','r']:
    mask = f'mri_vol2surf --src standard/test_1vox.nii.gz \
             --srcreg ers_comps/surf/anat_1mm_to_1mm.dat \
             --out ers_comps/surf/{hemi}h_combo_mask.mgh \
             --hemi {hemi}h --projfrac 1 --fwhm 0'
    os.system(mask)

data = np.concatenate((np.array(image.get_data(f'ers_comps/surf/lh_combo_mask.mgh').ravel()),np.array(image.get_data(f'ers_comps/surf/rh_combo_mask.mgh').ravel())))
data[np.where(data==0)] = np.nan

BR = cm.get_cmap('Greens')
BR.set_under(alpha=0)

mask_norm = matplotlib.colors.Normalize(vmin=.1,vmax=1)
data = BR(mask_norm(data))

verts = cortex.VertexRGB(data[:,0],data[:,1],data[:,2],subject='MNI2009c',alpha=data[:,3],description=f'group_mask')
cortex.webshow(verts)