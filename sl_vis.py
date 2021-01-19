import os
from fg_config import HOME
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from surfer import Brain, project_volume_data

cmaps = ['Purples','Greens']
reg_file = os.path.join(os.environ["FREESURFER_HOME"],
                    "average/mni152.register.dat")

#display ERS contrasts
# brain.scale_data_colormap(0, maxval/2 ,maxval, transparent=True, center=None)

def vis_ers_comp(group=None,phase=None,surf=None,cmap=None,split='split'):
    mri_file = f'{HOME}/Desktop/ers_comps/{group}_{phase}/{group}_{phase}_ClusterEffEst.nii.gz'
    surf_data_lh = project_volume_data(mri_file, "lh", reg_file, projarg=[0, 1, .5], smooth_fwhm=1)
    surf_data_rh = project_volume_data(mri_file, "rh", reg_file, projarg=[0, 1, .5], smooth_fwhm=1)
    minval = np.min([np.min(surf_data_lh[np.nonzero(surf_data_lh)]),np.min(surf_data_rh[np.nonzero(surf_data_rh)])])
    
    for view in ['med','lat']:
        brain = Brain('fsaverage', split, surf, cortex='low_contrast',size=(400,400),
                        views=view, background='white', foreground=None)
        
        brain.add_data(surf_data_lh, 0, .5, center=None, hemi='lh', thresh=minval, colorbar=False, colormap=cmap)
        brain.add_data(surf_data_rh, 0, .5, center=None, hemi='rh',thresh=minval, colorbar=False, colormap=cmap)
        
        dACC_coords = (1, 21, 27)
        brain.add_foci(dACC_coords, map_surface='pial', hemi='rh',color='orange')
        
        vmPFC_coords = (-4,34,-6)
        brain.add_foci(vmPFC_coords, map_surface='pial', hemi='lh',color='orange')
        input('press a key dipshit')

vis_ers_comp(group='healthy',phase='acquisition',surf='inflated_pre',cmap='Purples')
vis_ers_comp(group='healthy',phase='extinction',surf='inflated_pre',cmap='Greens')

vis_ers_comp(group='ptsd',phase='acquisition',surf='inflated_pre',cmap='Purples')
vis_ers_comp(group='ptsd',phase='extinction',surf='inflated_pre',cmap='Greens')

vis_ers_comp(group='ptsd',phase='acquisition',surf='inflated_pre',cmap='Purples',split='both')



    # cb2 = mpl.colorbar.ColorbarBase(ax[1], cmap=mpl.cm.Greens,
    #                                 norm=norm, orientation='horizontal')
    # cb1.set_label('CS+ - CS- reinstatement')
    # fig.show()


"""
This overlay represents resting-state correlations with a
seed in left angular gyrus. Let's plot that seed.
"""
seed_coords = (-45, -67, 36)
brain.add_foci(seed_coords, map_surface="white", hemi='lh')


'''
a priori rois for paper
'''
mri_file = f'{HOME}/Desktop/standard/vmPFC_mask.nii.gz'
surf_data_lh = project_volume_data(mri_file, "lh", reg_file, projarg=[0, 1, .5], smooth_fwhm=1)
surf_data_rh = project_volume_data(mri_file, "rh", reg_file, projarg=[0, 1, .5], smooth_fwhm=1)
minval = np.min([np.min(surf_data_lh[np.nonzero(surf_data_lh)]),np.min(surf_data_rh[np.nonzero(surf_data_rh)])])
    
for view in ['med','lat']:
    brain = Brain('fsaverage', 'split', surf, cortex='low_contrast',size=(400,400),
                    views=view, background='white', foreground=None)
    
    brain.add_data(surf_data_lh, 0, .5, center=None, hemi='lh', thresh=minval, colorbar=False, colormap=cmap)
    brain.add_data(surf_data_rh, 0, .5, center=None, hemi='rh',thresh=minval, colorbar=False, colormap=cmap)
        
    dACC_coords = (1, 21, 27)
    brain.add_foci(dACC_coords, map_surface='pial', hemi='rh',color='orange')