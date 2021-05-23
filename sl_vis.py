import os
import numpy as np
from surfer import Brain, utils, project_volume_data
from mayavi import mlab
mlab.options.offscreen = True
# from mne.viz import Brain
# from mne.transforms import apply_trans

# reg_file = os.path.join(os.environ["FREESURFER_HOME"],
#                     "average/mni152.register.dat")
os.environ['SUBJECTS_DIR'] = '/mnt/c/Users/ACH/Desktop'
reg_file = '/mnt/c/Users/ACH/Desktop/ers_comps/surf/anat_3mm_to_1mm.dat'

def vis_ers_comp(group=None,phase=None,surf=None,cmap=None,split='lh'):
    mri_file = f'/mnt/c/Users/ACH/Desktop/ers_comps/{group}_{phase}/{group}_{phase}_ClusterEffEst.nii.gz' #find the file containing stats
    surf_data_lh = project_volume_data(mri_file, "lh", reg_file, projarg=[0, 1, .01], smooth_fwhm=1) #project to lh

    _max = .55 if phase == 'acquisition' else .3
    for view in ['med','lat']: #lateral and medial views
        brain = Brain('MNI2009c', split, surf, cortex='low_contrast',size=1000,
                        views=view, background='white', foreground=None) #initialize the brain object
        
        brain.add_data(surf_data_lh, 0, _max, center=None, hemi='lh', thresh=None,
             colorbar=False, colormap=cmap, transparent=True) #add lh data
        
        for vert, color in zip([115262,135014],['white','black']): #add focal ROIs
            brain.add_foci(vert,coords_as_verts=True,color=color,alpha=1)

        fname = f'/mnt/c/Users/ACH/Documents/gPPI/paper/pysurfer/{group}_{phase}_{view}.png'
        os.system(f'rm {fname}')
        brain.save_image(fname,antialiased=True)


vis_ers_comp(group='healthy',phase='acquisition',surf='inflated',cmap='Reds')
vis_ers_comp(group='healthy',phase='extinction',surf='inflated',cmap='Blues')

vis_ers_comp(group='ptsd',phase='acquisition',surf='inflated',cmap='Reds')
vis_ers_comp(group='ptsd',phase='extinction',surf='inflated',cmap='Blues')



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