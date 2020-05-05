from fg_config import *
import nibabel as nib
from nilearn.image import get_data, new_img_like

atlas_str = os.path.join('../../../Desktop/BNA-maxprob-thr0-1mm.nii.gz')
atlas_img = nib.load(atlas_str)
atlas = get_data(atlas_str)

display = np.zeros(atlas.shape)

rois = {
         'A32sg':[187,188],
         'A32p':[179,180],
         'A24cd':[183,184],
         'A24rv':[177,178],
         'A14m':[41,42],
         'A11m':[47,48],
         'A13':[49,50],
         'A10m':[13,14],
         'A9m':[11,12],
         'A8m':[1,2],
         'A6m':[9,10]
         }

for roi in rois:
    lthr = rois[roi][0]
    uthr = rois[roi][1]
    display[np.where((atlas == lthr) | (atlas == uthr))] = np.random.randint(100)

display = new_img_like(atlas_img,display)
nib.save(display,os.path.join('../../../Desktop/display_test.nii.gz'))