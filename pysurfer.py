import os
import numpy as np
import nibabel as nib
from surfer import Brain
import matplotlib
matplotlib.use('TKAgg')

def bnsurf(data=None,val=None,cmap=None,min_val=None,max_val=None,mid_val=None,tail='greater',out=None):

    subjects_dir = '/Users/ACH/Documents/freesurfer/subjects'
    subject_id = 'fsaverage'
    hemi = 'lh'
    surf = 'inflated'

    brain = Brain(subject_id, hemi, surf, size=(1200, 1200), background='w', subjects_dir=subjects_dir,
                  interaction='terrain', cortex='low_contrast', units='mm', title=out)
    
    aparc_file = os.path.join(subjects_dir,
                              subject_id, "label",
                              hemi + 'BN_pfc.annot')
                              # hemi + ".BN_Atlas.annot")

    labels, ctab, names = nib.freesurfer.read_annot(aparc_file)
    # [187,179,183,177,41,47,49,13,11,1,9]
    # pfc_names = [b'Unknown',b'A32sg_L',b'A32p_L',b'A24cd_L',b'A24rv_L',b'A14m_L',b'A11m_L',b'A13_L',b'A10m_L',b'A9m_L',b'A8m_L',b'A6m_L']
    # roi_data = rs.uniform(.5, .8, size=len(names))
    rois = {
             'A32sg':{'name':b'A32sg_L','label':187},
             'A32p':{'name':b'A32p_L','label':179},
             'A24cd':{'name':b'A24cd_L','label':183},
             'A24rv':{'name':b'A24rv_L','label':177},
             'A14m':{'name':b'A14m_L','label':41},
             'A11m':{'name':b'A11m_L','label':47},
             'A13':{'name':b'A13_L','label':49},
             'A10m':{'name':b'A10m_L','label':13},
             'A9m':{'name':b'A9m_L','label':11},
             'A8m':{'name':b'A8m_L','label':1},
             'A6m':{'name':b'A6m_L','label':9}
             }
    #initialize the data with -1 so ROIs we don't care about are transparent
    if tail == 'greater':
        roi_data = np.repeat([-1.],len(names))
    elif tail == 'less':
        roi_data = np.ones(len(names))
    elif tail == 'two':
        roi_data = np.zeros(len(names))


    #this is where we apply the value we want
    for roi in rois:
        roi_data[rois[roi]['label']] = data.loc[roi,val]

    #this maps the roi values onto the vertex data using repeated sampling
    vtx_data = roi_data[labels]

    #again cancelling out ROIs we don't care about
    if tail == 'greater':
        vtx_data[labels == -1] = -1
    if tail == 'less':
        vtx_data[labels == -1] = -100
    #this is where we actually put it on the brain
    if tail == 'greater':
        #get min/max values if not specified
        if min_val is None: min_val = data[val].astype(float).drop_duplicates().nsmallest(2)[1]
        if max_val is None: max_val = data[val].max()
    
        brain.add_data(vtx_data, min=min_val, max=max_val, mid=mid_val, colormap=cmap, thresh=0.00001, alpha=1)
    
    elif tail == 'less':
        brain.add_data(vtx_data, data[val].min(), 0, colormap=cmap, thresh=data[val].min(), alpha=1)
    elif tail == 'two':
        brain.add_data(vtx_data, data[val].min(), data[val].max(), center=0, colormap=cmap,)
    # brain.show_view('m')
    brain.show_view({'azimuth':40,'elevation':100})

    if out is not None:
        brain.save_image('brainnetome_maps/%s.png'%(out),antialiased=True)