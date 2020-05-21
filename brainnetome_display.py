from fg_config import *
import nibabel as nib
from nilearn.image import get_data, new_img_like
from roi_rsa import *
from graphing_functions import *

atlas_str = os.path.join('../../../Desktop/BNA-maxprob-thr0-1mm.nii.gz')
atlas_img = nib.load(atlas_str)
atlas = get_data(atlas_str)

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

c = group_roi_rsa(group='control',ext_split=True,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=True,fs=True,hemi=False)

groups = ['control','ptsd']
memcon = ['encoding','retrieval']
levels = ['item','set']
mems = ['hit','miss']
cons = ['CS+','CS-']
phases = {'baseline':24,'acquisition':24,'early_extinction':8,'extinction':16}
subs = range(24)
conds = ['CSp','CSm']

cdf = c.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = p.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
mdf = pd.concat((cdf,pdf))
mdf = (mdf.loc['CS+'] - mdf.loc['CS-']).reset_index()
mdf['group'] = mdf.subject.apply(lgroup)
mdf['level'] = 'item'
phase3 = ['baseline','acquisition','extinction']
mdf = mdf.set_index(['group','roi','encode_phase']).sort_index()

stats = pd.DataFrame(columns=['t','p'],index=pd.MultiIndex.from_product([['healthy','ptsd'],phases,rois],names=['group','encode_phase','roi']))

for group in ['healthy','ptsd']:
    for phase in phases:
        for roi in rois:
            tres = pg.ttest(mdf.loc[(group,roi,phase),'rsa'],0,tail='greater').values
            stats.loc[(group,phase,roi),'t'] = tres[0,0]
            stats.loc[(group,phase,roi),'p'] = tres[0,3]

        stats.loc[(group,phase),'p'] = pg.multicomp(list(stats.loc[(group,phase),'p'].values),method='fdr_bh')[1]
stats.p = 1 - stats.p
stats.p = stats.p.apply(lambda x: 0 if x < .95 else x)
stats['p_mask'] = stats.p.apply(lambda x: 0 if x <.95 else 1)
stats['t_disp'] = stats.t * stats.p_mask


dummy = np.zeros(atlas.shape)
for group in ['healthy','ptsd']:
    for phase in phases:
        display = np.zeros(atlas.shape)
        for roi in rois:
            lthr = rois[roi][0]
            uthr = rois[roi][1]
            display[np.where((atlas == lthr) | (atlas == uthr))] = stats.loc[(group,phase,roi),'t_disp']
            dummy[np.where((atlas == lthr) | (atlas == uthr))] = .01
        
        display = new_img_like(atlas_img,display)
        nib.save(display,os.path.join('brainnetome_maps/%s_%s.nii.gz'%(group,phase)))

dummy = new_img_like(atlas_img,dummy)
nib.save(dummy,os.path.join('brainnetome_maps/bg_img.nii.gz'))
