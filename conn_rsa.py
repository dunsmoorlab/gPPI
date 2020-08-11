from fg_config import *
import nibabel as nib
from nilearn.image import get_data, new_img_like
from roi_rsa import *
from graphing_functions import *
from pysurfer import bnsurf

rois = ['A32sg','A32p','A24cd','A24rv','A14m','A11m','A13','A10m','A9m','A8m','A6m']


c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)

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
mdf = mdf.set_index(['group','roi','encode_phase','subject']).sort_index()

conn = pd.read_csv('extracted_mem_gPPI.csv')
conn['group'] = conn.subject.apply(lgroup)
conn = conn.set_index(['group','seed','target','cope','subject'])

for roi in ['sgACC','rACC']:
for seed in ['amyg_cem','amyg_bla','hc_head','hc_body','hc_tail']:
    print(pg.corr(conn.loc[('ptsd',seed,'rACC','ext_csp_csm'),'conn'],mdf.loc[('ptsd','rACC','extinction'),'rsa']))




print(pg.corr(conn.loc[('healthy','amyg_cem','rACC','acq_csp_csm'),'conn'],mdf.loc[('healthy','rACC','acquisition'),'rsa']))
print(pg.corr(conn.loc[('ptsd','hc_tail','sgACC','ext_acq'),'conn'],mdf.loc[('ptsd','sgACC','extinction'),'rsa']))
