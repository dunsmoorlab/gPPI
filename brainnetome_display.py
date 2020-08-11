from fg_config import *
import nibabel as nib
from nilearn.image import get_data, new_img_like
from roi_rsa import *
from graphing_functions import *
from pysurfer import bnsurf
from scipy.stats import wilcoxon

rois = ['A32sg','A32p','A24cd','A24rv','A14m','A11m','A13','A10m','A9m','A8m','A6m']


c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)

groups = ['healthy','ptsd']
memcon = ['encoding','retrieval']
levels = ['item','set']
mems = ['hit','miss']
cons = ['CS+','CS-']
phases = {'baseline':24,'acquisition':24,'early_extinction':8,'extinction':16}
subs = range(24)
conds = ['CSp','CSm']

cdf = c.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = p.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
mdf = pd.concat((cdf,pdf)).reset_index()

mdf = mdf.set_index('subject').drop([20,120]).reset_index()

# mdf = (mdf.loc['CS+'] - mdf.loc['CS-']).reset_index()
mdf['group'] = mdf.subject.apply(lgroup)
mdf['level'] = 'item'
phase3 = ['baseline','acquisition','extinction']
# mdf = mdf.set_index(['group','roi','encode_phase']).sort_index()
mdf = mdf.set_index(['group','roi','encode_phase','trial_type','subject']).sort_index()


stats = pd.DataFrame(columns=['w','p','cles'],index=pd.MultiIndex.from_product([['healthy','ptsd'],phase3,rois],names=['group','encode_phase','roi']))

for group in ['healthy','ptsd']:
    for phase in phase3:
        for roi in rois:
            # tres = pg.ttest(mdf.loc[(group,roi,phase),'rsa'],0,tail='greater').values
            wres = pg.wilcoxon(mdf.loc[(group,roi,phase,'CS+'),'rsa'],mdf.loc[(group,roi,phase,'CS-'),'rsa'],tail='greater').values
            stats.loc[(group,phase,roi),'w']    = wres[0,0]
            stats.loc[(group,phase,roi),'p']    = wres[0,2]
            stats.loc[(group,phase,roi),'cles'] = wres[0,-1]
        stats.loc[(group,phase),'p'] = pg.multicomp(list(stats.loc[(group,phase),'p'].values),method='fdr_bh')[1]
# stats.p = 1 - stats.ps
# stats.p = stats.p.apply(lambda x: 0 if x < .95 else x)
stats['p_mask'] = stats.p.apply(lambda x: 0 if x >.05 else 1)
stats['cles_disp'] = stats.cles * stats.p_mask


cmaps = {'baseline':'Greys',
         'acquisition':'Purples',
         'early_extinction':sns.light_palette('seagreen',as_cmap=True),
         'extinction':'Greens'}
# dummy = np.zeros(atlas.shape)
for group in ['healthy','ptsd']:
    for phase in phase3:
        disp = stats.loc[group,phase].copy()
        if disp.cles_disp.max() == 0: 
            pass
        else:
            bnsurf(data=disp,val='cles_disp',min_val=.5,max_val=stats.cles_disp.max(),cmap=cmaps[phase],out='test/%s_%s'%(group,phase))

#############memory effect
mem = 'high_confidence_accuracy'
cm = c.df.groupby(['trial_type',mem,'encode_phase','roi','subject']).mean()
pm = p.df.groupby(['trial_type',mem,'encode_phase','roi','subject']).mean()
memdf = pd.concat((cm,pm)).reset_index()
memdf = memdf.set_index('subject')
memdf = memdf.drop(index=[18,20,120]).reset_index()
memdf = memdf.set_index(['trial_type',mem,'encode_phase','roi','subject'])
memdf = (memdf.loc['CS+'] - memdf.loc['CS-']).reset_index()

# memdf = (memdf.loc[1] - memdf.loc[0]).reset_index()
memdf['group'] = memdf.subject.apply(lgroup)
if 'high' in mem: memdf[mem] = memdf[mem].apply(lambda x: 'hit' if x == 1 else 'miss')
memdf = memdf.set_index(['group','roi','encode_phase',mem,'subject']).sort_index()


mstats = pd.DataFrame(columns=['t','p'],index=pd.MultiIndex.from_product([['healthy','ptsd'],phases,rois],names=['group','encode_phase','roi']))

for group in ['healthy','ptsd']:
    for phase in phases:
        for roi in rois:
            # tres = pg.ttest(memdf.loc[(group,roi,phase,'hit'),'rsa'],memdf.loc[(group,roi,phase,'miss'),'rsa'],paired=True).values
            tres = pg.ttest(memdf.loc[(group,roi,phase,'hit'),'rsa'],0,tail='greater').values
            mstats.loc[(group,phase,roi),'t'] = tres[0,0]
            mstats.loc[(group,phase,roi),'p'] = tres[0,3]

        # mstats.loc[(group,phase),'p'] = pg.multicomp(list(mstats.loc[(group,phase),'p'].values),method='fdr_bh')[1]
# stats.p = 1 - stats.ps
# stats.p = stats.p.apply(lambda x: 0 if x < .95 else x)
mstats['p_mask'] = mstats.p.apply(lambda x: 0 if x >.05 else 1)
mstats['t_disp'] = mstats.t * mstats.p_mask
cmaps = {'baseline':'Greys',
         'acquisition':'Purples',
         'early_extinction':sns.light_palette('seagreen',as_cmap=True),
         'extinction':'Greens'}
# dummy = np.zeros(atlas.shape)
for group in ['healthy','ptsd']:
    for phase in phases:
        disp = mstats.loc[group,phase].copy()
        if disp.t_disp.min() == 0 and disp.t_disp.max() == 0:
            pass
        else:
            if disp.t_disp.max() > 0:
                tail = 'greater'
            else:
                tail = 'less'
            bnsurf(disp,'t_disp',cmaps[phase],tail=tail,out='rsa/mem/hit_vs_0/%s_%s'%(group,phase))
            # bnsurf(disp,'t_disp',cmaps[phase],tail=tail,out='rsa/mem/hit_vs_miss/%s_%s'%(group,phase))

######################
#rm ANOVAs
cdf = c.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = p.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
mdf = pd.concat((cdf,pdf)).reset_index()
mdf['group'] = mdf.subject.apply(lgroup)
mdf=mdf[mdf.encode_phase != 'baseline'] #drop baseline phase
mdf=mdf[mdf.encode_phase != 'early_extinction'] #drop early_extinction
# phase3 = ['baseline','acquisition','extinction']
mdf = mdf.set_index(['group','roi']).sort_index()

stats = pd.DataFrame(columns=['F','p'],index=pd.MultiIndex.from_product([['healthy','ptsd'],['condition','phase','interaction'],rois],names=['group','effect','roi']))
for group in ['healthy','ptsd']:
    for roi in rois:
        rm_res = pg.rm_anova(data=mdf.loc[group,roi],dv='rsa',within=['trial_type','encode_phase'],subject='subject')
        stats.loc[(group,'condition',roi),['F','p']] = rm_res.loc[0,['F','p-unc']].values
        stats.loc[(group,'phase',roi),['F','p']] = rm_res.loc[1,['F','p-unc']].values
        stats.loc[(group,'interaction',roi),['F','p']] = rm_res.loc[2,['F','p-unc']].values
stats['p_mask'] = stats.p.apply(lambda x: 0 if x >.05 else 1)
stats['F_disp'] = stats.F * stats.p_mask

cmaps = {'condition':'hot',
         'phase':'BuPu',
         'interaction':'cool',
         }
# dummy = np.zeros(atlas.shape)
for group in ['healthy','ptsd']:
    for effect in ['condition','phase','interaction']:
        disp = stats.loc[group,effect].copy()
        if disp.F_disp.max() == 0: 
            pass
        else:
            bnsurf(disp,'F_disp',cmaps[effect],out='rsa/rm_anova/%s_%s'%(group,effect))

########
#conjunction analysis


########
#memory CS+ lmm
lmm_res = pd.read_csv('rsa_lmm_results.csv')
lmm_res = lmm_res.drop(columns=['Unnamed: 0'])
lmm_res = lmm_res.set_index(['group','phase','condition','roi']).sort_index()
lmm_res['p_mask'] = lmm_res.pval.apply(lambda x: 0 if x >.05 else 1)
lmm_res['B_disp'] = lmm_res.beta * lmm_res.p_mask
# lmm_res.B_disp = lmm_res.B_disp.astype(object)
cmaps = {'baseline':'Greys',
         'acquisition':'Purples',
         'early_extinction':sns.light_palette('seagreen',as_cmap=True),
         'extinction':'Greens'}
# dummy = np.zeros(atlas.shape)
for group in ['healthy','ptsd']:
    for phase in phases:
        disp = lmm_res.loc[group,phase,'CS+'].copy()
        if disp.B_disp.min() == 0 and disp.B_disp.max() == 0:
            pass
        else:
            if disp.B_disp.max() > 0:
                tail = 'greater'
            else:
                tail = 'less'
            bnsurf(disp,'B_disp',cmaps[phase],tail=tail,out='rsa/lmm_mem/%s_%s'%(group,phase))
            # bnsurf(disp,'t_disp',cmaps[phase],tail=tail,out='rsa/mem/hit_vs_miss/%s_%s'%(group,phase))





'''CONNECTIVITY'''
seeds = ['hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']
targets = bn_rois
copes = ['ext_acq','ext_csp_csm','acq_csp_csm']


df = pd.read_csv('extracted_mem_gPPI.csv'
    ).set_index('subject').drop([20,120]).reset_index()
df['group'] = df.subject.apply(lgroup)
df = df.set_index(['cope','group','seed','target','subject'])

stats = pd.DataFrame(columns={'w':0.0,'p':0.0,'p_fdr':0.0,'p_mask':0.0,'cles_disp':0.0},
                     index=pd.MultiIndex.from_product([groups,seeds,copes,targets],
                     names=['group','seed','cope','target'])
                     ).sort_index()
for group in groups:
    for seed in seeds:
        for cope in copes:
            for target in targets:
                wres = wilcoxon(df.loc[(cope,group,seed,target),'conn'].values,correction=True)
                stats.loc[(group,seed,cope,target),['w','p']] = wres[0], wres[1]
            # stats.p = stats.p.astype(float)
            stats.loc[(group,seed,cope),'p_fdr'] = pg.multicomp(list(stats.loc[(group,seed,cope),'p'].values),method='fdr_bh')[1]

stats.p_mask = stats.p_fdr.apply(lambda x: 0 if x >.05 else 1)
# stats['t_disp'] = stats.t * stats.p_mask

# for group in groups:
#     for seed in seeds:
#         for cope in copes:
#             disp = stats.loc[group,seed,cope]
#             if disp.t_disp.min() == 0 and disp.t_disp.max() == 0:
#                 pass
#             else:
#                 if disp.t_disp.max() > 0:
#                     cmap = 'Reds'
#                     tail = 'greater'
#                 else:
#                     cmap = 'Blues_r'
#                     tail = 'less'
#                 bnsurf(disp,'t_disp',cmap,tail=tail,out='conn/%s_%s_%s'%(group,seed,cope))
# #bnsurf(data,val,cmap,tail='greater',out=None):

# stats = stats.reset_index()
# stats.loc[np.where(stats.p < 0.05)[0]]

'''Relating connectivity to RSA'''
conn = pd.read_csv('extracted_mem_gPPI.csv'
    ).set_index('subject').drop([20,120]).reset_index()
conn['group'] = conn.subject.apply(lgroup)
conn = conn.set_index(['group','seed','cope','target','subject']).sort_index()

c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)

cdf = c.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = p.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
mdf = pd.concat((cdf,pdf))
mdf = (mdf.loc['CS+'] - mdf.loc['CS-']).reset_index()

mdf = mdf.set_index('subject').drop([20,120]).reset_index()
mdf = mdf.drop(columns=['low_confidence_accuracy','high_confidence_accuracy','CSp_trial','CSm_trial'])
mdf['group'] = mdf.subject.apply(lgroup)

mdf = mdf.set_index(['group','roi','encode_phase','subject']).sort_index()

stats = stats = pd.DataFrame(columns={'r':0.0,'p':0.0,'p_fdr':0.0,'p_mask':0.0,'cles_disp':0.0},
                     index=pd.MultiIndex.from_product([groups,seeds,copes,phase3,targets],
                     names=['group','seed','cope','rsa_phase','target'])
                     ).sort_index()
for group in groups:
    for seed in seeds:
        for cope in copes:
            for phase in phase3:
                for target in bn_rois:
                    rres = pg.corr(conn.loc[(group,seed,cope,target),'conn'], mdf.loc[(group,target,phase),'rsa'])
                    stats.loc[(group,seed,cope,phase,target),['r','p']] = rres[['r','p-val']].values
                stats.loc[(group,seed,cope,phase),'p_fdr'] = pg.multicomp(list(stats.loc[(group,seed,cope,phase),'p'].values),method='fdr_bh')[1]

stats.p_mask = stats.p_fdr.apply(lambda x: 0 if x >.05 else 1)
stats = stats.reset_index().set_index('rsa_phase').drop('baseline').reset_index().set_index(['group','seed','cope','target'])
stats[stats.p < .05]