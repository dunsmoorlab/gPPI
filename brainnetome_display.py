from fg_config import *
import nibabel as nib
from nilearn.image import get_data, new_img_like
from roi_rsa import *
from graphing_functions import *
from pysurfer import bnsurf
from scipy.stats import wilcoxon, pearsonr

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

c.df = c.df.dropna(subset=['response'])
p.df = p.df.dropna(subset=['response'])

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


stats = pd.DataFrame(columns=['w','p','cles'],index=pd.MultiIndex.from_product([['healthy','ptsd'],phase3,bn_rois],names=['group','encode_phase','roi']))

for group in ['healthy','ptsd']:
    for phase in phase3:
        for roi in rois:
            # tres = pg.ttest(mdf.loc[(group,roi,phase),'rsa'],0,tail='greater').values
            wres = pg.wilcoxon(mdf.loc[(group,roi,phase,'CS+'),'rsa'],mdf.loc[(group,roi,phase,'CS-'),'rsa'],tail='greater').values
            stats.loc[(group,phase,roi),'w']    = wres[0,0]
            stats.loc[(group,phase,roi),'p']    = wres[0,2]
            stats.loc[(group,phase,roi),'cles'] = wres[0,-1]
        stats.loc[(group,phase,bn_rois),'p_corr'] = pg.multicomp(list(stats.loc[(group,phase,bn_rois),'p'].values),method='fdr_bh')[1]
# stats.p = 1 - stats.ps
# stats.p = stats.p.apply(lambda x: 0 if x < .95 else x)
stats['p_mask'] = stats.p_corr.apply(lambda x: 0 if x >.05 or pd.isna(x) else 1)
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

##########################basic effect group differences
cdf = c.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = p.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
mdf = pd.concat((cdf,pdf))
mdf = (mdf.loc['CS+'] - mdf.loc['CS-']).reset_index()

mdf = mdf.set_index('subject').drop([20,120]).reset_index()
mdf['group'] = mdf.subject.apply(lgroup)
mdf = mdf.set_index(['group','encode_phase','roi'])

stats = pd.DataFrame(columns=['w','p','cles'],index=pd.MultiIndex.from_product([phase3,bn_rois],names=['encode_phase','roi']))
for phase in phase3:
    for roi in bn_rois:
        wres = pg.mwu(mdf.loc[('healthy',phase,roi),'rsa'], mdf.loc[('ptsd',phase,roi),'rsa'], tail='two-sided')
        stats.loc[(phase,roi),['w','p','cles']] = wres[['U-val','p-val','CLES']].values
    stats.loc[(phase),'p_corr'] = pg.multicomp(list(stats.loc[(phase),'p'].values),method='fdr_bh')[1]
stats.p = stats.p.apply(lambda x: x[0] if type(x) == list else x)
stats.cles = stats.cles.apply(lambda x: x[0] if type(x) == list else x)

stats['p_mask'] = stats.p.apply(lambda x: 0 if x >.05 else 1)
stats['cles_disp'] = stats.cles * stats.p_mask

for phase in phase3:
    disp = stats.loc[phase].copy()
    if disp.cles_disp.max() == 0: 
        pass
    else:
        bnsurf(data=disp,val='cles_disp',min_val=.25,max_val=.75,mid_val=.5,cmap='RdBu',out='test/group_diff_%s'%(phase))


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
# seeds = ['rh_hc_tail','lh_hc_tail','rh_hc_body','lh_hc_body','rh_hc_head','lh_hc_head','rh_amyg_bla','lh_amyg_bla','rh_amyg_cem','lh_amyg_cem']
# targets = bn_rois
targets = ['rACC','sgACC']
copes = {'ext_acq':     ['CS+E','CS+A'],
         'ext_csp_csm': ['CS+E','CS-E'],
         'acq_csp_csm': ['CS+A','CS-A']}


df = pd.read_csv('extracted_mem_apriori_gPPI.csv')
# df = pd.read_csv('extracted_hemi_mem_gPPI.csv'
    # ).set_index('subject').drop([20,120]).reset_index()

df['group'] = df.subject.apply(lgroup)
df = df.set_index(['cope','group','seed','target','subject'])

stats = pd.DataFrame(columns={'w':0.0,'p':0.0,'p_fdr':0.0,'p_mask':0.0,'cles':0.0},
                     index=pd.MultiIndex.from_product([groups,seeds,copes,targets],
                     names=['group','seed','cope','target'])
                     ).sort_index()
for group in groups:
    for seed in seeds:
        for cope in copes:
            con1, con2 = copes[cope][0], copes[cope][1]
            for target in targets:
                wres = pg.wilcoxon(df.loc[(con1,group,seed,target),'conn'],df.loc[(con2,group,seed,target),'conn'])
                stats.loc[(group,seed,cope,target),['w','p','cles']] = wres[['W-val','p-val','CLES']].values

                #these lines for bilateral rois
                # wres = wilcoxon(df.loc[(cope,group,seed,target),'conn'].values,correction=True)
                # stats.loc[(group,seed,cope,target),['w','p']] = wres[0], wres[1]

            stats.loc[(group,seed,cope),'p_fdr'] = pg.multicomp(list(stats.loc[(group,seed,cope),'p'].values),method='fdr_bh')[1]

stats.p_mask = stats.p_fdr.apply(lambda x: 0 if x >.05 else 1)
stats['cles_disp'] = stats.cles * stats.p_mask

cmaps = {'acquisition':'Purples',
         'extinction':'Greens'}
for group in groups:
    for seed in seeds:
        for cope in copes:
            disp = stats.loc[group,seed,cope].copy()
            if disp.cles_disp.max() == 0: 
                pass
            else:
                bnsurf(data=disp,val='cles_disp',min_val=.1,max_val=.5,cmap='Blues_r',out='final_ish/conn_%s_%s_%s'%(group,seed,cope))


###########Group differences in connectivity############
stats = pd.DataFrame(columns={'w':0.0,'p':0.0,'p_fdr':0.0,'p_mask':0.0,'cles':0.0},
                     index=pd.MultiIndex.from_product([seeds,copes,targets],
                     names=['seed','cope','target'])
                     ).sort_index()
for seed in seeds:
    for cope in copes:
        con1, con2 = copes[cope][0], copes[cope][1]
        for target in targets:

            H = df.loc[(con1,'healthy',seed,target),'conn'] - df.loc[(con2,'healthy',seed,target),'conn']
            P = df.loc[(con1,'ptsd',seed,target),'conn'] - df.loc[(con2,'ptsd',seed,target),'conn']

            wres = pg.mwu(H,P)
            stats.loc[(seed,cope,target),['w','p','cles']] = wres[['U-val','p-val','CLES']].values
        
        stats.loc[(seed,cope),'p_fdr'] = pg.multicomp(list(stats.loc[(seed,cope),'p'].values),method='fdr_bh')[1]
stats.p_mask = stats.p_fdr.apply(lambda x: 0 if x >.05 else 1)
stats['cles_disp'] = stats.cles * stats.p_mask

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
# conn = pd.read_csv('extracted_mem_gPPI.csv'
df = pd.read_csv('extracted_mem_apriori_gPPI.csv')
    # ).set_index('subject').drop([20,120]).reset_index()
df['group'] = df.subject.apply(lgroup)
df = df.set_index(['cope','group','seed','target','subject']).sort_index()

conn = (df.loc['CS+E'] - df.loc['CS-E']).rename(columns={'conn':'ext_csp_csm'})
conn['acq_csp_csm'] = (df.loc['CS+A'] - df.loc['CS-E'])
conn['ext_acq'] = (df.loc['CS+E'] - df.loc['CS+A'])
conn['ext_bsl'] = (df.loc['CS+E'] - df.loc['CS+B'])
conn['acq_bsl'] = (df.loc['CS+A'] - df.loc['CS+B'])


conn = conn.reset_index().melt(id_vars=['group','seed','target','subject'],
                        value_vars=['ext_csp_csm', 'acq_csp_csm','ext_acq','ext_bsl','acq_bsl'],
                        value_name='conn',var_name='cope')
conn = conn.set_index(['group','seed','cope','target','subject']).sort_index()


c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)

cdf = c.df.dropna(subset=['response'])
pdf = p.df.dropna(subset=['response'])
cdf = cdf.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = pdf.groupby(['trial_type','encode_phase','roi','subject']).mean()

mdf = pd.concat((cdf,pdf))
mdf = (mdf.loc['CS+'] - mdf.loc['CS-']).reset_index()

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
                for target in targets:
                    # rres = pg.corr(conn.loc[(group,seed,cope,target),'conn'], mdf.loc[(group,target,phase),'rsa'])
                    # stats.loc[(group,seed,cope,phase,target),['r','p']] = rres[['r','p-val']].values

                    rres = pearsonr(conn.loc[(group,seed,cope,target),'conn'], mdf.loc[(group,target,phase),'rsa'])
                    stats.loc[(group,seed,cope,phase,target),['r','p']] = rres

                stats.loc[(group,seed,cope,phase),'p_fdr'] = pg.multicomp(list(stats.loc[(group,seed,cope,phase),'p'].values),method='fdr_bh')[1]

stats.p_mask = stats.p_fdr.apply(lambda x: 0 if x >.05 else 1)
stats = stats.reset_index().set_index('rsa_phase').drop('baseline').reset_index().set_index(['group','seed','cope','target'])
stats[stats.p_mask == 1]






#within correlations



within = pd.read_csv('within_phase2_similarity.csv').set_index(['condition','group','roi','phase','subject'])
within = within[within.memory_phase == 'encoding']
within = within.drop(columns='memory_phase')
within = (within.loc['CS+'] - within.loc['CS-'])

stats = stats = pd.DataFrame(columns={'r':0.0,'p':0.0,'p_fdr':0.0,'p_mask':0.0,'cles_disp':0.0},
                     index=pd.MultiIndex.from_product([groups,phase2,rois],
                     names=['group','phase','roi'])
                     ).sort_index()

for group in groups:
    for phase in phase2:
        for roi in rois:
            rres = pearsonr(within.loc[(group,roi,phase),'rsa'], mdf.loc[(group,roi,phase),'rsa'])
            stats.loc[(group,phase,roi),['r','p']] = rres
        stats.loc[(group,phase,bn_rois),'p_fdr'] = pg.multicomp(list(stats.loc[(group,phase,bn_rois),'p'].values),method='fdr_bh')[1]
stats.p_mask = stats.p_fdr.apply(lambda x: 0 if x >.05 or pd.isna(x) else 1)
# stats = stats.reset_index().set_index('rsa_phase').drop('baseline').reset_index().set_index(['group','seed','cope','target'])
stats[stats.p_mask == 1]
stats['r_disp'] = stats.r * stats.p_mask

#this is just a quick fix bc it doesn't like being passed a negative value
# stats.loc[('healthy','acquisition','hc_head'),'r'] *= -1

cmaps = {'acquisition':'Purples',
         'extinction':'Greens'}
for group in groups:
    for phase in phase2:
        MAX = stats.loc[(slice('healthy','ptsd'),phase),'r_disp'].max()
        disp = stats.loc[group,phase].copy()
        if disp.r_disp.max() == 0: 
            pass
        else:
            bnsurf(data=disp,val='r_disp',min_val=.4,max_val=MAX,cmap=cmaps[phase],out='correlations/within_phase_corr_%s_%s'%(group,phase))






elrsa = pd.read_csv('elrsa_labar.csv').set_index(['el','condition','group','roi','subject'])
elrsa = elrsa[elrsa.memory_phase == 'encoding']
elrsa = elrsa.drop(columns='memory_phase')
elrsa = (elrsa.loc['early'] - elrsa.loc['late'])
elrsa = (elrsa.loc['CS+'] - elrsa.loc['CS-'])

stats = stats = pd.DataFrame(columns={'r':0.0,'p':0.0,'p_fdr':0.0,'p_mask':0.0,'cles_disp':0.0},
                     index=pd.MultiIndex.from_product([groups,phase2,rois],
                     names=['group','phase','roi'])
                     ).sort_index()

for group in groups:
    for phase in phase2:
        for roi in rois:
            rres = pearsonr(elrsa.loc[(group,roi),'rsa'], mdf.loc[(group,roi,phase),'rsa'])
            stats.loc[(group,phase,roi),['r','p']] = rres
        stats.loc[(group,phase,bn_rois),'p_fdr'] = pg.multicomp(list(stats.loc[(group,phase,bn_rois),'p'].values),method='fdr_bh')[1]
stats.p_mask = stats.p_fdr.apply(lambda x: 0 if x >.05 or pd.isna(x) else 1)
# stats = stats.reset_index().set_index('rsa_phase').drop('baseline').reset_index().set_index(['group','seed','cope','target'])
stats[stats.p_mask == 1]



############Extinction to acquisition similarity############

eta = pd.read_csv('ext_to_acq_within_phase.csv').set_index(['condition','group','roi','subject'])
eta = eta[eta.memory_phase == 'retrieval']
eta = eta.drop(columns='memory_phase')
eta = (eta.loc['CS+'] - eta.loc['CS-'])
# eta = eta.loc['CS+']


stats = stats = pd.DataFrame(columns={'r':0.0,'p':0.0,'p_fdr':0.0,'p_mask':0.0,'cles_disp':0.0},
                     index=pd.MultiIndex.from_product([groups,phase2,rois],
                     names=['group','phase','roi'])
                     ).sort_index()
for group in groups:
    for phase in phase2:
        for roi in rois:
            rres = pearsonr(eta.loc[(group,roi),'rsa'], mdf.loc[(group,roi,phase),'rsa'])
            stats.loc[(group,phase,roi),['r','p']] = rres
        stats.loc[(group,phase,bn_rois),'p_fdr'] = pg.multicomp(list(stats.loc[(group,phase,bn_rois),'p'].values),method='fdr_bh')[1]
stats.p_mask = stats.p_fdr.apply(lambda x: 0 if x >.05 or pd.isna(x) else 1)
# stats = stats.reset_index().set_index('rsa_phase').drop('baseline').reset_index().set_index(['group','seed','cope','target'])
stats[stats.p_mask == 1]
stats['r_disp'] = stats.r * stats.p_mask


cmaps = {'acquisition':'Purples',
         'extinction':'Greens'}
for group in groups:
    for phase in phase2:
        MAX = stats.loc[(slice('healthy','ptsd'),phase),'r_disp'].max()
        disp = stats.loc[group,phase].copy()
        if disp.cles_disp.max() == 0: 
            pass
        else:
            bnsurf(data=disp,val='r_disp',min_val=.5,max_val=1,cmap=cmaps[phase],out='correlations/Ext_Acq_corr_%s_%s'%(group,phase))






within = pd.read_csv('within_phase2_similarity.csv').set_index(['condition','phase','group','roi','subject'])
within = within[within.memory_phase == 'encoding']
within = within.drop(columns='memory_phase')
within = (within.loc['CS+'] - within.loc['CS-'])
within = (within.loc['acquisition'] - within.loc['extinction'])

stats = stats = pd.DataFrame(columns={'r':0.0,'p':0.0,'p_fdr':0.0,'p_mask':0.0,'cles_disp':0.0},
                     index=pd.MultiIndex.from_product([groups,phase2,rois],
                     names=['group','phase','roi'])
                     ).sort_index()

for group in groups:
    for phase in phase2:
        for roi in rois:
            rres = pearsonr(within.loc[(group,roi),'rsa'], mdf.loc[(group,roi,phase),'rsa'])
            stats.loc[(group,phase,roi),['r','p']] = rres
        stats.loc[(group,phase,bn_rois),'p_fdr'] = pg.multicomp(list(stats.loc[(group,phase,bn_rois),'p'].values),method='fdr_bh')[1]
stats.p_mask = stats.p_fdr.apply(lambda x: 0 if x >.05 or pd.isna(x) else 1)
# stats = stats.reset_index().set_index('rsa_phase').drop('baseline').reset_index().set_index(['group','seed','cope','target'])
stats[stats.p_mask == 1]