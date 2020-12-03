from fg_config import *
from graphing_functions import *
from roi_rsa import *
from robust_corr import *

pfc_rois = ['rACC','sgACC']
hc_rois = ['hc_head','hc_body','hc_tail']
amyg_rois = ['amyg_cem','amyg_bla']

paper_rois = pfc_rois + hc_rois + amyg_rois

roi_list = [pfc_rois,hc_rois,amyg_rois]

phase3 = ['baseline','acquisition','extinction']

###############################################################

c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)


cdf = c.df.dropna(subset=['response'])
pdf = p.df.dropna(subset=['response'])
cdf = cdf.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = pdf.groupby(['trial_type','encode_phase','roi','subject']).mean()
mdf = pd.concat((cdf,pdf)).reset_index()

# mdf = mdf.set_index('subject').drop([20,120]).reset_index()

# mdf = (mdf.loc['CS+'] - mdf.loc['CS-']).reset_index()
mdf['group'] = mdf.subject.apply(lgroup)

mdf = mdf.set_index('roi')
mdf = mdf.loc[paper_rois].reset_index()


mdf = mdf.set_index(['group','encode_phase','roi','trial_type','subject']).sort_index()

stats = pd.DataFrame(columns=['w','p','cles','p_fdr'],
                         index=pd.MultiIndex.from_product([groups,paper_rois,phase3],
                         names=['group','roi','phase']))

for group in groups:
    for phase in phase3:
        for roi in paper_rois:
            wres = pg.wilcoxon(mdf.loc[(group,phase,roi,'CS+'),'rsa'], mdf.loc[(group,phase,roi,'CS-'),'rsa'], tail='greater')
            stats.loc[(group,roi,phase),['w','p','cles']] = wres[['W-val','p-val','CLES']].values
    for rlist in roi_list:
        stats.loc[(group,rlist),'p_fdr'] = pg.multicomp(list(stats.loc[(group,rlist),'p'].values),method='fdr_bh')[1]

diff = mdf.reset_index().set_index(['trial_type','group','roi','encode_phase','subject'])
diff = (diff.loc['CS+'] - diff.loc['CS-']).reset_index()

diff.roi = diff.roi.apply(pfc_rename).apply(amyg_rename)
stats = stats.reset_index()
stats.roi = stats.roi.apply(pfc_rename).apply(amyg_rename)
stats = stats.set_index(['group','roi','phase'])
diff = diff.set_index(['group','roi','encode_phase']).sort_index()

cscomp('healthy',diff,['dACC','vmPFC'],stats,phases=phase3)
cscomp('ptsd',diff,['dACC','vmPFC'],stats,phases=phase3)

cscomp('healthy',diff,['Amyg. BLA','Amyg. CeM'],stats,phases=phase3)
cscomp('ptsd',diff,['Amyg. BLA','Amyg. CeM'],stats,phases=phase3)

pg.wilcoxon(diff.loc[('healthy','dACC','extinction'),'rsa'], diff.loc[('healthy','dACC','acquisition'),'rsa'])
pg.wilcoxon(diff.loc[('healthy','vmPFC','extinction'),'rsa'], diff.loc[('healthy','vmPFC','acquisition'),'rsa'])
pg.wilcoxon(diff.loc[('healthy','amyg_cem','extinction'),'rsa'], diff.loc[('healthy','amyg_cem','acquisition'),'rsa'])

pg.wilcoxon(diff.loc[('ptsd','dACC','extinction'),'rsa'], diff.loc[('ptsd','dACC','acquisition'),'rsa'])
pg.wilcoxon(diff.loc[('ptsd','vmPFC','extinction'),'rsa'], diff.loc[('ptsd','vmPFC','acquisition'),'rsa'])

pg.mwu(diff.loc[('healthy','vmPFC','extinction'),'rsa'], diff.loc[('ptsd','vmPFC','extinction'),'rsa'])
pg.mwu(diff.loc[('healthy','dACC','extinction'),'rsa'], diff.loc[('ptsd','dACC','extinction'),'rsa'])






################################some very quick connectivity stuff
seeds = ['hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']
# seeds = ['rh_hc_tail','lh_hc_tail','rh_hc_body','lh_hc_body','rh_hc_head','lh_hc_head','rh_amyg_bla','lh_amyg_bla','rh_amyg_cem','lh_amyg_cem']
targets = ['rACC','sgACC','hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']
copes = {'ext_acq':     ['CS+E','CS+A'],
         'ext_csp_csm': ['CS+E','CS-E'],
         'acq_csp_csm': ['CS+A','CS-A'],
         'ext_bsl':     ['CS+E','CS+B'],
         'acq_bsl':     ['CS+A','CS+B']}

df = pd.read_csv('extracted_mem_apriori_gPPI.csv')
# df = pd.read_csv('extracted_hemi_mem_gPPI.csv'
    # ).set_index('subject').drop([20,120]).reset_index()

df['group'] = df.subject.apply(lgroup)
df = df.set_index(['cope','group','seed','target','subject']).sort_index()

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
for col in stats.columns: stats[col] = stats[col].apply(lambda x: x[0] if type(x) == list else x)
stats.p_mask = stats.p_fdr.apply(lambda x: 0 if x >.05 else 1)
stats['cles_disp'] = stats.cles * stats.p_mask


roi = 'amyg_bla'
cem = df.reset_index()
cem = cem[cem.seed == roi]
cem = pd.concat((cem[cem.cope == 'CS+E'], cem[cem.cope == 'CS+A'])).set_index('cope')

cem = (df.loc['CS+E'] - df.loc['CS+A'])

pg.mwu(cem.loc[('ptsd',roi,'rACC'),'conn'], cem.loc[('healthy',roi,'rACC'),'conn'])
pg.wilcoxon(cem.loc[('ptsd',roi,'rACC'),'conn'], cem.loc[('ptsd',roi,'sgACC'),'conn'])
pg.wilcoxon(cem.loc[('healthy',roi,'rACC'),'conn'], cem.loc[('healthy',roi,'sgACC'),'conn'])


cem = cem.reset_index()
cem = cem[cem.seed == roi]

cem.target = cem.target.apply(pfc_rename)

qpal = sns.xkcd_palette(['windows blue','amber'])
fig, ax = plt.subplots(figsize=(10,8))
sns.barplot(data=cem, x='group', y='conn', hue='target',palette=qpal,ax=ax)
ax.set_xticklabels(['Healthy','PTSD'],fontsize=20)
ax.set_xlabel('')
ax.set_ylabel('∆ Connectivity')
legend_elements = [Patch(facecolor=qpal[0],edgecolor=None,label='dACC'),
                   Patch(facecolor=qpal[1],edgecolor=None,label='vmPFC')]
ax.legend(handles=legend_elements,loc='lower center')
ax.legend_.set_title('gPPI Target')
ax.set_title('gPPI seed = Amygdala CeM',fontsize=30)
plt.tight_layout()
# #lets look at this split out
# df = pd.read_csv('extracted_mem_apriori_gPPI.csv')
# df['group'] = df.subject.apply(lgroup)
# df['phase'] = df.cope.apply(lambda x: x[-1])
# df['condition'] = df.cope.apply(lambda x: x[:-1])
# df = df[df.phase != 'B']

# for seed in df.seed.unique():
#     sns.catplot(data=df[df.seed == 'amyg_cem'], x='condition', y='conn', hue='phase', palette='mako', hue_order=['A','E'], col='target', row='group', kind='bar')
#     plt.suptitle('amyg_cem')

'''subcortical univariate to crotical ERS'''
df = pd.read_csv('ers_subcort_betas_diff.csv').set_index(['group','phase'])#,'subject'])
df = df.reset_index().set_index(['group','phase','subject'])
legend_elements = [Patch(facecolor=gpal[0],edgecolor=None,label='Healthy'),
                   Patch(facecolor=gpal[1],edgecolor=None,label='PTSD')]

phase = 'extinction'
for seed in ['amyg_bla','hc_tail']:
    fig, ax = plt.subplots()
    skipped_corr(df.loc[('healthy',phase),seed],df.loc[('healthy',phase),'dACC_ers'],vis=True,ax=ax,color=gpal[0])
    skipped_corr(df.loc[('ptsd',phase),seed],df.loc[('ptsd',phase),'dACC_ers'],vis=True,ax=ax,color=gpal[1])
    ax.set_ylabel('dACC Extinction E-R Overlap (CS+ > CS-)')
    ax.set_xlabel(f'{seed} activity (CS+ > CS-)')
    ax.set_title(f'{seed} activity predicts dACC Extinction E-R overlap in PTSD')
    ax.legend(handles=legend_elements)
    plt.tight_layout()

fig, ax = plt.subplots()
skipped_corr(df.loc[('healthy',phase),'hc_head'],df.loc[('healthy',phase),'vmPFC_ers'],vis=True,ax=ax,color=gpal[0])
skipped_corr(df.loc[('ptsd',phase),'hc_head'],df.loc[('ptsd',phase),'vmPFC_ers'],vis=True,ax=ax,color=gpal[1])
ax.set_ylabel('vmPFC Extinction E-R Overlap (CS+ > CS-)')
ax.set_xlabel(f'Anterior HPC activity (CS+ > CS-)')
ax.set_title(f'Anterior HPC activity predicts vmPFC Extinction E-R overlap in Healthy')
ax.legend(handles=legend_elements)
plt.tight_layout()