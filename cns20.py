from fg_config import *
from roi_rsa import *
c = group_roi_rsa(group='control',ext_split=True,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=True,fs=True,hemi=False)

groups = ['control','ptsd']
memcon = ['encoding','retrieval']
levels = ['item','set']
mems = ['hit','miss']
cons = ['CS+','CS-']
rois = ['mOFC','dACC','amyg','hpc','ins','hc_head','hc_body','hc_tail','rh_hc_head','rh_hc_body','rh_hc_tail','lh_hc_head','lh_hc_body','lh_hc_tail','amyg_bla','amyg_cem','rh_amyg_bla','rh_amyg_cem','lh_amyg_bla','lh_amyg_cem']
phases = {'baseline':24,'acquisition':24,'early_extinction':8,'extinction':16}
# phases = {'baseline':24,'acquisition':24,'extinction':24}
subs = range(24)

from graphing_functions import *
##############Item level##############################
cdf = c.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = p.df.groupby(['trial_type','encode_phase','roi','subject']).mean()

ers_pvals = pd.DataFrame(index=pd.MultiIndex.from_product([groups,rois,phases,levels]))
for roi in rois:
    for phase in phases:
        ers_pvals.loc[('control',roi,phase,'item'),'W'] = pg.wilcoxon(cdf.loc[('CS+',phase,roi),'rsa'], cdf.loc[('CS-',phase,roi),'rsa'])['p-val'].values
        ers_pvals.loc[('ptsd',roi,phase,'item'),'W'] = pg.wilcoxon(pdf.loc[('CS+',phase,roi),'rsa'], pdf.loc[('CS-',phase,roi),'rsa'])['p-val'].values

mdf = pd.concat((cdf,pdf))
mdf = (mdf.loc['CS+'] - mdf.loc['CS-']).reset_index()
mdf['group'] = mdf.subject.apply(lgroup)
mdf['level'] = 'item'

mdf = mdf.set_index(['group','roi','encode_phase']).sort_index()
# cscomp('control',mdf,['mOFC','dACC'],phases=phases.keys())
# cscomp('ptsd',mdf,['mOFC','dACC'],phases=phases.keys())
# cscomp('control',mdf,['rh_amyg_bla','lh_amyg_bla'],phases=phases.keys())
# cscomp('ptsd',mdf,['rh_amyg_bla','lh_amyg_bla'],phases=phases.keys())


##############Set level##############################
csl = pd.DataFrame(index=pd.MultiIndex.from_product([cons,rois,phases,sub_args],names=['condition','roi','encode_phase','subject']))
psl = pd.DataFrame(index=pd.MultiIndex.from_product([cons,rois,phases,p_sub_args],names=['condition','roi','encode_phase','subject']))
for con in cons:
    for roi in rois:
        for phase in phases:
            for sub in subs:
                csl.loc[(con,roi,phase,sub_args[sub]),'rsa'] = c.mats[roi][sub,slices[con][phase]['encoding'],slices[con][phase]['retrieval']][np.eye(phases[phase]) == 0].mean()
                psl.loc[(con,roi,phase,p_sub_args[sub]),'rsa'] = p.mats[roi][sub,slices[con][phase]['encoding'],slices[con][phase]['retrieval']][np.eye(phases[phase]) == 0].mean()
for roi in rois:
    for phase in phases:
        ers_pvals.loc[('control',roi,phase,'set'),'W'] = pg.wilcoxon(csl.loc[('CS+',roi,phase),'rsa'], csl.loc[('CS-',roi,phase),'rsa'])['p-val'].values
        ers_pvals.loc[('ptsd',roi,phase,'set'),'W'] = pg.wilcoxon(psl.loc[('CS+',roi,phase),'rsa'], psl.loc[('CS-',roi,phase),'rsa'])['p-val'].values

# csl = (csl.loc['CS+'] - csl.loc['CS-']).reset_index()
# psl = (psl.loc['CS+'] - psl.loc['CS-']).reset_index()
sl = pd.concat([csl,psl])
sl = (sl.loc['CS+'] - sl.loc['CS-']).reset_index()
sl['group'] = sl.subject.apply(lgroup)
sl['level'] = 'set'


sl = sl.set_index(['group','roi','encode_phase']).sort_index()
# cscomp('control',sl,['mOFC','dACC'],phases=phases.keys())
# cscomp('ptsd',sl,['mOFC','dACC'],phases=phases.keys())
# cscomp('control',sl,['rh_amyg_bla','lh_amyg_bla'],phases=phases.keys())
# cscomp('ptsd',sl,['rh_amyg_bla','lh_amyg_bla'],phases=phases.keys())

#############item and set level#######################
lev = pd.concat([mdf,sl]).reset_index().set_index(['group','roi','encode_phase','level']).sort_index()
split_level(lev,'control',phases=phases,split='level')
split_level(lev,'ptsd',phases=phases,split='level')
# split_level(lev,'amyg_cem',phases=phases,split='level')
# split_level(lev,'amyg_bla',phases=phases,split='level')
# split_level(lev,'hc_tail',phases=phases,split='level')
# split_level(lev,'hc_body',phases=phases,split='level')
# split_level(lev,'hc_head',phases=phases,split='level')


##############Within session similarity###############
cws = pd.DataFrame(index=pd.MultiIndex.from_product([memcon,cons,rois,phases,sub_args],names=['memory_phase','condition','roi','encode_phase','subject']))
pws = pd.DataFrame(index=pd.MultiIndex.from_product([memcon,cons,rois,phases,p_sub_args],names=['memory_phase','condition','roi','encode_phase','subject']))
for mem in memcon:
    for con in cons:
        for roi in rois:
            for phase in phases:
                for sub in subs:
                    cws.loc[(mem,con,roi,phase,sub_args[sub]),'rsa'] = c.mats[roi][sub,slices[con][phase][mem],slices[con][phase][mem]][np.tril_indices(phases[phase],-1)].mean()
                    pws.loc[(mem,con,roi,phase,p_sub_args[sub]),'rsa'] = p.mats[roi][sub,slices[con][phase][mem],slices[con][phase][mem]][np.tril_indices(phases[phase],-1)].mean()

ws = pd.concat([cws,pws]).reset_index().set_index(['condition','subject','memory_phase','encode_phase','roi'])
ws = (ws.loc['CS+'] - ws.loc['CS-']).reset_index()
ws['group'] = ws.subject.apply(lgroup)
ws = ws.set_index(['group','roi','encode_phase','memory_phase']).sort_index()
#go to graphing_functions and figure out how to swap in memory_phase for levels
split_level(ws,'control',phases=phases,split='memory_phase')
split_level(ws,'ptsd',phases=phases,split='memory_phase')
# split_level(ws,'amyg_cem',phases=phases,split='memory_phase')
# split_level(ws,'amyg_bla',phases=phases,split='memory_phase')
# split_level(ws,'hc_tail',phases=phases,split='memory_phase')
# split_level(ws,'hc_body',phases=phases,split='memory_phase')
# split_level(ws,'hc_head',phases=phases,split='memory_phase')


pal = sns.color_palette(['orangered','slategrey','orange','lightgrey'])

def mem_phase(roi):
    g = sns.catplot(x='encode_phase',y='rsa',hue='gvar',col='group',data=ws.query('roi == @roi'),
                    kind='bar',aspect=1,palette=pal)
    g.set_titles("{col_name} {col_var}")
##############Group comp#################################
#save things out
for name, df in zip(['level','memory'],[lev,ws]):
    cout = df.loc['control']
    cout = cout.loc[['dACC','mOFC']]
    pout = df.loc['ptsd']
    pout = pout.loc[['dACC','mOFC']]
    cout['group'] = 'control'
    pout['group'] = 'ptsd'
    out = pd.concat((cout,pout)).reset_index()
    out.to_csv('./posters&talks/%s_df.csv'%(name),index=False)


##############ERS with memory############################
beh = pd.read_csv('memory_behavior.csv')
for i in beh.index:
    if beh.loc[i,'phase'] == 'fear_conditioning':
        beh.loc[i,'phase'] = 'acquisition'
    elif beh.loc[i,'phase'] == 'extinction' and beh.loc[i,'block'] <= 1:
        beh.loc[i,'phase'] = 'early_extinction'
beh = beh.rename(columns={'phase':'encode_phase','condition':'trial_type'})
beh = beh.groupby(['subject','encode_phase','trial_type']).mean().reset_index()
beh['group'] = beh.subject.apply(lgroup)
beh = beh.set_index(['trial_type','group','encode_phase']).sort_index()


mem = 'high_confidence_accuracy'
cm = c.df.groupby(['trial_type',mem,'encode_phase','roi','subject']).mean()
pm = p.df.groupby(['trial_type',mem,'encode_phase','roi','subject']).mean()
memdf = pd.concat((cm,pm)).reset_index()

memdf = memdf.set_index(['trial_type',mem,'encode_phase','roi','subject'])
memdf = (memdf.loc['CS+'] - memdf.loc['CS-'])
memdf = (memdf.loc[1] - memdf.loc[0]).reset_index()
memdf['group'] = memdf.subject.apply(lgroup)
memdf = memdf.set_index(['group','roi','encode_phase']).sort_index()
# cscomp('control',memdf,['mOFC'],phases=phases)
# cdf['group'] = cdf.subject.apply(lgroup)

# pdf = p.df.groupby(['trial_type','encode_phase','roi',mem,'subject']).mean()
# pdf = (pdf.loc['CS+'] - pdf.loc['CS-']).reset_index()
# pdf['group'] = pdf.subject.apply(lgroup)

# mdf = pd.concat((cdf,pdf))
# mdf.roi = mdf.roi.apply(lambda x: 'vmPFC' if x == 'mOFC' else x)
# mdf.encode_phase = mdf.encode_phase.apply(lambda x: 'conditioning' if x == 'fear_conditioning' else x)
# if 'high' in mem: mdf[mem] = mdf[mem].apply(lambda x: 'hit' if x == 1 else 'miss')
# else: mdf[mem] = mdf[mem].apply(lambda x: 'hit' if x in [1,2] else 'miss')
# mdf = mdf.set_index(['group','roi','encode_phase',mem])

# mem_cscomp('control',mdf,['fvmPFC','fdACC'],phases=phases)
# mem_cscomp('ptsd',mdf,['fvmPFC','fdACC'],phases=phases)
###############hits vs. miss samephase with cs#####################
# mem = 'high_confidence_accuracy'
# cdf = c.df.groupby([mem,'trial_type','encode_phase','roi','subject']).mean()
# cdf = (cdf.loc[1] - cdf.loc[0]).reset_index()
# cdf['group'] = cdf.subject.apply(lgroup)

# pdf = p.df.groupby([mem,'trial_type','encode_phase','roi','subject']).mean()
# pdf = (pdf.loc[1] - pdf.loc[0]).reset_index()
# pdf['group'] = pdf.subject.apply(lgroup)

# hm = pd.concat([cdf,pdf]).reset_index(drop=True)
# roi = 'fvmPFC'
# g = sns.catplot(x='encode_phase',y='rsa',hue='trial_type',col='group',data=hm.query('roi == @roi'),
#                 kind='bar',aspect=1,palette=cpal)
# g.set_titles("{col_name} {col_var}")

##############hit count####################################
# beh = pd.read_csv(os.path.join('all_template_df.csv'),index_col=0)
#             for cs in self.conditions:
#                 con = '%s_trial'%(self.conditions[cs])
#                 for i in range(1,9): 
#                     self.df.loc[ self.df[ self.df.encode_phase == 'extinction' ][ self.df[con] == i ].index,'encode_phase' ] = 'early_extinction'


# beh = beh.groupby(['subject','encode','trial_type']).sum().reset_index()
# beh['group'] = beh.subject.apply(lgroup)
# beh.encode = pd.Categorical(beh.encode,categories=['baseline','fear_conditioning','extinction'],ordered=True)
# beh[['hc_acc','acc']] /= 24

# for acc in ['hc_acc','acc']:
#     fig, ax = plt.subplots(1,2,sharey=True)
#     for i, group in enumerate(['control','ptsd']):
#         dat = beh.query('group == @group')
#         sns.boxplot(data=dat,x='encode',y=acc,hue='trial_type',palette=cpal,ax=ax[i])
#         sns.swarmplot(data=dat,x='encode',y=acc,hue='trial_type',palette=cpal,ax=ax[i],
#                     linewidth=2,edgecolor='black',dodge=True,size=5)
#     # ax.set_yticks(range(0,25,4))
#         ax[i].yaxis.set_major_locator(MultipleLocator(.2))
#         ax[i].yaxis.set_minor_locator(MultipleLocator(.1))
#         ax[i].set_ylim(0,1.01)
#         ax[i].legend_.remove()
#         ax[i].set_title(group)

# beh = beh.set_