from fg_config import *
from roi_rsa import *
c = group_roi_rsa(group='control',ext_split=True,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=True,fs=True,hemi=False)

memcon = ['encoding','retrieval']
mems = ['hit','miss']
cons = ['CS+','CS-']
rois = ['fvmPFC','fdACC','mOFC', 'dACC', 'amyg', 'hpc', 'ins']
phases = {'baseline':24,'acquisition':24,'early_extinction':8,'extinction':16}
# phases = {'baseline':24,'acquisition':24,'extinction':24}
subs = range(24)

from graphing_functions import *
##############Item level##############################
cdf = c.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
cdf = (cdf.loc['CS+'] - cdf.loc['CS-']).reset_index()
cdf['group'] = cdf.subject.apply(lgroup)

pdf = p.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = (pdf.loc['CS+'] - pdf.loc['CS-']).reset_index()

pdf['group'] = pdf.subject.apply(lgroup)

mdf = pd.concat((cdf,pdf))
# mdf.roi = mdf.roi.apply(lambda x: 'vmPFC' if x == 'mOFC' else x)

mdf = mdf.set_index(['group','roi','encode_phase'])
cscomp('control',mdf,['mOFC','fdACC'],phases=phases.keys())
cscomp('ptsd',mdf,['mOFC','fdACC'],phases=phases.keys())

##############Set level##############################
csl = pd.DataFrame(index=pd.MultiIndex.from_product([cons,rois,phases,subs],names=['condition','roi','encode_phase','subject']))
psl = pd.DataFrame(index=pd.MultiIndex.from_product([cons,rois,phases,subs],names=['condition','roi','encode_phase','subject']))
for con in cons:
    for roi in rois:
        for phase in phases:
            for sub in subs:
                csl.loc[(con,roi,phase,sub),'rsa'] = c.mats[roi][sub,slices[con][phase]['encoding'],slices[con][phase]['retrieval']][np.eye(phases[phase]) == 0].mean()
                psl.loc[(con,roi,phase,sub),'rsa'] = p.mats[roi][sub,slices[con][phase]['encoding'],slices[con][phase]['retrieval']][np.eye(phases[phase]) == 0].mean()

csl = (csl.loc['CS+'] - csl.loc['CS-']).reset_index()
psl = (psl.loc['CS+'] - psl.loc['CS-']).reset_index()
csl.subject = csl.subject.astype(int)
psl.subject = psl.subject.astype(int) + 101
sl = pd.concat([csl,psl])
sl['group'] = sl.subject.apply(lgroup)

sl = sl.set_index(['group','roi','encode_phase'])
cscomp('control',sl,['fvmPFC','fdACC'],phases=phases.keys())
cscomp('ptsd',sl,['fvmPFC','fdACC'],phases=phases.keys())

##############Within session similarity###############
cws = pd.DataFrame(index=pd.MultiIndex.from_product([memcon,cons,rois,phases,subs],names=['memory_phase','condition','roi','encode_phase','subject']))
pws = pd.DataFrame(index=pd.MultiIndex.from_product([memcon,cons,rois,phases,subs],names=['memory_phase','condition','roi','encode_phase','subject']))
for mem in memcon:
    for con in cons:
        for roi in rois:
            for phase in phases:
                for sub in subs:
                    cws.loc[(mem,con,roi,phase,sub),'rsa'] = c.mats[roi][sub,slices[con][phase][mem],slices[con][phase][mem]][np.tril_indices(phases[phase],-1)].mean()
                    pws.loc[(mem,con,roi,phase,sub),'rsa'] = p.mats[roi][sub,slices[con][phase][mem],slices[con][phase][mem]][np.tril_indices(phases[phase],-1)].mean()

cws = cws.reset_index()
pws = pws.reset_index()
cws.subject = cws.subject.astype(int)
pws.subject = pws.subject.astype(int) + 101
ws = pd.concat([cws,pws])
ws['group'] = ws.subject.apply(lgroup)
ws['gvar'] = ws.condition + ws.memory_phase
pal = sns.color_palette(['orangered','slategrey','orange','lightgrey'])

def mem_phase(roi):
    g = sns.catplot(x='encode_phase',y='rsa',hue='gvar',col='group',data=ws.query('roi == @roi'),
                    kind='bar',aspect=1,palette=pal)
    g.set_titles("{col_name} {col_var}")


##############ERS with memory############################
mem = 'high_confidence_accuracy'
cdf = c.df.groupby(['trial_type','encode_phase','roi',mem,'subject']).mean()
cdf = (cdf.loc['CS+'] - cdf.loc['CS-']).reset_index()
cdf['group'] = cdf.subject.apply(lgroup)

pdf = p.df.groupby(['trial_type','encode_phase','roi',mem,'subject']).mean()
pdf = (pdf.loc['CS+'] - pdf.loc['CS-']).reset_index()
pdf['group'] = pdf.subject.apply(lgroup)

mdf = pd.concat((cdf,pdf))
mdf.roi = mdf.roi.apply(lambda x: 'vmPFC' if x == 'mOFC' else x)
mdf.encode_phase = mdf.encode_phase.apply(lambda x: 'conditioning' if x == 'fear_conditioning' else x)
if 'high' in mem: mdf[mem] = mdf[mem].apply(lambda x: 'hit' if x == 1 else 'miss')
else: mdf[mem] = mdf[mem].apply(lambda x: 'hit' if x in [1,2] else 'miss')
mdf = mdf.set_index(['group','roi','encode_phase',mem])

mem_cscomp('control',mdf,['fvmPFC','fdACC'],phases=phases)
mem_cscomp('ptsd',mdf,['fvmPFC','fdACC'],phases=phases)
###############hits vs. miss samephase with cs#####################
mem = 'high_confidence_accuracy'
cdf = c.df.groupby([mem,'trial_type','encode_phase','roi','subject']).mean()
cdf = (cdf.loc[1] - cdf.loc[0]).reset_index()
cdf['group'] = cdf.subject.apply(lgroup)

pdf = p.df.groupby([mem,'trial_type','encode_phase','roi','subject']).mean()
pdf = (pdf.loc[1] - pdf.loc[0]).reset_index()
pdf['group'] = pdf.subject.apply(lgroup)

hm = pd.concat([cdf,pdf]).reset_index(drop=True)
roi = 'fvmPFC'
g = sns.catplot(x='encode_phase',y='rsa',hue='trial_type',col='group',data=hm.query('roi == @roi'),
                kind='bar',aspect=1,palette=cpal)
g.set_titles("{col_name} {col_var}")

##############hit count####################################
df = pd.read_csv(os.path.join('all_template_df.csv'),index_col=0)
df = df.groupby(['subject','encode','trial_type']).sum().reset_index()
df['group'] = df.subject.apply(lgroup)
df.encode = pd.Categorical(df.encode,categories=['baseline','fear_conditioning','extinction'],ordered=True)
df[['hc_acc','acc']] /= 24
for acc in ['hc_acc','acc']:
    fig, ax = plt.subplots(1,2,sharey=True)
    for i, group in enumerate(['control','ptsd']):
        dat = df.query('group == @group')
        sns.boxplot(data=dat,x='encode',y=acc,hue='trial_type',palette=cpal,ax=ax[i])
        sns.swarmplot(data=dat,x='encode',y=acc,hue='trial_type',palette=cpal,ax=ax[i],
                    linewidth=2,edgecolor='black',dodge=True,size=5)
    # ax.set_yticks(range(0,25,4))
        ax[i].yaxis.set_major_locator(MultipleLocator(.2))
        ax[i].yaxis.set_minor_locator(MultipleLocator(.1))
        ax[i].set_ylim(0,1.01)
        ax[i].legend_.remove()
        ax[i].set_title(group)
