from roi_rsa import *
from scipy.interpolate import interp1d
from fg_config import *

# sns.set_style({'axes.facecolor':'.9','figure.facecolor':'.9'})
sns.set_style({'axes.facecolor':'1','figure.facecolor':'1'})
idx = pd.IndexSlice
c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)

cdf = c.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
cdf = (cdf.loc['CS+'] - cdf.loc['CS-']).reset_index()
cdf['group'] = cdf.subject.apply(lgroup)

pdf = p.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = (pdf.loc['CS+'] - pdf.loc['CS-']).reset_index()

pdf['group'] = pdf.subject.apply(lgroup)

pfc = pd.concat((cdf,pdf))
pfc.roi = pfc.roi.apply(lambda x: 'vmPFC' if x == 'mOFC' else x)
# pfc.encode = pfc.encode.apply(lambda x: 'conditioning' if x == 'fear_conditioning' else x)

pfc = pfc.set_index(['group','roi','encode_phase'])

encode = ['baseline','acquisition','early_extinction','extinction']
# encode = ['baseline','acquisition','extinction']
def cscomp(group,df,rois,n_boot=1000):
    df = df.loc[group]
    out = {}
    for val in ['ci','dist']: out[val] = {}
    for roi in rois:
        out['ci'][roi] = {}
        out['dist'][roi] = {}
        for phase in encode:
            out['ci'][roi][phase] = {}
            out['dist'][roi][phase] = {}

            out['ci'][roi][phase], out['dist'][roi][phase] = pg.compute_bootci(df.loc[(roi,phase),'rsa'].values,
                                                                func='mean',n_boot=n_boot,return_dist=True,
                                                                method='cper',decimals=3,seed=42)
    
    ci = pd.DataFrame.from_dict(out['ci'],orient='index').reset_index().rename(columns={'index':'roi'})
    ci = ci.melt(id_vars='roi',var_name='encode_phase',value_name='ci').set_index(['roi','encode_phase']).sort_index()

    dist = pd.DataFrame.from_dict(out['dist'],orient='index').reset_index().rename(columns={'index':'roi'})
    dist = dist.melt(id_vars='roi',var_name='encode_phase',value_name='dist').set_index(['roi','encode_phase'])
    dist = dist.dist.apply(pd.Series).stack().reset_index(-1,drop=True)
    dist = dist.reset_index().rename(columns={0:'dist'})
    ci['point'] = dist.groupby(['roi','encode_phase']).mean()
    ci = ci.reset_index()
    
    print(ci.columns)
    ci.roi = pd.Categorical(ci.roi,rois,ordered=True)
    ci.encode_phase = pd.Categorical(ci.encode_phase,encode,ordered=True)
    ci = ci.sort_values(by='encode_phase')
    ci = ci.set_index('roi').sort_index()


    dist.roi = pd.Categorical(dist.roi,rois,ordered=True)
    dist.encode_phase = pd.Categorical(dist.encode_phase,encode,ordered=True)
    ci = ci.sort_values(by='encode_phase')
    dist = dist.set_index('roi').sort_index()

    phase_pal = sns.color_palette(['black','darkmagenta','lightgreen','seagreen'],desat=.75)
    fig, ax = plt.subplots(1,len(rois),sharey=True,figsize=(12,8))
    for i, roi in enumerate(rois):
        if len(rois) == 1: Ax = ax
        else: Ax = ax[i]
        sns.violinplot(data=dist.loc[roi],x='encode_phase',y='dist',
                        inner=None,ax=Ax,scale='count',palette=phase_pal)
        lower = ci.loc[roi,'ci'].apply(lambda x: x[0])
        upper = ci.loc[roi,'ci'].apply(lambda x: x[1])
        Y = ci.loc[roi,'point']
        X = Ax.get_xticks()
        Ax.vlines(X,lower,upper,linewidth=3,color='white').set_capstyle('round')
        Ax.scatter(X,Y,s=50,color='white')
        Ax.hlines(0,Ax.get_xlim()[0],Ax.get_xlim()[1],color='grey',linestyle='--',linewidth=3)
        Ax.set_xticklabels('',rotation=45)
        Ax.set_title(group+'_'+roi)
    if len(rois) > 1:
        ax[0].set_ylabel('∆ fisher z(r)')
        ax[1].set_ylabel('')
    else:
        ax.set_ylabel('∆ fisher z(r)')
        ax.set_ylabel('')

cscomp('control',pfc,['vmPFC','dACC'])
cscomp('ptsd',pfc,['vmPFC','dACC'])
cscomp('control',pfc,['ins'])
cscomp('ptsd',pfc,['ins'])
cscomp('control',pfc,['amyg','hpc'])
cscomp('ptsd',pfc,['amyg','hpc'])

cscomp('control',pfc,['rh_hpc','lh_hpc'])
cscomp('ptsd',pfc,['rh_hpc','lh_hpc'])

cscomp('control',pfc,['rh_amyg','lh_amyg'])
cscomp('ptsd',pfc,['rh_amyg','lh_amyg'])

# cscomp('control',pfc,['amyg_cem','amyg_bla'])
# cscomp('ptsd',pfc,['amyg_cem','amyg_bla'])
# cscomp('control',pfc,['hc_head','hc_body','hc_tail'])
# cscomp('ptsd',pfc,['hc_head','hc_body','hc_tail'])

##############################################################################
cdf = c.df.groupby(['encode_phase','trial_type','roi','subject']).mean()
cdf = (cdf.loc['acquisition'] - cdf.loc['extinction']).reset_index()
cdf['group'] = cdf.subject.apply(lgroup)

pdf = p.df.groupby(['encode_phase','trial_type','roi','subject']).mean()
pdf = (pdf.loc['acquisition'] - pdf.loc['extinction']).reset_index()
pdf['group'] = pdf.subject.apply(lgroup)

pfc = pd.concat((cdf,pdf))
pfc.roi = pfc.roi.apply(lambda x: 'vmPFC' if x == 'mOFC' else x)

pfc = pfc.set_index(['group','roi','trial_type'])

cons = ['CS+','CS-']
def phasecomp(group,df,rois,n_boot=1000):
    df = df.loc[group]
    out = {}
    for val in ['ci','dist']: out[val] = {}
    for roi in rois:
        out['ci'][roi] = {}
        out['dist'][roi] = {}
        for cs in cons:
            out['ci'][roi][cs] = {}
            out['dist'][roi][cs] = {}

            out['ci'][roi][cs], out['dist'][roi][cs] = pg.compute_bootci(df.loc[(roi,cs),'rsa'].values,
                                                                func='mean',n_boot=n_boot,return_dist=True,
                                                                method='cper',decimals=3,seed=420)
    ci = pd.DataFrame.from_dict(out['ci'],orient='index').reset_index().rename(columns={'index':'roi'})
    ci = ci.melt(id_vars='roi',var_name='trial_type',value_name='ci').set_index(['roi','trial_type']).sort_index()

    dist = pd.DataFrame.from_dict(out['dist'],orient='index').reset_index().rename(columns={'index':'roi'})
    dist = dist.melt(id_vars='roi',var_name='trial_type',value_name='dist').set_index(['roi','trial_type'])
    dist = dist.dist.apply(pd.Series).stack().reset_index(-1,drop=True)
    dist = dist.reset_index().rename(columns={0:'dist'})
    ci['point'] = dist.groupby(['roi','trial_type']).mean()
    ci = ci.reset_index()

    ci.roi = pd.Categorical(ci.roi,rois,ordered=True)
    ci.trial_type = pd.Categorical(ci.trial_type,cons,ordered=True)

    dist.roi = pd.Categorical(dist.roi,rois,ordered=True)
    dist.trial_type = pd.Categorical(dist.trial_type,cons,ordered=True)

    ci = ci.sort_values(by=['roi','trial_type'])

    fig, ax = plt.subplots()
    sns.violinplot(data=dist,x='roi',y='dist',hue='trial_type',split=True,
                    inner=None,ax=ax,scale='count',palette=cpal)
    lower = ci['ci'].apply(lambda x: x[0])
    upper = ci['ci'].apply(lambda x: x[1])
    Y = ci['point']
    X = [[x-.05,x+.05] for x in ax.get_xticks()]
    X = [x for t in X for x in t]
    ax.vlines(X,lower,upper,linewidth=3,color='white').set_capstyle('round')
    ax.scatter(X,Y,s=50,color='white')
    ax.hlines(0,ax.get_xlim()[0],ax.get_xlim()[1],color='black',linestyle='--',linewidth=3)
    ax.legend_.remove()
    ax.set_title(group)
    ax.set_ylabel('∆ fisher z(r)')

phasecomp('control',pfc,['vmPFC','dACC'])
phasecomp('ptsd',pfc,['vmPFC','dACC'])
phasecomp('control',pfc,['ins'])
phasecomp('ptsd',pfc,['ins'])
phasecomp('control',pfc,['amyg','hpc'])
phasecomp('ptsd',pfc,['amyg','hpc'])


# phasecomp('control',pfc,['amyg_cem','amyg_bla'])
# phasecomp('ptsd',pfc,['amyg_cem','amyg_bla'])
# phasecomp('control',pfc,['hc_head','hc_body','hc_tail'])
# phasecomp('ptsd',pfc,['hc_head','hc_body','hc_tail'])

#################################################################################
c = group_roi_rsa(group='control',ext_split=True,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=True,fs=True,hemi=False)

mem = 'low_confidence_accuracy'
cdf = c.df.groupby(['trial_type','encode_phase','roi',mem,'subject']).mean()
cdf = (cdf.loc['CS+'] - cdf.loc['CS-']).reset_index()
cdf['group'] = cdf.subject.apply(lgroup)

pdf = p.df.groupby(['trial_type','encode_phase','roi',mem,'subject']).mean()
pdf = (pdf.loc['CS+'] - pdf.loc['CS-']).reset_index()
pdf['group'] = pdf.subject.apply(lgroup)

pfc = pd.concat((cdf,pdf))
pfc.roi = pfc.roi.apply(lambda x: 'vmPFC' if x == 'mOFC' else x)
pfc.encode_phase = pfc.encode_phase.apply(lambda x: 'conditioning' if x == 'fear_conditioning' else x)
if 'high' in mem: pfc[mem] = pfc[mem].apply(lambda x: 'hit' if x == 1 else 'miss')
else: pfc[mem] = pfc[mem].apply(lambda x: 'hit' if x in [1,2] else 'miss')

pfc = pfc.set_index(['group','roi','encode_phase',mem])

mems = ['hit','miss']
def mem_cscomp(group,df,rois,n_boot=1000):

    df = df.loc[group]
    out = {}
    for val in ['ci','dist']: out[val] = {}
    for roi in rois:
        out['ci'][roi] = {}
        out['dist'][roi] = {}
        for phase in encode:
            out['ci'][roi][phase] = {}
            out['dist'][roi][phase] = {}
            for mem in mems:
                out['ci'][roi][phase][mem] = {}
                out['dist'][roi][phase][mem] = {}
                
                data = df.loc[(roi,phase,mem),'rsa'].dropna().values
                out['ci'][roi][phase][mem], out['dist'][roi][phase][mem] = pg.compute_bootci(data,
                                                                func='mean',n_boot=n_boot,return_dist=True,
                                                                method='cper',decimals=3,seed=42)
    ci = pd.DataFrame.from_dict(out['ci'],orient='index').reset_index().rename(columns={'index':'roi'})
    ci = ci.melt(id_vars='roi',var_name='encode_phase',value_name='ci').set_index(['roi','encode_phase']).sort_index()
    ci = ci.ci.apply(pd.Series).reset_index()
    ci = ci.melt(id_vars=['roi','encode_phase'],var_name='mem',value_name='ci').set_index(['roi','encode_phase','mem'])
    
    dist = pd.DataFrame.from_dict(out['dist'],orient='index').reset_index().rename(columns={'index':'roi'})
    dist = dist.melt(id_vars='roi',var_name='encode_phase',value_name='dist').set_index(['roi','encode_phase'])
    dist = dist.dist.apply(pd.Series).reset_index()
    dist = dist.melt(id_vars=['roi','encode_phase'],var_name='mem',value_name='dist').set_index(['roi','encode_phase','mem'])
    dist = dist.dist.apply(pd.Series).stack().reset_index(-1,drop=True)

    dist = dist.reset_index().rename(columns={0:'dist'})
    ci['point'] = dist.groupby(['roi','encode_phase','mem']).mean()
    ci = ci.reset_index()

    ci.roi = pd.Categorical(ci.roi,rois,ordered=True)
    ci.encode_phase = pd.Categorical(ci.encode_phase,encode,ordered=True)
    ci.mem  = pd.Categorical(ci.mem,mems,ordered=True)
    ci = ci.set_index('roi').sort_index()

    dist.roi = pd.Categorical(dist.roi,rois,ordered=True)
    dist.encode_phase = pd.Categorical(dist.encode_phase,encode,ordered=True)
    dist.mem  = pd.Categorical(dist.mem,mems,ordered=True)
    dist = dist.set_index('roi').sort_index()

    ci = ci.sort_values(by=['roi','encode_phase','mem'])


    # hit_pal = sns.color_palette(['black','darkmagenta','lightgreen','seagreen'],desat=.75)
    # miss_pal = sns.color_palette(['midnightblue','darkmagenta','lightgreen','seagreen'],desat=.25)
    # phase_pal = np.empty(8,dtype=object)
    # phase_pal[0::2] = hit_pal
    # phase_pal[1::2] = miss_pal
    mem_pal = ['darkblue','lightblue']

    fig, ax = plt.subplots(1,len(rois),sharey=True)
    for i, roi in enumerate(rois):
        if len(rois) == 1: Ax = ax
        else: Ax = ax[i]
        sns.violinplot(data=dist.loc[roi],x='encode_phase',y='dist',hue='mem',split=True,
                        inner=None,ax=Ax,scale='count',palette=mem_pal,scale_hue=True)
        lower = ci.loc[roi,'ci'].apply(lambda x: x[0])
        upper = ci.loc[roi,'ci'].apply(lambda x: x[1])
        Y = ci.loc[roi,'point']
        X = [[x-.05,x+.05] for x in Ax.get_xticks()]
        X = [x for t in X for x in t]
        Ax.vlines(X,lower,upper,linewidth=3,color='white').set_capstyle('round')
        Ax.scatter(X,Y,s=50,color='white')
        Ax.hlines(0,Ax.get_xlim()[0],Ax.get_xlim()[1],color='grey',linestyle='--',linewidth=3)
        Ax.set_xticklabels('',rotation=45)
        Ax.set_title(group+'_'+roi)
        Ax.legend_.remove()
    if len(rois) > 1:
        ax[0].set_ylabel('∆ fisher z(r)')
        ax[1].set_ylabel('')
    else:
        ax.set_ylabel('∆ fisher z(r)')
        ax.set_ylabel('')
mem_cscomp('healthy',pfc,['sgACC','rSMA'])
mem_cscomp('ptsd',pfc,['sgACC','rSMA'])
mem_cscomp('control',pfc,['ins'])
mem_cscomp('ptsd',pfc,['ins'])
mem_cscomp('control',pfc,['amyg','hpc'])
mem_cscomp('ptsd',pfc,['amyg','hpc'])


# mem_cscomp('control',pfc,['amyg_cem','amyg_bla'])
# mem_cscomp('ptsd',pfc,['amyg_cem','amyg_bla'])
# mem_cscomp('control',pfc,['hc_head','hc_body','hc_tail'])
# mem_cscomp('ptsd',pfc,['hc_head','hc_body','hc_tail'])

#################################################################################
cdf = c.df.groupby(['encode_phase','trial_type','roi','high_confidence_accuracy','subject']).mean()
cdf = (cdf.loc['acquisition'] - cdf.loc['extinction']).reset_index()
cdf['group'] = cdf.subject.apply(lgroup)

pdf = p.df.groupby(['encode_phase','trial_type','roi','high_confidence_accuracy','subject']).mean()
pdf = (pdf.loc['acquisition'] - pdf.loc['extinction']).reset_index()
pdf['group'] = pdf.subject.apply(lgroup)

pfc = pd.concat((cdf,pdf))
pfc.roi = pfc.roi.apply(lambda x: 'vmPFC' if x == 'mOFC' else x)
pfc.high_confidence_accuracy = pfc.high_confidence_accuracy.apply(lambda x: 'hit' if x == 1 else 'miss')
pfc = pfc.set_index(['group','roi','trial_type','high_confidence_accuracy'])

def mem_phasecomp(group,df,rois,n_boot=1000):
    df = df.loc[group]
    out = {}
    for val in ['ci','dist']: out[val] = {}
    for roi in rois:
        out['ci'][roi] = {}
        out['dist'][roi] = {}
        for cs in cons:
            out['ci'][roi][cs] = {}
            out['dist'][roi][cs] = {}
            for mem in mems:
                out['ci'][roi][cs][mem] = {}
                out['dist'][roi][cs][mem] = {}
                
                data = df.loc[(roi,cs,mem),'rsa'].dropna().values
                out['ci'][roi][cs][mem], out['dist'][roi][cs][mem] = pg.compute_bootci(data,
                                                                func='mean',n_boot=n_boot,return_dist=True,
                                                                method='cper',decimals=3,seed=42)
    ci = pd.DataFrame.from_dict(out['ci'],orient='index').reset_index().rename(columns={'index':'roi'})
    ci = ci.melt(id_vars='roi',var_name='trial_type',value_name='ci').set_index(['roi','trial_type']).sort_index()
    ci = ci.ci.apply(pd.Series).reset_index()
    ci = ci.melt(id_vars=['roi','trial_type'],var_name='mem',value_name='ci').set_index(['roi','trial_type','mem'])
    
    dist = pd.DataFrame.from_dict(out['dist'],orient='index').reset_index().rename(columns={'index':'roi'})
    dist = dist.melt(id_vars='roi',var_name='trial_type',value_name='dist').set_index(['roi','trial_type'])
    dist = dist.dist.apply(pd.Series).reset_index()
    dist = dist.melt(id_vars=['roi','trial_type'],var_name='mem',value_name='dist').set_index(['roi','trial_type','mem'])
    dist = dist.dist.apply(pd.Series).stack().reset_index(-1,drop=True)

    dist = dist.reset_index().rename(columns={0:'dist'})
    ci['point'] = dist.groupby(['roi','trial_type','mem']).mean()
    ci = ci.reset_index()

    ci.roi = pd.Categorical(ci.roi,rois,ordered=True)
    ci.trial_type = pd.Categorical(ci.trial_type,cons,ordered=True)
    ci.mem  = pd.Categorical(ci.mem,mems,ordered=True)
    ci.reset_index()

    dist.roi = pd.Categorical(dist.roi,rois,ordered=True)
    dist.trial_type = pd.Categorical(dist.trial_type,cons,ordered=True)
    dist.mem  = pd.Categorical(dist.mem,mems,ordered=True)
    ci.reset_index()

    ci = ci.sort_values(by=['roi','trial_type','mem']).reset_index(drop=True)
    dist = dist.sort_values(by=['roi','trial_type','mem']).reset_index(drop=True)

    ci['roi_cs'] = ''
    for i in ci.index: ci.loc[i,'roi_cs'] = ci.loc[i,'roi'] + '_' + ci.loc[i,'trial_type'] 
    dist['roi_cs'] = ''
    for i in dist.index: dist.loc[i,'roi_cs'] = dist.loc[i,'roi'] + '_' + dist.loc[i,'trial_type'] 

    mem_pal = ['darkblue','lightblue']
    fig, ax = plt.subplots()
    sns.violinplot(data=dist,x='roi_cs',y='dist',hue='mem',split=True,
                    inner=None,ax=ax,scale='count',palette=mem_pal)
    lower = ci['ci'].apply(lambda x: x[0])
    upper = ci['ci'].apply(lambda x: x[1])
    Y = ci['point']
    X = [[x-.05,x+.05] for x in ax.get_xticks()]
    X = [x for t in X for x in t]
    ax.vlines(X,lower,upper,linewidth=3,color='white').set_capstyle('round')
    ax.scatter(X,Y,s=50,color='white')
    ax.hlines(0,ax.get_xlim()[0],ax.get_xlim()[1],color='black',linestyle='--',linewidth=3)
    ax.legend_.remove()
    ax.set_title(group)
    ax.set_ylabel('∆ fisher z(r)')

mem_phasecomp('control',pfc,['vmPFC','dACC'])
mem_phasecomp('ptsd',pfc,['vmPFC','dACC'])
# mem_phasecomp('control',pfc,['insula'])
# mem_phasecomp('ptsd',pfc,['insula'])
#####################################################################3#
cdf = c.df.groupby(['trial_type','encode_phase','roi','high_confidence_accuracy','subject']).mean()
cdf = (cdf.loc['CS+'] - cdf.loc['CS-']).reset_index()
cdf['group'] = cdf.subject.apply(lgroup)

pdf = p.df.groupby(['trial_type','encode_phase','roi','high_confidence_accuracy','subject']).mean()
pdf = (pdf.loc['CS+'] - pdf.loc['CS-']).reset_index()
pdf['group'] = pdf.subject.apply(lgroup)

pfc = pd.concat((cdf,pdf))
pfc.roi = pfc.roi.apply(lambda x: 'vmPFC' if x == 'mOFC' else x)
pfc.encode_phase = pfc.encode_phase.apply(lambda x: 'conditioning' if x == 'fear_conditioning' else x)
pfc.high_confidence_accuracy = pfc.high_confidence_accuracy.apply(lambda x: 'hit' if x == 1 else 'miss')

pfc = pfc.set_index(['roi','group','encode_phase','low_confidence_accuracy'])
groups = ['healthy','ptsd']
mems = ['hit','miss']
def mem_group_comp(roi,df,n_boot=1000):
    roi = 'sgACC'
    phase = 'extinction'
    df = df.loc[roi]
    out = {}
    for val in ['ci','dist']: out[val] = {}
    for group in groups:
        out['ci'][group] = {}
        out['dist'][group] = {}
        for mem in mems:
            out['ci'][group][mem] = {}
            out['dist'][group][mem] = {}
                
            data = df.loc[(group,phase,mem),'rsa'].dropna().values
            out['ci'][group][mem], out['dist'][group][mem] = pg.compute_bootci(data,
                                                                func='mean',n_boot=n_boot,return_dist=True,
                                                                method='cper',decimals=3,seed=42)
    ci = pd.DataFrame.from_dict(out['ci'],orient='index').reset_index().rename(columns={'index':'group'})
    ci = ci.melt(id_vars='group',var_name='mem',value_name='ci').set_index(['group','mem']).sort_index()
    # ci = ci.ci.apply(pd.Series).reset_index()
    # ci = ci.melt(id_vars=['group','mem'],var_name='mem',value_name='ci').set_index(['roi','encode','mem'])
    
    dist = pd.DataFrame.from_dict(out['dist'],orient='index').reset_index().rename(columns={'index':'group'})
    dist = dist.melt(id_vars='group',var_name='mem',value_name='dist').set_index(['group','mem'])
    # dist = dist.dist.apply(pd.Series).reset_index()
    # dist = dist.melt(id_vars=['roi','encode'],var_name='mem',value_name='dist').set_index(['roi','encode','mem'])
    dist = dist.dist.apply(pd.Series).stack().reset_index(-1,drop=True)

    dist = dist.reset_index().rename(columns={0:'dist'})
    ci['point'] = dist.groupby(['group','mem']).mean()
    ci = ci.reset_index()

    ci.group = pd.Categorical(ci.group,groups,ordered=True)
    # ci.encode = pd.Categorical(ci.encode,encode,ordered=True)
    ci.mem  = pd.Categorical(ci.mem,mems,ordered=True)
    # ci = ci.set_index('roi').sort_index()

    dist.group = pd.Categorical(dist.group,groups,ordered=True)
    # dist.encode = pd.Categorical(dist.encode,encode,ordered=True)
    dist.mem  = pd.Categorical(dist.mem,mems,ordered=True)
    # dist = dist.set_index('roi').sort_index()

    ci = ci.sort_values(by=['group','mem'])
    ci = ci.set_index('group')


    # hit_pal = sns.color_palette(['black','darkmagenta','lightgreen','seagreen'],desat=.75)
    # miss_pal = sns.color_palette(['midnightblue','darkmagenta','lightgreen','seagreen'],desat=.25)
    # phase_pal = np.empty(8,dtype=object)
    # phase_pal[0::2] = hit_pal
    # phase_pal[1::2] = miss_pal
    mem_pal = ['darkblue','lightblue']

    fig, Ax = plt.subplots()
    sns.violinplot(data=dist,x='group',y='dist',hue='mem',split=True,
                    inner=None,ax=Ax,scale='count',palette=mem_pal,scale_hue=True)
    lower = ci['ci'].apply(lambda x: x[0])
    upper = ci['ci'].apply(lambda x: x[1])
    Y = ci['point']
    X = [[x-.05,x+.05] for x in Ax.get_xticks()]
    X = [x for t in X for x in t]
    Ax.vlines(X,lower,upper,linewidth=3,color='white').set_capstyle('round')
    Ax.scatter(X,Y,s=50,color='white')
    Ax.hlines(0,Ax.get_xlim()[0],Ax.get_xlim()[1],color='grey',linestyle='--',linewidth=3)
    Ax.set_xticklabels('',rotation=45)
    Ax.set_title(group+'_'+roi)
    Ax.legend_.remove()
    Ax.set_ylabel('∆ fisher z(r)')
    Ax.set_xlabel('')







###################################################################
c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
c.cstrial()
cdf = c.csdf
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)
p.cstrial()
pdf = p.csdf

pfc = pd.concat((cdf,pdf)).reset_index()
pfc['group'] = pfc.subject.apply(lgroup)
pfc.roi = pfc.roi.apply(lambda x: 'vmPFC' if x == 'mOFC' else x)
# pfc.encode = pfc.encode.apply(lambda x: 'conditioning' if x == 'fear_conditioning' else x)
pfc = pfc.set_index(['roi','encode_phase','block','group']).sort_index()


groups = ['ptsd','control']
cspal = [gpal[1],gpal[0]] 
# groups = ['control']
# cspal = [gpal[0]] 
blocks = np.array(range(1,7))
bint = np.linspace(blocks.min(),blocks.max(),1000,endpoint=True)
phases = ['baseline','acquisition','extinction']
def interp(x,point,ci):
    if 'series' in str(type(ci)):
      lower = ci.apply(pd.Series)[0].values
      upper = ci.apply(pd.Series)[1].values

    if 'series' in str(type(point)): point = point.values
    if type(x) is range: x = np.array(x)    
    x_interp = np.linspace(x.min(),x.max(),1000,endpoint=True)

    up_spliner = interp1d(x,upper,kind='cubic')
    lo_spliner = interp1d(x,lower,kind='cubic')
    mean_spliner = interp1d(x,point,kind='cubic')

    upper_interp = up_spliner(x_interp)
    lower_interp = lo_spliner(x_interp)
    mean_interp = mean_spliner(x_interp)

    return mean_interp, lower_interp, upper_interp
def cstrial(roi,df,n_boot=1000):
    df = df.loc[roi]
    out = {}
    for val in ['ci','dist']: out[val] = {}
    for group in groups:
        out['ci'][group] = {}
        out['dist'][group] = {}
        for phase in phases:
            out['ci'][group][phase] = {}
            out['dist'][group][phase] = {}
            for block in blocks:
                out['ci'][group][phase][block] = {}
                out['dist'][group][phase][block] = {}
                
                out['ci'][group][phase][block], out['dist'][group][phase][block] = pg.compute_bootci(df.loc[(phase,block,group),'rsa'],
                                                                func='mean',n_boot=n_boot,return_dist=True,
                                                                method='cper',decimals=3,seed=42)
    ci = pd.DataFrame.from_dict(out['ci'],orient='index').reset_index().rename(columns={'index':'group'})
    ci = ci.melt(id_vars='group',var_name='encode_phase',value_name='ci').set_index(['group','encode_phase']).sort_index()
    ci = ci.ci.apply(pd.Series).reset_index()
    ci = ci.melt(id_vars=['group','encode_phase'],var_name='block',value_name='ci').set_index(['group','encode_phase','block'])
    
    dist = pd.DataFrame.from_dict(out['dist'],orient='index').reset_index().rename(columns={'index':'group'})
    dist = dist.melt(id_vars='group',var_name='encode_phase',value_name='dist').set_index(['group','encode_phase'])
    dist = dist.dist.apply(pd.Series).reset_index()
    dist = dist.melt(id_vars=['group','encode_phase'],var_name='block',value_name='dist').set_index(['group','encode_phase','block'])
    dist = dist.dist.apply(pd.Series).stack().reset_index(-1,drop=True)

    dist = dist.reset_index().rename(columns={0:'dist'})
    ci['point'] = dist.groupby(['group','encode_phase','block']).mean()
    ci = ci.reset_index()

    ci.group = pd.Categorical(ci.group,groups,ordered=True)
    ci.encode_phase = pd.Categorical(ci.encode_phase,phases,ordered=True)
    ci.block  = pd.Categorical(ci.block,blocks,ordered=True)

    dist.group = pd.Categorical(dist.group,groups,ordered=True)
    dist.encode_phase = pd.Categorical(dist.encode_phase,phases,ordered=True)
    dist.block  = pd.Categorical(dist.block,blocks,ordered=True)

    ci = ci.sort_values(by=['group','encode_phase','block']).reset_index(drop=True)
    dist = dist.sort_values(by=['group','encode_phase','block']).reset_index(drop=True)

    ci = ci.set_index(['group','encode_phase']).sort_index()
 
    fig, ax = plt.subplots(1,3,sharey=True)
    for i, phase in enumerate(phases):
        gp, gl, gu = {}, {}, {}
        for g, group in enumerate(groups):
            dat = ci.loc[(group,phase)]
            # get mean/sem per TR for current class
            gp[group],gl[group],gu[group] = interp(blocks,dat.point,dat.ci) 
            # draw
            # ax[i].plot(x_interp,gp[group],color=cspal[g],alpha=.8,linewidth=3)
            ax[i].fill_between(bint,gl[group],gu[group],
                            color=cspal[g],alpha=.5)
            ax[i].hlines(0,ax[i].get_xlim()[0],ax[i].get_xlim()[1],color='black',linestyle='--',linewidth=1)
            ax[i].set_xlim(1,6)
            ax[i].set_title(phase)
        # du = np.zeros(bint.shape[0])
        # np.where(gu['control'] > gu['ptsd'] > gl['control'])[0].shape
        # du = np.array([i for i in range(1000) if gu['control'][i] > gu['ptsd'][i] and gu['ptsd'][i] > gl['control'][i]])
        # dl = gl['control']
        # ax[i].fill_between(bint,dl,du,color='purple',alpha=1)

    ax[0].set_ylabel('CS+ - CS- relative ERS')
    plt.suptitle(roi)
cstrial('vmPFC',pfc,1000)
cstrial('dACC',pfc,1000)



#retroactive effect. baseline to conditioning overlap as a factor of encoding to retrieval
for group in ['control','ptsd']:
    r = group_roi_rsa('control',ext_split=False)
    r.mat_stats()

    df = pd.DataFrame.from_dict(r.stats,orient='index')
    df = df.apply(pd.Series).stack(
      ).apply(pd.Series).stack(
      ).reset_index()
    df = df.rename(columns={'level_0':'roi','level_1':'cs','level_2':'subject',0:'rsa'})
    df.roi = pd.Categorical(df.roi,['vmPFC','dACC','amyg_cem','amyg_bla','hc_head','hc_body','hc_tail'],ordered=True)
    fig, ax = plt.subplots()
    sns.pointplot(data=df,x='roi',y='rsa',hue='cs',ax=ax,dodge=True,join=False,palette=cpal)
    zero(ax)
    diff = df.set_index(['cs','subject','roi'])
    diff = (diff.loc['bc_CS+'] - diff.loc['bc_CS-']).reset_index()
    fig, ax = plt.subplots()
    sns.pointplot(data=diff,x='roi',y='rsa',ax=ax,join=False,color='purple')
    zero(ax)

    r.bc_logreg()
    gdf = pd.DataFrame.from_dict(r.bc_betas,orient='index').reset_index(
        ).melt(id_vars='index'
        ).rename(columns={'index':'roi','variable':'trial_type','value':'beta'})

    gdf = gdf.set_index(['roi', 'trial_type'])['beta'].apply(pd.Series).stack(
        ).reset_index().drop(columns='level_2').rename(columns={0:'beta'})
    # if not self.hemi: gdf.roi = pd.Categorical(gdf.roi,self.rois,ordered=True)
    gstats = gdf.groupby(['trial_type','roi']).mean()
    gstats['conf']  = gdf.groupby(['trial_type','roi']).apply(lambda x: np.percentile(x,[5,95]))
    fig, ax = plt.subplots()
    sns.pointplot(data=gdf,x='roi',y='beta',hue='trial_type',
              palette=cpal,dodge=True,capsize=.08,join=False,ax=ax,ci=None)
    zero(ax)
    ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
    for j, roi in enumerate(gdf.roi.unique()):
        x = ax.get_xticks()[j]
        ax.vlines(x-.015,ymin=gstats.loc[('CS+',roi),'conf'][0], ymax=gstats.loc[('CS+',roi),'conf'][1], color=cpal[0], linewidth=3).set_capstyle('round')
        ax.vlines(x+.015,ymin=gstats.loc[('CS-',roi),'conf'][0], ymax=gstats.loc[('CS-',roi),'conf'][1], color=cpal[1], linewidth=3).set_capstyle('round')
        # ax[i].errorbar(x+.05,gstats.loc[('CS-'),'beta'], yerr=np.concatenate(gstats.loc[('CS-'),'conf'].values).reshape(5,2).transpose(), color='black', linewidth=3,capsize=5,capthick=3)

