from roi_rsa import *
from scipy.interpolate import interp1d
from fg_config import *

sns.set_context('talk')
sns.set_style('ticks')

def cscomp(group,df,rois,n_boot=1000,phases=None):
    df = df.loc[group]
    out = {}
    for val in ['ci','dist']: out[val] = {}
    for roi in rois:
        out['ci'][roi] = {}
        out['dist'][roi] = {}
        for phase in phases:
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
    ci.encode_phase = pd.Categorical(ci.encode_phase,phases,ordered=True)
    ci = ci.sort_values(by='encode_phase')
    ci = ci.set_index('roi').sort_index()


    dist.roi = pd.Categorical(dist.roi,rois,ordered=True)
    dist.encode_phase = pd.Categorical(dist.encode_phase,phases,ordered=True)
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

def mem_cscomp(group,df,rois,n_boot=1000,phases=None):
    mems = ['hit','miss']

    df = df.loc[group]
    out = {}
    for val in ['ci','dist']: out[val] = {}
    for roi in rois:
        out['ci'][roi] = {}
        out['dist'][roi] = {}
        for phase in phases:
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
    ci.encode_phase = pd.Categorical(ci.encode_phase,phases,ordered=True)
    ci.mem  = pd.Categorical(ci.mem,mems,ordered=True)
    ci = ci.set_index('roi').sort_index()

    dist.roi = pd.Categorical(dist.roi,rois,ordered=True)
    dist.encode_phase = pd.Categorical(dist.encode_phase,phases,ordered=True)
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