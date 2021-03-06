from roi_rsa import *
from scipy.interpolate import interp1d
from fg_config import *
from matplotlib.ticker import MultipleLocator, ScalarFormatter

sns.set_context('talk')
sns.set_style('ticks', {'axes.spines.right':False, 'axes.spines.top':False})
sns.set_style({'axes.facecolor':'1','figure.facecolor':'1'})

def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)


def pconvert(p):
    if p < .001:
        return '3'
    elif p < .01:
        return '2'
    elif p < .05:
        return '1'
    else:
        return ''


def cscomp_simp(group,df,rois,phases=None):
    df = df.loc[(group,rois,phases)].reset_index()
    df.encode_phase = pd.Categorical(df.encode_phase.values,categories=df.encode_phase.unique(),ordered=True)

    pvals = pd.DataFrame(columns=['p'],index=pd.MultiIndex.from_product([rois,phases],names = ['roi','encode_phase']))
    df = df.set_index(['roi','encode_phase'])
    for roi in rois:
        for phase in phases:
            pvals.loc[(roi,phase),'p'] = pg.ttest(df.loc[(roi,phase),'rsa'].values,0)['p-val'].values[0]
    pvals['corrp'] = pg.multicomp(list(pvals.p.values),method='fdr_bh')[1]
    pvals['sig'] = pvals.corrp.apply(pconvert)
    df = df.reset_index()

    if len(phases) == 4: phase_pal = sns.color_palette(['darkgrey','darkmagenta','lightgreen','seagreen'],desat=1)
    elif len(phases) == 3:phase_pal = sns.color_palette(['darkgrey','darkmagenta','seagreen'],desat=1)

    fig, ax = plt.subplots(1,len(rois),sharey=True,figsize=(9,3.5))
    for i, roi in enumerate(rois):
        if len(rois) == 1: Ax = ax
        else: Ax = ax[i]
        dat = df[df.roi == roi]
        sns.barplot(data=dat,x='encode_phase',y='rsa',
                    palette=phase_pal,ax=Ax,errcolor='black')
        X = Ax.get_xticks()

        # Ax.set_ylim(ymin,ymax)            
        Ax.yaxis.set_major_locator(MultipleLocator(.1))
        Ax.yaxis.set_minor_locator(MultipleLocator(.05))
        
        sns.despine(ax=Ax,trim=True)

        for i, x in enumerate(X): Ax.annotate(pvals.loc[roi,'sig'].values[i], xy=(x-.05,Ax.get_ylim()[1]-.1))

    print(yval)
    if len(rois) > 1:
        ax[0].set_ylabel('∆ fisher z(r)') if yval=='rsa' else ax[0].set_ylabel('CS+ > CS- activity')
        ax[1].set_ylabel('')
    else:
        ax.set_ylabel('∆ fisher z(r)')
        ax.set_ylabel('')
    # plt.savefig('./plots/group_comp_%s.eps'%(split),format='eps')
    # plt.close()

def cscomp(group,df,rois,statsdf,n_boot=10000,phases=None,yval='rsa'):
    df = df.loc[group]
    out = {}
    for val in ['ci','dist']: out[val] = {}
    for roi in rois:
        out['ci'][roi] = {}
        out['dist'][roi] = {}
        for phase in phases:
            out['ci'][roi][phase] = {}
            out['dist'][roi][phase] = {}

            out['ci'][roi][phase], out['dist'][roi][phase] = pg.compute_bootci(df.loc[(roi,phase),yval].values,
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
    
    ci.roi = pd.Categorical(ci.roi,rois,ordered=True)
    ci.encode_phase = pd.Categorical(ci.encode_phase,phases,ordered=True)
    ci = ci.sort_values(by='encode_phase')
    ci = ci.set_index('roi').sort_index()


    dist.roi = pd.Categorical(dist.roi,rois,ordered=True)
    dist.encode_phase = pd.Categorical(dist.encode_phase,phases,ordered=True)
    ci = ci.sort_values(by='encode_phase')
    dist = dist.set_index('roi').sort_index()

    pvals = statsdf.loc[group]
    pvals['sig'] = pvals.p.apply(pconvert)

    ymax = dist.dist.max() + .05
    ymin = dist.dist.min() - .05

    # phase_pal = sns.color_palette(['black','darkmagenta','lightgreen','seagreen'],desat=.75)
    fig, ax = plt.subplots(1,len(rois),sharey=False,figsize=(len(rois)*5,8))
    for i, roi in enumerate(rois):
        if len(rois) == 1: Ax = ax
        else: Ax = ax[i]
        dat = df.loc[roi].reset_index()
        

        pal_colors = ['darkgrey','darkmagenta','seagreen']
        stable_pal = ['black','darkmagenta','seagreen']

        # point_colors = ['white','white','white']
        for pi, pval in enumerate(pvals.loc[roi,'sig'].values):
            if pval == '': pal_colors[pi] = 'white'

        phase_pal = sns.color_palette(pal_colors,desat=.75)

        #dist violins
        sns.violinplot(data=dist.loc[roi],x='encode_phase',y='dist',
                        inner=None,ax=Ax,scale='count',palette=phase_pal,saturation=1)
        
        for viol in Ax.collections: viol.set_edgecolor('black')
        
        # sns.barplot(data=dat,x='encode_phase',y='rsa',hue='encode_phase',palette=phase_pal,ax=Ax,
                    # order=phases,seed=42,errcolor='black')
        # change_width(Ax,.75)
        # sns.pointplot(data=dat,x='encode_phase',y='rsa',hue='encode_phase',palette=phase_pal,ax=Ax,
                        # order=phases,seed=42)
        X = Ax.get_xticks()

        # sns.pointplot(data=dat,x='encode_phase',y='rsa',ci=None,units='subject',order=phases,ax=Ax)
        # sdat = dat.set_index(['subject','encode_phase'])
        # for sub in dat.subject.unique():
            # Ax.plot(X,sdat.loc[(sub,phases),'rsa'].values,color='grey',lw=1,alpha=.5)


        lower = ci.loc[roi,'ci'].apply(lambda x: x[0])
        upper = ci.loc[roi,'ci'].apply(lambda x: x[1])
        Y = ci.loc[roi,'point']


        Ax.vlines(X,lower,upper,linewidth=3,color='black').set_capstyle('round')
        Ax.scatter(X,Y,s=50,color='black')
        Ax.hlines(0,Ax.get_xlim()[0],Ax.get_xlim()[1],color='black',linestyle=':',linewidth=3)


        Ax.set_xticklabels(['Pre acq.','Acquisition','Extinction'],rotation=45,ha='center',fontsize=20)
        for labeli, t in enumerate(Ax.xaxis.get_ticklabels()): t.set_color(stable_pal[labeli])
        Ax.set_xlabel('')


        Ax.set_title(roi,fontsize=30)
        # Ax.legend_.remove()

        Ax.set_ylim(ymin,ymax)            
        # Ax.yaxis.set_major_locator(MultipleLocator(.1))
        # Ax.yaxis.set_minor_locator(MultipleLocator(.05))
        

        # for j, x in enumerate(X): Ax.annotate(pvals.loc[roi,'sig'].values[j], xy=(x-.05,Ax.get_ylim()[1]-.1))


    if len(rois) > 1:
        ax[0].set_ylabel('∆ fisher z(r)') if yval=='rsa' else ax[0].set_ylabel('CS+ > CS- univariate activity')
        ax[1].set_ylabel('')
        sns.despine(ax=ax[1],left=True)
        # ax[1].set_ylim(ax[0].get_ylim()[0],ax[0].get_ylim()[1])
        ax[1].set_yticks([])


    else:
        ax.set_ylabel('∆ fisher z(r)')
        ax.set_ylabel('')

    plt.tight_layout()

def split_level(df,group,rois=['rSMA','sgACC'],phases=None,split=None):
    if split == 'level':
        splits = ['item','set']
        # pal = sns.palettes.color_palette('Set2',n_colors=2)
        pal = ['blue','lightblue']

    elif split == 'memory_phase':
        splits = ['encoding','retrieval']
        # pal = sns.palettes.color_palette('Set2',n_colors=4)[2:]
        pal = ['firebrick','salmon']

    df = df.loc[group]
    out = {}

    for val in ['ci','dist']: out[val] = {}
    for roi in rois:
        out['ci'][roi] = {}
        out['dist'][roi] = {}
        for phase in phases:
            out['ci'][roi][phase] = {}
            out['dist'][roi][phase] = {}
            for level in splits:
                out['ci'][roi][phase][level] = {}
                out['dist'][roi][phase][level] = {}

                out['ci'][roi][phase][level], out['dist'][roi][phase][level] = pg.compute_bootci(df.loc[(roi,phase,level),'rsa'].values,
                                                                func='mean',n_boot=10000,return_dist=True,
                                                                method='cper',decimals=3,seed=42)
    
    ci = pd.DataFrame.from_dict(out['ci'],orient='index').reset_index().rename(columns={'index':'roi'})
    ci = ci.melt(id_vars='roi',var_name='encode_phase',value_name='ci').set_index(['roi','encode_phase']).sort_index()
    ci = ci.ci.apply(pd.Series).reset_index()
    ci = ci.melt(id_vars=['roi','encode_phase'],var_name='level',value_name='ci').set_index(['roi','encode_phase','level'])
    
    dist = pd.DataFrame.from_dict(out['dist'],orient='index').reset_index().rename(columns={'index':'roi'})
    dist = dist.melt(id_vars='roi',var_name='encode_phase',value_name='dist').set_index(['roi','encode_phase'])
    dist = dist.dist.apply(pd.Series).reset_index()
    dist = dist.melt(id_vars=['roi','encode_phase'],var_name='level',value_name='dist').set_index(['roi','encode_phase','level'])
    dist = dist.dist.apply(pd.Series).stack().reset_index(-1,drop=True)
    dist = dist.reset_index().rename(columns={0:'dist'})
    
    ci['point'] = dist.groupby(['roi','encode_phase','level']).mean()
    ci = ci.reset_index()
    
    ci.roi = pd.Categorical(ci.roi,rois,ordered=True)
    ci.encode_phase = pd.Categorical(ci.encode_phase,phases,ordered=True)
    ci.level = pd.Categorical(ci.level,splits,ordered=True)
    ci = ci.sort_values(by=['encode_phase','level'])
    ci = ci.set_index('roi').sort_index()


    dist.roi = pd.Categorical(dist.roi,rois,ordered=True)
    dist.encode_phase = pd.Categorical(dist.encode_phase,phases,ordered=True)
    dist.level = pd.Categorical(dist.level,splits,ordered=True)
    dist = dist.sort_values(by=['encode_phase','level'])
    dist = dist.set_index('roi').sort_index()
    
    #collect 2 tailed p-values
    # pvals = dist.reset_index().groupby(['roi','encode_phase','level']
                # ).apply(lambda x : (1-np.mean(x>0))/2
                # ).reset_index().set_index('roi').rename(columns={'dist':'p'})
    # for roi in rois: pvals.loc[roi,'corrp'] = pg.multicomp(pvals.loc[roi,'p'].values,method='fdr_bh')[1]
    pvals = pd.DataFrame(columns=['p'],index=pd.MultiIndex.from_product([rois,phases,splits],
                                                            names = ['roi','encode_phase','level']))
    for roi in rois:
        for phase in phases:
            for s in splits:
                pvals.loc[(roi,phase,s),'p'] = pg.ttest(df.loc[(roi,phase,s),'rsa'].values,0)['p-val'].values[0]
    # for roi in rois: pvals.loc[roi,'corrp'] = pg.multicomp(list(pvals.loc[roi,'p'].values),method='fdr_bh')[1]
    pvals['corrp'] = pg.multicomp(list(pvals.p.values),method='fdr_bh')[1]
    pvals['sig'] = pvals.corrp.apply(pconvert)

    phase_pal = sns.color_palette(['black','darkmagenta','lightgreen','seagreen'],desat=.75)

    ymax = dist.dist.max() + .05
    ymin = dist.dist.min() - .05
    fig, ax = plt.subplots(1,2,sharey=True,figsize=(8,3))
    for i, roi in enumerate(rois):
        dat = df.loc[roi].reset_index()
        Ax = ax[i]
        
        '''
        #dist violins
        sns.violinplot(data=dist.loc[roi],x='encode_phase',y='dist',cut=0,
                        inner=None,ax=Ax,scale='count',hue='level',split=True,palette=pal)#,palette=['white','white'])
        '''
        sns.barplot(data=dat,x='encode_phase',y='rsa',hue=split,palette=pal,ax=Ax,
                    order=phases,seed=42,errcolor='black')
        change_width(Ax,.35)
        '''
        sns.violinplot(data=dat,x='encode_phase',y='rsa',hue=split,
                        split=True,palette=pal,cut=0,ax=Ax,inner=None)
        '''

        '''
        #jointplot
        sns.pointplot(data=dat,x='encode_phase',y='rsa',hue=split,
                     dodge=.3,n_boot=1000,seed=42,markers=['o','D'],ax=Ax,palette=pal,
                     join=False)
        '''

        lower = ci.loc[roi,'ci'].apply(lambda x: x[0])
        upper = ci.loc[roi,'ci'].apply(lambda x: x[1])
        Y = ci.loc[roi,'point']
        X_together = [[x-.1,x+.1] for x in Ax.get_xticks()]
        X = [x for t in X_together for x in t]
        
        # Ax.vlines(X,lower,upper,linewidth=3,color='black').set_capstyle('butt')
        # Ax.scatter(X,lower,s=240,marker='_',color='black',linewidths=5)
        # Ax.scatter(X,upper,s=240,marker='_',color='black',linewidths=5)
        # Ax.scatter(X,Y,s=480,marker='_',color='black',edgecolors='black') #set level



        # #stripplots
        # sns.stripplot(data=dat,x='encode_phase',y='rsa',hue=split,palette=pal,
        #             ax=Ax,dodge=.3,order=phases,alpha=.3,edgecolor='black',linewidth=1)



        '''
        #chicken scratch
        dat = dat.set_index(['encode_phase',split])
        for p, phase in enumerate(phases):
            left = dat.loc[(phase,splits[0]),'rsa']
            right = dat.loc[(phase,splits[1]),'rsa']
            for xy in range(left.shape[0]):
                xvals = X_together[p]
                xvals[0] -= .01
                xvals[1] += .01
                Ax.plot(X_together[p],[left[xy],right[xy]],color='grey',lw=1,alpha=.5)
        '''

        Ax.hlines(0,Ax.get_xlim()[0],Ax.get_xlim()[1],color='grey',linestyle='--',linewidth=3)
        Ax.set_xticklabels('',rotation=45)
        Ax.set_title(roi+'_'+group+'_'+split)
        
        Ax.legend_.remove()

        Ax.set_ylim(ymin,ymax)            
        Ax.yaxis.set_major_locator(MultipleLocator(.2))
        Ax.yaxis.set_minor_locator(MultipleLocator(.1))
        
        sns.despine(ax=Ax,trim=True)

        for i, x in enumerate(X): Ax.annotate(pvals.loc[roi,'sig'].values[i], xy=(x-.05,Ax.get_ylim()[1]-.1))

    
    if len(rois) > 1:
        ax[0].set_ylabel('∆ fisher z(r)')
        ax[1].set_ylabel('')
    else:
        ax.set_ylabel('∆ fisher z(r)')
        ax.set_ylabel('')
    plt.savefig('./plots/%s_%s.eps'%(group,split),format='eps')
    # plt.close()

def group_comp_simp(df,phases,rois=['rACC','sgACC']):
    comp = pd.DataFrame(index=pd.MultiIndex.from_product([rois,phases],names=['roi','encode_phase']))
    phase_pal = sns.color_palette(['black','darkmagenta','lightgreen','seagreen'],desat=1)

    for roi in rois:
        for phase in phases:
                
                c = df.loc[('healthy',roi,phase),'rsa'].values
                cci, cdist = pg.compute_bootci(c,func='mean',n_boot=1000,return_dist=True,method='cper',decimals=3,seed=42)
                
                p = df.loc[('ptsd',roi,phase),'rsa'].values
                cci, pdist = pg.compute_bootci(p,func='mean',n_boot=1000,return_dist=True,method='cper',decimals=3,seed=42)
                
                dif = cdist - pdist

                comp.loc[(roi,phase),'point'] = dif.mean()
                
                # tres = pg.ttest(c,p,paired=False,correction=False)
                comp.loc[(roi,phase),'ci_l'] = np.percentile(dif,5)#tres['CI95%'][0][0]
                comp.loc[(roi,phase),'ci_u'] = np.percentile(dif,100)#tres['CI95%'][0][1]
                comp.loc[(roi,phase),'p'] = 1-np.mean(dif > 0)#tres['p-val'][0]
        comp.loc[roi,'corrp'] = pg.multicomp(list(comp.loc[roi,'p'].values),method='fdr_bh')[1]
    
    comp['sig'] = comp.corrp.apply(pconvert)
    
    fig, ax = plt.subplots(1,2,sharey=True,sharex=True,figsize=(10.5,4.5))
    for i, roi in enumerate(rois):
        dat = comp.loc[roi].reset_index()
        Ax = ax[i]
        
        sns.barplot(data=dat,y='encode_phase',x='point',palette=phase_pal,ax=Ax,
                    seed=42,errcolor='black',orient='h')
        # X_together = [[x-.2,x+.2] for x in Ax.get_xticks()]
        # X = [x for t in X_together for x in t]
        X = Ax.get_yticks()
        Ax.hlines(X,dat.ci_l.values,dat.ci_u.values,linewidth=5,color='black').set_capstyle('round')
        Ax.vlines(0,Ax.get_ylim()[0],Ax.get_ylim()[1],color='grey',linestyle='--',linewidth=3)

        # Ax.legend_.remove()

        # Ax.set_ylim(ymin,ymax)            
        Ax.xaxis.set_major_locator(MultipleLocator(.1))
        Ax.xaxis.set_minor_locator(MultipleLocator(.05))
        
        sns.despine(ax=Ax,trim=True)

        for i, x in enumerate(X): Ax.annotate(comp.loc[roi,'sig'].values[i], xy=(x-.05,Ax.get_ylim()[1]-.1))

    

    ax[0].set_xlabel('∆ fisher z(r)')
    ax[1].set_xlabel('∆ fisher z(r)')


def group_comp(df,phases,split):
    if split == 'level':
        splits = ['item','set']
        # pal = sns.palettes.color_palette('Set2',n_colors=2)
        pal = ['blue','lightblue']

    elif split == 'memory_phase':
        splits = ['encoding','retrieval']
        # pal = sns.palettes.color_palette('Set2',n_colors=4)[2:]
        pal = ['firebrick','salmon']
    rois = ['rSMA','sgACC']
    comp = pd.DataFrame(index=pd.MultiIndex.from_product([rois,phases,splits],names=['roi','encode_phase',split]))
    for roi in rois:
        for phase in phases:
            for level in splits:
                c = df.loc[('control',roi,phase,level),'rsa'].values
                p = df.loc[('ptsd',roi,phase,level),'rsa'].values
                comp.loc[(roi,phase,level),'point'] = c.mean() - p.mean()
                tres = pg.ttest(c,p,paired=False,correction=False)
                comp.loc[(roi,phase,level),'ci_l'] = tres['CI95%'][0][0]
                comp.loc[(roi,phase,level),'ci_u'] = tres['CI95%'][0][1]
                comp.loc[(roi,phase,level),'p'] = tres['p-val'][0]
    for roi in rois: comp.loc[roi,'corrp'] = pg.multicomp(list(comp.loc[roi,'p'].values),method='fdr_bh')[1]
    comp['sig'] = comp.corrp.apply(pconvert)

    ymax = comp.ci_u.max() + .05
    ymin = comp.ci_l.min() - .05
    fig, ax = plt.subplots(1,2,sharey=True,figsize=(24,12))
    for i, roi in enumerate(rois):
        dat = comp.loc[roi].reset_index()
        Ax = ax[i]
        
        sns.barplot(data=dat,x='encode_phase',y='point',hue=split,palette=pal,ax=Ax,
                    order=phases,seed=42,errcolor='black')
        X_together = [[x-.2,x+.2] for x in Ax.get_xticks()]
        X = [x for t in X_together for x in t]
        Ax.vlines(X,dat.ci_l.values,dat.ci_u.values,linewidth=5,color='black').set_capstyle('round')
        Ax.hlines(0,Ax.get_xlim()[0],Ax.get_xlim()[1],color='grey',linestyle='--',linewidth=3)

        Ax.legend_.remove()

        Ax.set_ylim(ymin,ymax)            
        Ax.yaxis.set_major_locator(MultipleLocator(.1))
        Ax.yaxis.set_minor_locator(MultipleLocator(.05))
        
        sns.despine(ax=Ax,trim=True)

        for i, x in enumerate(X): Ax.annotate(comp.loc[roi,'sig'].values[i], xy=(x-.05,Ax.get_ylim()[1]-.1))

    
    if len(rois) > 1:
        ax[0].set_ylabel('∆ fisher z(r)')
        ax[1].set_ylabel('')
    else:
        ax.set_ylabel('∆ fisher z(r)')
        ax.set_ylabel('')
 
    plt.savefig('./plots/group_comp_%s.eps'%(split),format='eps')
    # plt.close()


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

    pvals = pd.DataFrame(columns=['p'],index=pd.MultiIndex.from_product([rois,phases],names = ['roi','encode_phase']))
    pvals = pd.DataFrame(columns=['p'],index=pd.MultiIndex.from_product([rois,phases,mems],names = ['roi','encode_phase','mem']))
    for roi in rois:
        for phase in phases:
            for mem in mems:
                pvals.loc[(roi,phase,mem),'p'] = pg.ttest(df.loc[(roi,phase,mem),'rsa'].values,0)['p-val'].values[0]
            # pvals.loc[(roi,phase),'p'] = pg.ttest(df.loc[(roi,phase,'hit'),'rsa'], df.loc[(roi,phase,'miss'),'rsa'])['p-val'].values[0]
    pvals['corrp'] = pg.multicomp(list(pvals.p.values),method='fdr_bh')[1]
    pvals['sig'] = pvals.corrp.apply(pconvert)

    print(pvals)


    # hit_pal = sns.color_palette(['black','darkmagenta','lightgreen','seagreen'],desat=.75)
    # miss_pal = sns.color_palette(['midnightblue','darkmagenta','lightgreen','seagreen'],desat=.25)
    # phase_pal = np.empty(8,dtype=object)
    # phase_pal[0::2] = hit_pal
    # phase_pal[1::2] = miss_pal
    mem_pal = ['darkblue','lightblue']

    fig, ax = plt.subplots(1,len(rois),sharey=True,figsize=(6.5,5))
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

    # for i, x in enumerate(X): Ax.annotate(pvals.loc[roi,'sig'].values[i], xy=(x-.05,Ax.get_ylim()[1]-.1))
