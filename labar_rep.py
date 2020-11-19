'''prettier labar replication graphs'''
df = pd.read_csv('labar_replication_data.csv')
stats = pd.read_csv('labar_replication_stats.csv')

df.roi = df.roi.apply(pfc_rename)
df['rds'] = df.rsa.apply(lambda x: 1-np.tanh(x))
df = df.set_index(['roi','group']).sort_index()


pfc_rois = ['dACC','vmPFC']
amyg = ['amyg_cem','amyg_bla']
hpc = ['hc_head','hc_body','hc_tail']
df = df.loc[rois].reset_index().set_index(['memory_phase','group','roi','condition'])
group = 'healthy'

roi_dct = {'pfc':pfc_rois,
           'amyg':amyg,
           'hpc':hpc}

pals = [sns.xkcd_palette(['windows blue','amber']), sns.color_palette('rocket',2), sns.color_palette('mako',3)]
rcParams['axes.titlepad'] = 0
for mem_phase in ['encoding','retrieval']:
    print(mem_phase)
    for group in groups:
        for r, _rois_ in enumerate(roi_dct):            
            fig, ax = plt.subplots(1,2,sharey=True,figsize=(8,4))
            for i, con in enumerate(['CS+','CS-']):
                Ax = ax[i]
                dat = df.loc[mem_phase,group,roi_dct[_rois_],con].reset_index()

                sns.pointplot(data=dat,x='el',y='rds',hue='roi',palette=pals[r],dodge=True,ax=Ax)
                
                Ax.set_title(con)#,fontsize=10)
                Ax.set_xticklabels(['Early\nextinction','Late\nextinction'])#,fontsize=8)
                Ax.set_xlabel('')

            # ax[0].set_ylim(.65,1.1)
            ax[0].set_ylabel('Dissimilarity with\nlate acquisition')#,fontsize=10)
            ax[1].set_ylabel('')
            ax[1].legend_.remove()
            plt.suptitle(f'Group = {group}\nmemory phase = {memory_phase}',ha='center')#, fontsize=10, y=1)
            # plt.tight_layout()
            
            # plt.savefig(f'plots/labar_replication/{group}_{memory_phase}_{_rois_}.png')


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


    Ax.set_xticklabels(['Baseline','Acquisition','Extinction'],rotation=45,ha='center',fontsize=20)
    for labeli, t in enumerate(Ax.xaxis.get_ticklabels()): t.set_color(stable_pal[labeli])
    Ax.set_xlabel('')


    Ax.set_title(roi,fontsize=30)
    # Ax.legend_.remove()

    Ax.set_ylim(ymin,ymax)
