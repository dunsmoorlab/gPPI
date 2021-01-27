from fg_config import *
import pingouin as pg
from matplotlib.ticker import MultipleLocator
sns.set_context('paper')
rcParams['savefig.dpi'] = 900
def mm2inch(*tupl):
    inch = 25.4
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


df = pd.read_csv('all_data_lmm.csv').set_index(['group','phase','condition','subject','stimulus']).sort_index()
df = df.loc[(['healthy','ptsd'],['acquisition','extinction']),].reset_index().drop(columns='stimulus')
df = df.groupby(['group','condition','subject']).mean().reset_index()
# df['phase_con'] = df.phase + '_' + df.condition

vals = {'CS+':{'est': 0.1387,
               'ci': [0.0692, 0.208],
               'int': -0.1297,
               'int_ci': [-0.2035,-0.0558]},
        'CS-':{'est': 0.0479,
               'ci': [-0.0214, 0.117],
               'int': -0.0167,
               'int_ci': [-0.0901,0.0566]}}

fig, ax = plt.subplots()
sns.scatterplot(data=df,x='ev',y='pfc_diff_ers',hue='condition',palette=cpal,ax=ax)
ax.legend(loc='center left',bbox_to_anchor=(1,.5))


# ax.scatter(df.groupby('subject').mean()['ev'],df.groupby('subject').mean()['pfc_diff_ers'])
for c, cs in enumerate(['CS+','CS-']):
    ax.plot(np.linspace(0,1,100),(np.linspace(0,1,100)*vals[cs]['est'])+vals[cs]['int'],color=cpal[c])
    #ax.fill_between(np.linspace(0,1,100),(np.linspace(0,1,100)*vals[cs]['est'])+vals[cs]['int_ci'][0],(np.linspace(0,1,100)*vals[cs]['est'])+vals[cs]['int_ci'][1],color=cpal[c],alpha=.5)
# ax.set_ylabel('pfc_diff_ers');ax.set_xlabel('ev')


'''PFC ERS'''
df = pd.read_csv('all_data_lmm.csv').set_index(['group','phase','condition','subject','stimulus']).sort_index()
df = df.loc[(['healthy','ptsd'],['acquisition','baseline','extinction']),].reset_index().drop(columns='stimulus')
df = df.groupby(['condition','group','phase','subject']).mean()
diff = df.loc['CS+'] - df.loc['CS-']

Wpal = ['white','white','white']

for roi in ['vmPFC','dACC']:
    stats = diff[f'{roi}_ers'].groupby(['group','phase']).apply(pg.compute_bootci,None,'mean',n_boot=10000,decimals=4,seed=42,return_dist=True).reset_index()
    stats.phase = pd.Categorical(stats.phase, categories=phase3, ordered=True)
    stats = stats.set_index(['group','phase']).sort_index()

    sig = stats[f'{roi}_ers'].apply(lambda x: 1 if np.sign(x[0][0]) == np.sign(x[0][1]) else 0).reset_index()[f'{roi}_ers'].values
    sigpal = phase_pal+phase_pal
    # for i, val in enumerate(sig):
    #     if not val: sigpal[i] = (1.0,1.0,1.0)

    stats = stats[f'{roi}_ers'].apply(pd.Series).drop(columns=0)
    stats = stats[1].apply(pd.Series).stack().reset_index().drop(columns='level_2').rename(columns={0:'est'})

    fig, ax = plt.subplots(figsize=mm2inch(90,80))
    sns.violinplot(data=stats,x='phase',y='est',hue='group',order=phase3,cut=0,saturation=1,
                    inner='mean_ci',palette=cpal,ax=ax,scale='count',zorder=1,edgecolor='white')

    sns.violinplot(data=stats,x='phase',y='est',hue='group',order=phase3,cut=0,saturation=1,
                    inner='mean_ci',palette=cpal,ax=ax,scale='count',zorder=10,edgecolor='white')
    # for v, viol in enumerate(ax.collections[6:]):
    #     viol.set_edgecolor((phase_pal+phase_pal)[v])
    #     viol.set_facecolor(sigpal[v])
    ax.hlines(0,*ax.get_xlim(),linestyle=':',color='black',zorder=0)
    ax.set_ylim(-.25,.45)
    ax.yaxis.set_major_locator(MultipleLocator(.2))   
    ax.legend_.remove()
    sns.despine(left=True,ax=ax)
    ax.set_xlabel('')
    ax.set_xticklabels(['Healthy','PTSS'])
    ax.set_ylabel('Encoding-retrieval similarity')
    ax.set_title(roi)
    plt.tight_layout()

df = pd.read_csv('all_data_lmm.csv').set_index(['group','phase','condition','subject','stimulus']).sort_index()
df = df.loc[(['healthy','ptsd'],['acquisition','baseline','extinction']),].reset_index().drop(columns='stimulus')
df = df.groupby(['group','phase','subject']).mean()

hc_rois = {'hc_tail':'Tail (Posterior)', 'hc_head':'Head (Anterior)'}
fig, ax = plt.subplots(1,2,figsize=mm2inch(120,80),sharey=True)
for r, roi in enumerate(hc_rois):
    stats = df[f'{roi}_ers'].groupby(['group','phase']).apply(pg.compute_bootci,None,'mean',n_boot=10000,decimals=4,seed=42,return_dist=True).reset_index()
    stats.phase = pd.Categorical(stats.phase, categories=phase3, ordered=True)
    stats = stats.set_index(['group','phase']).sort_index()

    sig = stats[f'{roi}_ers'].apply(lambda x: 1 if np.sign(x[0][0]) == np.sign(x[0][1]) else 0).reset_index()[f'{roi}_ers'].values
    # sigpal = phase_pal+phase_pal
    sigalpha = [1,1,1,1,1,1]
    for i, val in enumerate(sig):
        # if not val: sigpal[i] = (1.0,1.0,1.0)
        if val == 0: sigalpha[i] = .3

    stats = stats[f'{roi}_ers'].apply(pd.Series).drop(columns=0)
    stats = stats[1].apply(pd.Series).stack().reset_index().drop(columns='level_2').rename(columns={0:'est'})

    sns.violinplot(data=stats,x='group',y='est',hue='phase',cut=0,saturation=1,
                    inner='mean_ci',palette=Wpal,ax=ax[r],scale='count',zorder=1,edgecolor='white')

    sns.violinplot(data=stats,x='group',y='est',hue='phase',order=phase3,cut=0,saturation=1,
                    inner='mean_ci',palette=phase_pal,ax=ax[r],scale='count',zorder=10,edgecolor='white')
    for v, viol in enumerate(ax[r].collections[6:]):
        viol.set_edgecolor((phase_pal+phase_pal)[v])
        viol.set_alpha(sigalpha[v])

    ax[r].hlines(0,*ax[r].get_xlim(),linestyle=':',color='black',zorder=0)
    ax[r].yaxis.set_major_locator(MultipleLocator(.05))   
    ax[r].legend_.remove()
    sns.despine(left=True,ax=ax[r])
    ax[r].set_xlabel('')
    ax[r].set_xticklabels(['Healthy','PTSS'])
    ax[r].set_ylabel('Encoding-retrieval similarity') if r == 0 else ax[r].set_ylabel('')
    ax[r].set_title(hc_rois[roi])
plt.tight_layout()   