from fg_config import *
from paper_graphics_style import *
import pingouin as pg
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D




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

Wpal = np.repeat('white',6)

pfc_rois = {'dACC':{'stars':['','','***','***','','**']},
           'vmPFC':{'stars':['','*','*','','***','']}}

''' using ROI as columns'''
fig, ax = plt.subplots(1,2,figsize=mm2inch(200,90),sharey=False)
for r, roi in enumerate(pfc_rois):
    stats = diff[f'{roi}_ers'].groupby(['group','phase']).apply(pg.compute_bootci,None,'mean',n_boot=10000,decimals=4,seed=42,return_dist=True).reset_index()
    stats.phase = pd.Categorical(stats.phase, categories=phase3, ordered=True)
    stats = stats.set_index(['group','phase']).sort_index()

    sig = stats[f'{roi}_ers'].apply(lambda x: 1 if np.sign(x[0][0]) == np.sign(x[0][1]) else 0).reset_index()[f'{roi}_ers'].values
    # for i, val in enumerate(sig):
    #     if not val: sigpal[i] = (1.0,1.0,1.0)

    stats = stats[f'{roi}_ers'].apply(pd.Series).drop(columns=0)
    stats = stats[1].apply(pd.Series).stack().reset_index().drop(columns='level_2').rename(columns={0:'est'})

    sns.violinplot(data=stats,x='phase',y='est',hue='group',order=phase3,cut=0,saturation=1,
                    inner='mean_ci',palette=Wpal,ax=ax[r],scale='count',zorder=1,edgecolor='white')

    sns.violinplot(data=stats,x='phase',y='est',hue='group',order=phase3,cut=0,saturation=1,
                    inner='mean_ci',palette=double_pal,ax=ax[r],scale='count',zorder=10,edgecolor='white')
    for v, viol in enumerate(ax[r].collections[:6]):viol.set_edgecolor(double_pal[v])
    for v, viol in enumerate(ax[r].collections[6:]):
        viol.set_edgecolor(double_pal[v])
        viol.set_facecolor(double_pal[v])
    for viol in ax[r].collections[7::2]: viol.set_alpha(.2)
    ax[r].hlines(0,*ax[r].get_xlim(),linestyle=':',color='black',zorder=0)
    ax[r].set_ylim(-.275,.45)
    ax[r].yaxis.set_major_locator(MultipleLocator(.2))   
    ax[r].legend_.remove()
    sns.despine(left=True,ax=ax[r])
    ax[r].set_xlabel('')
    ax[r].set_xticklabels(['Pre-\nconditioning','Conditioning','Extinction'],ha='center')
    if r == 0:
        ax[r].set_ylabel('Encoding-retrieval similarity')
    else:
        ax[r].set_ylabel('')
        ax[r].set_yticks([])

    ax[r].set_title(roi)

    top = stats.groupby(['phase','group']).est.apply(np.max).values
    bottom = stats.groupby(['phase','group']).est.apply(np.min).values
    if roi == 'vmPFC': top[1] = bottom[1] - .05
    xwhere = [(i-.2,i+.2) for i in range(3)]
    xwhere = [i for j in xwhere for i in j]
    for star, x, y in zip(pfc_rois[roi]['stars'],xwhere,top): ax[r].text(x,y,star,ha='center',va='bottom',fontsize=10)

legend_elements = [Patch(facecolor='black',edgecolor='black',label='Healthy'),
                   Patch(facecolor='lightgrey',edgecolor='black',label='PTSS')]
fig.text(.05,.85,'Shown as CS+ - CS-',ha='left',va='bottom',fontsize=8)
fig.legend(handles=legend_elements,ncol=1,frameon=False,loc='upper right',bbox_to_anchor=(1,.9))
plt.tight_layout()
# plt.tight_layout(rect=(0,0,.9,1))

# pfc_groups = {'healthy':{'stars':['','','***','*','','***']},
#                  'ptsd':{'stars':['','*','***','','**','']}}

# ''' using ROI as columns'''
# df = pd.read_csv('pfc_ers_cleaned_lmm.csv')
# df.roi = df.roi.apply(lambda x: x[:-4])
# df = df.set_index(['group','phase','condition','subject','stimulus','roi']).sort_index()
# df = df.reset_index().drop(columns='stimulus').groupby(['condition','group','phase','roi','subject']).mean()
# diff = df.loc['CS+'] - df.loc['CS-']

# fig, ax = plt.subplots(1,2,figsize=mm2inch(200,90),sharey=False)
# for r, group in enumerate(pfc_groups):
#     stats = diff.loc[group].reset_index().groupby(['phase','roi'])['ers'].apply(pg.compute_bootci,None,'mean',n_boot=10000,decimals=4,seed=42,return_dist=True).reset_index()
#     stats.phase = pd.Categorical(stats.phase, categories=phase3, ordered=True)
#     stats = stats.set_index(['phase','roi']).sort_index()

#     stats = stats['ers'].apply(pd.Series).drop(columns=0)
#     stats = stats[1].apply(pd.Series).stack().reset_index().drop(columns='level_2').rename(columns={0:'est'})

#     sns.violinplot(data=stats,x='phase',y='est',hue='group',order=phase3,cut=0,saturation=1,
#                     inner='mean_ci',palette=Wpal,ax=ax[r],scale='count',zorder=1,edgecolor='white')

#     sns.violinplot(data=stats,x='phase',y='est',hue='group',order=phase3,cut=0,saturation=1,
#                     inner='mean_ci',palette=double_pal,ax=ax[r],scale='count',zorder=10,edgecolor='white')
  


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


df = pd.read_csv('pfc_ers_cleaned.csv').set_index(['condition','group','phase','roi','subject']).sort_index()
df = df.loc['CS+'] - df.loc['CS-']

fig, ax = plt.subplots()
sns.pointplot(data=df.loc['healthy'].reset_index(),x='phase',y='rsa',hue='roi',order=phase3,ax=ax)
sns.stripplot(data=df.loc['healthy'].reset_index(),x='phase',y='rsa',hue='roi',order=phase3,ax=ax)

fig, ax = plt.subplots()
sns.violinplot(data=df.loc['healthy'].reset_index(),x='phase',y='rsa',hue='roi',split=True,cut=0,order=phase3,ax=ax)







