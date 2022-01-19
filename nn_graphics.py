from fg_config import *
from paper_graphics_style import *
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from matplotlib import cm, gridspec
from itertools import product
import matplotlib as mpl
import pingouin as pg

'''PFC barplots with subject dots'''
df = pd.read_csv('pfc_csdif_stats.csv')
df.roi = df.roi.apply(lambda x: x[:-4])
df = df[df.phase != 'baseline']
df = df.set_index(['group','phase','roi']).sort_index()

hpc_groups = {'healthy':{'name':'Healthy','stars':['***','*','','***'],'subs':sub_args,'marker':'D','line':'-','c':'k'},
                 'ptsd':{'name':'PTSS','stars':['***','','**',''],'subs':p_sub_args,'marker':'s','line':'--','c':'grey'}}

xvals = [(i-.2,i+.2) for i in range(2)]
xvals = [i for j in xvals for i in j]
pal2 = sns.color_palette(['red','dodgerblue'],desat=.8)
pal26 = [pal2[0],pal2[0],pal2[1],pal2[1]]

subdat = pd.read_csv('all_data.csv')
subdf = pd.read_csv('pfc_ers_cleaned.csv').set_index(['condition','phase','group','roi','subject']).sort_index()
subdf = (subdf.loc['CS+'] - subdf.loc['CS-']).drop('baseline')

subxvals = [(i+.1,i-.1) for i in range(2)]
subxvals = [i for j in subxvals for i in j]

fig, ax = plt.subplots(2,1,figsize=(mm2inch(86,150)),sharey=True)
for r, group in enumerate(hpc_groups):

    #subject data points
    for bar, (phase, roi) in enumerate(product(['acquisition','extinction'],['vmPFC','dACC'])):
        # edge = pal26[bar] if not bar & 1 else 'k'
        # face = 'white' if not bar & 1 else pal26[bar]
        edge=pal26[bar];face='white';
        ax[r].scatter(np.repeat(subxvals[bar],24),subdf.loc[(phase,group,roi),'rsa'],
            color=face,zorder=10,alpha=.8,marker=hpc_groups[group]['marker'],
            edgecolors=edge,s=rcParams['lines.markersize']**1.7)

    #subject chicken scratch
    for phase, p_ in zip(['acquisition','extinction'],[0,1]):
        for sub in hpc_groups[group]['subs']:
            yvals = subdf.loc[(phase,group,slice('dACC','vmPFC'),sub),'rsa'].values
            ax[r].plot([p_-.1,p_+.1],yvals,color='grey',alpha=.2)

    #bars
    bars = ax[r].bar(xvals,df.loc[group,'estimate'].values,width=.4,color=pal26)
    for b, bar in enumerate(bars): 
        if b & 1: bar.set_hatch('////')
        bar.set_edgecolor('k')

    #adding the CI from the LMM estimates
    ax[r].vlines(xvals,df.loc[group,'asymp.LCL'].values,df.loc[group,'asymp.UCL'].values,color='k',zorder=4,capstyle='round')

    #labels and such
    ax[r].yaxis.set_major_locator(MultipleLocator(.15))   
    ax[r].set_xticks(range(2))
    ax[r].set_ylabel('Memory reinstatement',labelpad=1)
    ax[r].set_xticklabels(['Conditioning','Extinction'],ha='center')
    ax[r].set_title(hpc_groups[group]['name'],pad=4)


    #putting in the CS+/- descriptive arrows
    ax[r].set_xlim(-.7,ax[r].get_xlim()[1])
    ax[r].arrow(-.55,0,0,.05,head_width=0.06, head_length=0.025, fc='grey', ec='grey')
    ax[r].arrow(-.55,0,0,-.05,head_width=0.06, head_length=0.025, fc='grey', ec='grey')
    ax[r].text(-.55,.075,'more\nCS+',ha='center',va='bottom',fontsize=8)
    ax[r].text(-.55,-.075,'more\nCS-',ha='center',va='top',fontsize=8)

    top = df.loc[group,'asymp.UCL'].values.copy()
    bottom = df.loc[group,'asymp.LCL'].values.copy()
    if group == 'ptsd': top[1] = bottom[1] - .04
    for star, x, y in zip(hpc_groups[group]['stars'],xvals,top): ax[r].text(x,y,star,ha='center',va='bottom',fontsize=10)

    sns.despine(ax=ax[r],trim=True,bottom=True)
    ax[r].yaxis.set_major_formatter(FuncFormatter(no_leading_zeros))

legend_elements = [Patch(facecolor='k',edgecolor='k',label='dACC'),
                   Patch(facecolor='w',edgecolor='k',hatch='////',label='vmPFC')]

fig.legend(handles=legend_elements,ncol=1,frameon=False,loc='center left',bbox_to_anchor=(0.025,.5))

#on diag is [**,*,*]
paired_barplot_annotate_brackets('**', 0, df.loc[('healthy','conditioning'),'asymp.UCL'].values, ax[0].get_ylim(), xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[0])
paired_barplot_annotate_brackets('*', 1, df.loc[('healthy','extinction'),'asymp.UCL'].values, ax[0].get_ylim(), xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[0])
paired_barplot_annotate_brackets('*', 0, df.loc[('ptsd','conditioning'),'asymp.UCL'].values, ax[1].get_ylim(), xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[1])

ax[0].scatter(.5,.7,marker='o',s=35,c='white',edgecolors='k')
ax[0].scatter(.5,.7,marker='x',s=30,c='k')

plt.tight_layout()

'''behavioral plots'''
df = pd.read_csv('day1_behavior.csv')
df = df[df.phase.isin(['fear','early_ext','late_ext'])]
df.group = df.group.apply(lambda x: 'ptsd' if x == 'ptss' else x)
df.phase = pd.Categorical(df.phase,['fear','early_ext','late_ext'],ordered=True)
est = df.groupby(['group','phase']).mean()
err = df.groupby(['group','phase']).sem() * 1.96

memdf = pd.read_csv('mem_hit_rate_sub_means.csv')
memdf = memdf[memdf.phase != 'baseline']
memdf = memdf.set_index(['condition','group','phase','subject'])
memdf = memdf.loc['CS+'] - memdf.loc['CS-']
memest = memdf.groupby(['group','phase']).mean()
memerr = memdf.groupby(['group','phase']).sem() * 1.96

xvals = [[.95,1.95,2.95],[1.05,2.05,3.05]]
memxvals = [[.95,1.95],[1.05,2.05]]

fig, ax = plt.subplots(1,3,figsize=(mm2inch(173,60)),sharey=False)
for g, group in enumerate(hpc_groups):
    for v, val in enumerate(['scr','exp']):
        ax[v].scatter(xvals[g],est.loc[group,val],marker=hpc_groups[group]['marker'],color=hpc_groups[group]['c'],zorder=10)
        ax[v].plot(xvals[g],est.loc[group,val],color=hpc_groups[group]['c'],zorder=1,linestyle=hpc_groups[group]['line'])
        ax[v].vlines(xvals[g],est.loc[group,val] - err.loc[group,val],est.loc[group,val] + err.loc[group,val],
            color=hpc_groups[group]['c'],linestyle=hpc_groups[group]['line'],zorder=9)

        ax[v].set_xticks(range(1,4))
        ax[v].set_xticklabels(['Conditioning','Early','Late'])

    ax[2].scatter(memxvals[g],memest.loc[group,'mem_acc'],marker=hpc_groups[group]['marker'],color=hpc_groups[group]['c'],zorder=10)
    ax[2].plot(memxvals[g],memest.loc[group,'mem_acc'],color=hpc_groups[group]['c'],zorder=1,linestyle=hpc_groups[group]['line'],marker=hpc_groups[group]['marker'],label=hpc_groups[group]['name'])
    ax[2].vlines(memxvals[g],memest.loc[group,'mem_acc'] - memerr.loc[group,'mem_acc'],memest.loc[group,'mem_acc'] + memerr.loc[group,'mem_acc'],
            color=hpc_groups[group]['c'],linestyle=hpc_groups[group]['line'],zorder=9)
    ax[2].set_xticks(range(1,3))
    ax[2].set_xticklabels(['Conditioning','Extinction'])

for a, Ax in enumerate(ax):
    if a == 0:
        Ax.set_xlim(0,Ax.get_xlim()[1])
        y2, y1 = Ax.get_ylim()
        Ax.arrow(.5,0,0,y1*.03,head_width=.06, head_length=y1*0.02, fc='grey', ec='grey')
        Ax.arrow(.5,0,0,y1*-.03,head_width=.06, head_length=y1*0.02, fc='grey', ec='grey')
        Ax.text(.5,y1*.055,'more\nCS+',ha='center',va='bottom',fontsize=8)
        Ax.text(.5,y1*-.05,'more\nCS-',ha='center',va='top',fontsize=8)

    Ax.set_ylim(Ax.get_ylim()[0]*1.15,Ax.get_ylim()[1])
    sns.despine(ax=Ax,trim=True,bottom=True)
    Ax.hlines(0,*Ax.get_xlim(),color='k')  if a != 0 else Ax.hlines(0,.5,Ax.get_xlim()[1],color='k')
    Ax.yaxis.set_major_locator(MultipleLocator(.10))   
    Ax.yaxis.set_major_formatter(FuncFormatter(no_leading_zeros))
    
#labels and such
ax[0].set_title('Autonomic Arousal\n(Day 1)',pad=.1)
ax[0].set_ylabel('Sqrt( SCR )')
ax[1].set_title('Shock expectancy\n(Day 1)',pad=.1)
ax[1].set_ylabel('Mean expectancy')
ax[2].set_title('Recognition memory\n(Day 2)',pad=.1)
ax[2].set_ylabel('Hit rate')
#ax[2].legend(frameon=False)

#align the axes
align_yaxis(ax[0],0,ax[1],0,0,ax[1].get_ylim()[1])
align_yaxis(ax[0],0,ax[2],0,0,ax[2].get_ylim()[1])
plt.tight_layout()


'''hippocampal plot just like the pfc bars'''
df = pd.read_csv('hpc_phasedif_stats.csv').rename(columns={'emmean':'estimate'})
df.roi = df.roi.apply(lambda x: x[:-4])
df = df[df.phase != 'baseline'][df.roi != 'hc_body']
df.roi = pd.Categorical(df.roi,['hc_tail','hc_head'],ordered=True)
df = df.set_index(['group','phase','roi']).sort_index()

hpc_groups = {'healthy':{'name':'Healthy','subs':sub_args,'marker':'D','line':'-','c':'k'},
                 'ptsd':{'name':'PTSS','subs':p_sub_args,'marker':'s','line':'--','c':'grey'}}

xvals = [(i-.2,i+.2) for i in range(2)]
xvals = [i for j in xvals for i in j]
pal2 = sns.color_palette(['red','dodgerblue'],desat=.8)
pal26 = [pal2[0],pal2[0],pal2[1],pal2[1]]

subdf = pd.read_csv('subcort_ers_cleaned_lmm.csv').groupby(['phase','group','roi','subject']).mean()
subdf = subdf.drop(index='baseline',columns='uni').reset_index()
subdf = subdf[subdf.roi.isin(['hc_head','hc_tail'])]
subdf.roi = pd.Categorical(subdf.roi,['hc_tail','hc_head'],ordered=True)
subdf = subdf.set_index(['phase','group','roi','subject']).sort_index()

subxvals = [(i+.1,i-.1) for i in range(2)]
subxvals = [i for j in subxvals for i in j]

fig, ax = plt.subplots(2,1,figsize=(mm2inch(86,150)),sharey=True)
for r, group in enumerate(hpc_groups):

    #subject data points
    for bar, (phase, roi) in enumerate(product(['acquisition','extinction'],['hc_head','hc_tail'])):
        # edge = pal26[bar] if not bar & 1 else 'k'
        # face = 'white' if not bar & 1 else pal26[bar]
        edge=pal26[bar];face='white';
        ax[r].scatter(np.repeat(subxvals[bar],24),subdf.loc[(phase,group,roi),'ers'],
            color=face,zorder=10,alpha=.8,marker=hpc_groups[group]['marker'],
            edgecolors=edge,s=rcParams['lines.markersize']**1.7)

    #subject chicken scratch
    for phase, p_ in zip(['acquisition','extinction'],[0,1]):
        for sub in hpc_groups[group]['subs']:
            yvals = subdf.loc[(phase,group,slice('hc_tail','hc_head'),sub),'ers'].values
            ax[r].plot([p_-.1,p_+.1],yvals,color='grey',alpha=.2)

    #bars
    bars = ax[r].bar(xvals,df.loc[group,'estimate'].values,width=.4,color=pal26)
    for b, bar in enumerate(bars): 
        if b & 1: bar.set_hatch('////')
        bar.set_edgecolor('k')

    #adding the CI from the LMM estimates
    ax[r].vlines(xvals,df.loc[group,'asymp.LCL'].values,df.loc[group,'asymp.UCL'].values,color='k',zorder=4,capstyle='round')

    #labels and such
    ax[r].yaxis.set_major_locator(MultipleLocator(.15))   
    ax[r].set_xticks(range(2))
    ax[r].set_ylabel('Memory reinstatement',labelpad=1)
    ax[r].set_xticklabels(['Conditioning','Extinction'],ha='center')
    ax[r].set_title(hpc_groups[group]['name'],pad=4)

    sns.despine(ax=ax[r],trim=True,bottom=True)
    ax[r].yaxis.set_major_formatter(FuncFormatter(no_leading_zeros))

legend_elements = [Patch(facecolor='k',edgecolor='k',label='pHC'),
                   Patch(facecolor='w',edgecolor='k',hatch='////',label='aHC')]

fig.legend(handles=legend_elements,ncol=1,frameon=False,loc='center left',bbox_to_anchor=(0.015,.525))
fig.text(.05,.475,'Collapsed across CS',ha='left',va='bottom',fontsize=8)

paired_barplot_annotate_brackets('*', 0, df.loc[('healthy','conditioning'),'asymp.UCL'].values, ax[0].get_ylim(),xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[0])
paired_barplot_annotate_brackets('***', 1, df.loc[('healthy','extinction'),'asymp.UCL'].values, ax[0].get_ylim(),xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[0])
paired_barplot_annotate_brackets('*', 0, df.loc[('ptsd','conditioning'),'asymp.UCL'].values, ax[1].get_ylim(),xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[1])

ax[0].scatter(.5,.3,marker='o',s=35,c='white',edgecolors='k')
ax[0].scatter(.5,.3,marker='x',s=30,c='k')

plt.tight_layout()
