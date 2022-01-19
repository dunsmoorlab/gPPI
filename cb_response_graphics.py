from fg_config import *
from paper_graphics_style import *
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from matplotlib import cm, gridspec
from itertools import product
import matplotlib as mpl
import pingouin as pg

df = pd.read_csv('pfc_csdif_off_stats.csv')
df = df[df.phase.isin(['conditioning','extinction'])]
df = df.set_index(['group','phase','roi']).sort_index()

pfc_groups = {'healthy':{'name':'Healthy','stars':['***','***','**','***'],'subs':sub_args,'marker':'D','line':'-','c':'k'},
                 'ptsd':{'name':'PTSS','stars':['***','','***','**'],'subs':p_sub_args,'marker':'s','line':'--','c':'grey'}}

xvals = [(i-.2,i+.2) for i in range(2)]
xvals = [i for j in xvals for i in j]
pal2 = sns.color_palette(['red','dodgerblue'],desat=.8)
pal26 = [pal2[0],pal2[0],pal2[1],pal2[1]]

subdf = pd.read_csv('cb_response_rsa.csv').groupby(['condition','phase','group','roi','subject']).mean()
subdf = (subdf.loc['CS+'] - subdf.loc['CS-']).reset_index()
subdf = subdf[subdf.roi.isin(['dACC','vmPFC'])].set_index(['phase','group','roi','subject']).sort_index()

subxvals = [(i+.1,i-.1) for i in range(2)]
subxvals = [i for j in subxvals for i in j]

fig, ax = plt.subplots(2,1,figsize=(mm2inch(86,150)),sharey=True)
for r, group in enumerate(pfc_groups):

    #subject data points
    for bar, (phase, roi) in enumerate(product(['conditioning','extinction'],['vmPFC','dACC'])):
        # edge = pal26[bar] if not bar & 1 else 'k'
        # face = 'white' if not bar & 1 else pal26[bar]
        edge=pal26[bar];face='white';
        ax[r].scatter(np.repeat(subxvals[bar],24),subdf.loc[(phase,group,roi),'off_ers'],
            color=face,zorder=10,alpha=.8,marker=pfc_groups[group]['marker'],
            edgecolors=edge,s=rcParams['lines.markersize']**1.7)

    #subject chicken scratch
    for phase, p_ in zip(['conditioning','extinction'],[0,1]):
        for sub in pfc_groups[group]['subs']:
            yvals = subdf.loc[(phase,group,slice('dACC','vmPFC'),sub),'off_ers'].values
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
    ax[r].set_title(pfc_groups[group]['name'],pad=4)


    #putting in the CS+/- descriptive arrows
    ax[r].set_xlim(-.7,ax[r].get_xlim()[1])
    ax[r].arrow(-.55,0,0,.05,head_width=0.06, head_length=0.025, fc='grey', ec='grey')
    ax[r].arrow(-.55,0,0,-.05,head_width=0.06, head_length=0.025, fc='grey', ec='grey')
    ax[r].text(-.55,.075,'more\nCS+',ha='center',va='bottom',fontsize=8)
    ax[r].text(-.55,-.075,'more\nCS-',ha='center',va='top',fontsize=8)

    top = df.loc[group,'asymp.UCL'].values.copy()
    bottom = df.loc[group,'asymp.LCL'].values.copy()
    # if group == 'ptsd': top[1] = bottom[1] - .04
    for star, x, y in zip(pfc_groups[group]['stars'],xvals,top): ax[r].text(x,y,star,ha='center',va='bottom',fontsize=10)

    sns.despine(ax=ax[r],trim=True,bottom=True)
    ax[r].yaxis.set_major_formatter(FuncFormatter(no_leading_zeros))

legend_elements = [Patch(facecolor='k',edgecolor='k',label='dACC'),
                   Patch(facecolor='w',edgecolor='k',hatch='////',label='vmPFC')]

fig.legend(handles=legend_elements,ncol=1,frameon=False,loc='center left',bbox_to_anchor=(0.025,.5))

paired_barplot_annotate_brackets('***', 0, df.loc[('healthy','conditioning'),'asymp.UCL'].values, ax[0].get_ylim(), xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[0])
paired_barplot_annotate_brackets('**', 1, df.loc[('healthy','extinction'),'asymp.UCL'].values, ax[0].get_ylim(), xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[0])

paired_barplot_annotate_brackets('***', 0, df.loc[('ptsd','conditioning'),'asymp.UCL'].values, ax[1].get_ylim(), xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[1])

ax[0].scatter(.5,.6,marker='o',s=35,c='white',edgecolors='k')
ax[0].scatter(.5,.6,marker='x',s=30,c='k')

plt.tight_layout()



'''insula and precuneus graphs'''
df = pd.read_csv('ant_ins_precun_stats.csv')
df = df.set_index(['group','roi','phase']).sort_index()

subdf = pd.read_csv('group_ant_ins_precun_ers.csv').groupby(['condition','roi','group','phase','subject']).mean()
subdf = (subdf.loc['CS+'] - subdf.loc['CS-']).sort_index()

xvals = [(i-.2,i+.2) for i in range(2)] 
xvals = [i for j in xvals for i in j]
pal2 = sns.color_palette(['red','dodgerblue'],desat=.8)
pal26 = [pal2[0],pal2[1],pal2[0],pal2[1]]

subxvals = [(i-.1,i+.1) for i in range(2)]
subxvals = [i for j in subxvals for i in j]

cb_groups = {'healthy':{'name':'Healthy','stars':['***','***','~','***'],'subs':sub_args,'marker':'D','line':'-','c':'k'},
                 'ptsd':{'name':'PTSS','stars':['***','**','','***'],'subs':p_sub_args,'marker':'s','line':'--','c':'grey'}}


fig, ax = plt.subplots(2,1,figsize=(mm2inch(86,150)),sharey=True)
for r, group in enumerate(cb_groups):

    #subject data points
    for bar, (roi, phase) in enumerate(product(['ant_ins','precun'],['conditioning','extinction'])):
        edge = pal26[bar] if not bar & 1 else 'k'
        face = 'white' if not bar & 1 else pal26[bar]
        edge=pal26[bar];face='white';
        ax[r].scatter(np.repeat(subxvals[bar],24),subdf.loc[(roi,group,phase),'ers'],
            color=face,zorder=10,alpha=.8,marker=pfc_groups[group]['marker'],
            edgecolors=edge,s=rcParams['lines.markersize']**1.7)
    #subject chicken scratch
    for roi, p_ in zip(['ant_ins','precun'],[0,1]):
        for sub in pfc_groups[group]['subs']:
            yvals = subdf.loc[(roi,group,slice('conditioning','extinction'),sub),'ers'].values
            ax[r].plot([p_-.1,p_+.1],yvals,color='grey',alpha=.2)

    #bars
    bars = ax[r].bar(xvals,df.loc[group,'estimate'].values,width=.4,color=pal26)
    for b, bar in enumerate(bars): 
        # if b & 1: bar.set_hatch('////')
        bar.set_edgecolor('k')
    
    ax[r].vlines(xvals,df.loc[group,'asymp.LCL'].values,df.loc[group,'asymp.UCL'].values,color='k',zorder=4,capstyle='round')

    #labels and such
    ax[r].yaxis.set_major_locator(MultipleLocator(.15))   
    ax[r].set_xticks(range(2))
    ax[r].set_ylabel('Memory reinstatement',labelpad=1)
    ax[r].set_xticklabels(['Anterior insula','Precuneus'],ha='center')
    ax[r].set_title(pfc_groups[group]['name'],pad=4)


    #putting in the CS+/- descriptive arrows
    ax[r].set_xlim(-.7,ax[r].get_xlim()[1])
    ax[r].arrow(-.55,0,0,.05,head_width=0.06, head_length=0.025, fc='grey', ec='grey')
    ax[r].arrow(-.55,0,0,-.05,head_width=0.06, head_length=0.025, fc='grey', ec='grey')
    ax[r].text(-.55,.075,'more\nCS+',ha='center',va='bottom',fontsize=8)
    ax[r].text(-.55,-.075,'more\nCS-',ha='center',va='top',fontsize=8)

    top = df.loc[group,'asymp.UCL'].values.copy()
    bottom = df.loc[group,'asymp.LCL'].values.copy()
    # if group == 'ptsd': top[1] = bottom[1] - .04
    for star, x, y in zip(cb_groups[group]['stars'],xvals,top): ax[r].text(x,y,star,ha='center',va='bottom',fontsize=10)

    sns.despine(ax=ax[r],trim=True,bottom=True)
    ax[r].yaxis.set_major_formatter(FuncFormatter(no_leading_zeros))

legend_elements = [Patch(facecolor=pal26[0],edgecolor='k',label='Conditioning'),
                   Patch(facecolor=pal26[1],edgecolor='k',label='Extinction')]

fig.legend(handles=legend_elements,ncol=1,frameon=False,loc='center left',bbox_to_anchor=(0.025,.475))

paired_barplot_annotate_brackets('**', 0, df.loc[('healthy','ant_ins'),'asymp.UCL'].values, ax[1].get_ylim(), xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[0])
paired_barplot_annotate_brackets('**', 0, df.loc[('ptsd','ant_ins'),'asymp.UCL'].values, ax[1].get_ylim(), xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[1])
paired_barplot_annotate_brackets('~', 1, df.loc[('ptsd','precun'),'asymp.UCL'].values, ax[1].get_ylim(), xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[1])

plt.tight_layout()


'''cross phase similarity'''
df = pd.read_csv('pfc_cross_stats.csv')
df = df[df.phase.isin(['conditioning','extinction'])]
df = df.set_index(['group','phase','roi']).sort_index()

pfc_groups = {'healthy':{'name':'Healthy','stars':['***','***','***','***'],'subs':sub_args,'marker':'D','line':'-','c':'k'},
                 'ptsd':{'name':'PTSS','stars':['***','***','***','***'],'subs':p_sub_args,'marker':'s','line':'--','c':'grey'}}

xvals = [(i-.2,i+.2) for i in range(2)]
xvals = [i for j in xvals for i in j]
pal2 = sns.color_palette(['red','dodgerblue'],desat=.8)
pal26 = [pal2[0],pal2[0],pal2[1],pal2[1]]

subdf = pd.read_csv('cb_cross_rsa.csv').groupby(['condition','phase','group','roi','subject']).mean()
subdf = (subdf.loc['CS+'] - subdf.loc['CS-']).reset_index()
subdf = subdf[subdf.roi.isin(['dACC','vmPFC'])].set_index(['phase','group','roi','subject']).sort_index()

subxvals = [(i+.1,i-.1) for i in range(2)]
subxvals = [i for j in subxvals for i in j]

fig, ax = plt.subplots(2,1,figsize=(mm2inch(86,150)),sharey=True)
for r, group in enumerate(pfc_groups):

    #subject data points
    for bar, (phase, roi) in enumerate(product(['conditioning','extinction'],['vmPFC','dACC'])):
        # edge = pal26[bar] if not bar & 1 else 'k'
        # face = 'white' if not bar & 1 else pal26[bar]
        edge=pal26[bar];face='white';
        ax[r].scatter(np.repeat(subxvals[bar],24),subdf.loc[(phase,group,roi),'cross_ers'],
            color=face,zorder=10,alpha=.8,marker=pfc_groups[group]['marker'],
            edgecolors=edge,s=rcParams['lines.markersize']**1.7)

    #subject chicken scratch
    for phase, p_ in zip(['conditioning','extinction'],[0,1]):
        for sub in pfc_groups[group]['subs']:
            yvals = subdf.loc[(phase,group,slice('dACC','vmPFC'),sub),'cross_ers'].values
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
    ax[r].set_ylabel('Similarity',labelpad=1)
    ax[r].set_xticklabels(['Conditioning','Extinction'],ha='center')
    ax[r].set_title(pfc_groups[group]['name'],pad=4)


    #putting in the CS+/- descriptive arrows
    ax[r].set_xlim(-.7,ax[r].get_xlim()[1])
    ax[r].arrow(-.55,0,0,.05,head_width=0.06, head_length=0.025, fc='grey', ec='grey')
    ax[r].arrow(-.55,0,0,-.05,head_width=0.06, head_length=0.025, fc='grey', ec='grey')
    ax[r].text(-.55,.075,'more\nCS+',ha='center',va='bottom',fontsize=8)
    ax[r].text(-.55,-.075,'more\nCS-',ha='center',va='top',fontsize=8)

    top = df.loc[group,'asymp.UCL'].values.copy()
    bottom = df.loc[group,'asymp.LCL'].values.copy()
    # if group == 'ptsd': top[1] = bottom[1] - .04
    for star, x, y in zip(pfc_groups[group]['stars'],xvals,top): ax[r].text(x,y,star,ha='center',va='bottom',fontsize=10)

    sns.despine(ax=ax[r],trim=True,bottom=True)
    ax[r].yaxis.set_major_formatter(FuncFormatter(no_leading_zeros))

legend_elements = [Patch(facecolor='k',edgecolor='k',label='dACC'),
                   Patch(facecolor='w',edgecolor='k',hatch='////',label='vmPFC')]

fig.legend(handles=legend_elements,ncol=1,frameon=False,loc='center left',bbox_to_anchor=(0.025,.5))

paired_barplot_annotate_brackets('***', 0, df.loc[('ptsd','conditioning'),'asymp.UCL'].values, ax[1].get_ylim(), xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[1])
paired_barplot_annotate_brackets('**', 1, df.loc[('ptsd','extinction'),'asymp.UCL'].values, ax[1].get_ylim(), xtick_spread=.2, dh=.08, barh=.025, fs=10, ax=ax[1])

ax[0].vlines(.5,-.15,.75,colors='k',linestyles='--')
ax[1].vlines(.5,-.15,.75,colors='k',linestyles='--')

ax[0].text(0,.55,'conditioning encoding\nsimilarity to\nextinction retrieval',ha='center',va='center',fontsize=8)
ax[0].text(1,.55,'extinction encoding\nsimilarity to\nconditioning retrieval',ha='center',va='center',fontsize=8)


ax[1].scatter(.5,.6,marker='o',s=35,c='white',edgecolors='k')
ax[1].scatter(.5,.6,marker='x',s=30,c='k')

plt.tight_layout()


'''on vs off diagonal comparison'''
df = pd.read_csv('roi_on_off_diff_stats.csv').sort_values(by='emmean',ascending=False)
xvals = list(range(10))
pal = sns.color_palette('crest',n_colors=10)

fig, ax = plt.subplots(figsize=(mm2inch(173,60)))
ax.bar(xvals,df.emmean.values,color=pal)
ax.vlines(xvals,df['asymp.LCL'].values,df['asymp.UCL'].values,color='k',zorder=4,capstyle='round')
sns.despine(ax=ax,trim=True,bottom=True)
ax.yaxis.set_major_formatter(FuncFormatter(no_leading_zeros))
ax.set_xlim(-1,ax.get_xlim()[1])
ax.arrow(-.7,0,0,.005,head_width=0.15, head_length=0.0025, fc='grey', ec='grey')
ax.arrow(-.7,0,0,-.005,head_width=0.15, head_length=0.0025, fc='grey', ec='grey')
ax.text(-.7,.008,'more\nItem',ha='center',va='bottom',fontsize=8)
ax.text(-.7,-.008,'more\nSet',ha='center',va='top',fontsize=8)
ax.set_ylabel('Item vs. set level reinstatement',labelpad=1)
ax.set_xticks(xvals)
ax.set_xticklabels(['Fusiform','BLA','CeM','HC body','vmPFC','dACC','pHC','aHC','Precun.','Ant. ins.'],ha='center')
ax.text(0,df['asymp.UCL'].values[0],'**',ha='center',va='bottom',fontsize=10)
plt.tight_layout()

'''revisting a plot of the subcortical ers predictions'''
fig, ax1 = plt.subplots(figsize=mm2inch(55,60))
ers = pd.read_csv('ers_slopes_pfc_pred.csv').drop(columns='Unnamed: 0').rename(index={0:'est',1:'LCL',2:'UCL'}).swapaxes('index','columns').reset_index().rename(columns={'index':'roi'})
ers['stat'] = 'ers'

yvals = [i for i in range(5)]
ax1.scatter(ers.est.values,yvals,s=8,c='k')
ax1.hlines(yvals,ers.LCL.values,ers.UCL.values,color='k',capstyle='round')
ax1.set_xlabel('Predictiveness on mPFC\ndifference in reinstatement')

ax1.set_yticks(range(5))
sns.despine(ax=ax1,left=True)
ax1.set_yticklabels(['pHC','HC Body','aHC','BLA','CeM'],ha='right')

ers_stars = [('*',-1),('',1),('***',1),('',1),('**',-1)]


for s, pair in enumerate(ers_stars):
    if pair[1] == -1:
        ax1.text(ers.LCL.values[s]-.005,yvals[s]-.18,pair[0],ha='right',va='bottom',fontsize=10)        
    else:
        ax1.text(ers.UCL.values[s]+.005,yvals[s]-.18,pair[0],ha='left',va='bottom',fontsize=10)        

ax1.set_ylim(ax1.get_ylim()[0],4.7)
ax1.vlines(0,*ax1.get_ylim(),linestyle=':',color='black',zorder=0)
ax1.arrow(0,4.45,0.012,0,head_width=0.15, head_length=0.0045, fc='grey', ec='grey')
ax1.arrow(0,4.45,-0.012,0,head_width=0.15, head_length=0.0045, fc='grey', ec='grey')
ax1.text(.018,4.45,'predicts more\nvmPFC',ha='left',va='center',fontsize=8)
ax1.text(-.018,4.45,'predicts more\ndACC',ha='right',va='center',fontsize=8)

ax1.set_ylabel('MTL reinstatement')

plt.tight_layout()

'''additional behavioral plot   for the supplement'''
df = pd.read_csv('day1_behavior.csv')
df.group = df.group.apply(lambda x: 'ptsd' if x == 'ptss' else x)
df.phase = pd.Categorical(df.phase,['fear','early_ext','late_ext','early_rnw'],ordered=True)
est = df.groupby(['group','phase']).mean()
err = df.groupby(['group','phase']).sem() * 1.96


xvals = [[.95,1.95,2.95,3.95],[1.05,2.05,3.05,4.05]]

hpc_groups = {'healthy':{'name':'Healthy','subs':sub_args,'marker':'D','line':'-','c':'k'},
                 'ptsd':{'name':'PTSS','subs':p_sub_args,'marker':'s','line':'--','c':'grey'}}

fig, ax = plt.subplots(2,2,figsize=(mm2inch(173,200)),sharey=False)
for g, group in enumerate(hpc_groups):
    for v, val in enumerate(['scr','exp']):
        ax[0][v].scatter(xvals[g],est.loc[group,val],marker=hpc_groups[group]['marker'],color=hpc_groups[group]['c'],zorder=10)
        ax[0][v].plot(xvals[g],est.loc[group,val],color=hpc_groups[group]['c'],zorder=1,linestyle=hpc_groups[group]['line'])
        ax[0][v].vlines(xvals[g],est.loc[group,val] - err.loc[group,val],est.loc[group,val] + err.loc[group,val],
            color=hpc_groups[group]['c'],linestyle=hpc_groups[group]['line'],zorder=9)

        ax[0][v].set_xticks(range(1,5))
        ax[0][v].set_xticklabels(['Conditioning','Early','Late','Early\nRenewal'])


for a, Ax in enumerate(ax[0]):
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
align_yaxis(ax[0][0],0,ax[0][1],0,0,ax[0][1].get_ylim()[1])
ax[0][0].set_title('Autonomic Arousal\n(Day 1)',pad=.1)
ax[0][0].set_ylabel('Sqrt( SCR )')
ax[0][1].set_title('Shock expectancy\n(Day 1)',pad=.1)
ax[0][1].set_ylabel('Mean expectancy')


memdf = pd.read_csv('mem_hit_rate_sub_means.csv')
memdf.phase = pd.Categorical(memdf.phase,['baseline','acquisition','extinction'],ordered=True)
memdf = memdf.set_index(['condition','group','phase','subject'])
memdf = memdf.loc['CS+'] - memdf.loc['CS-']
memest = memdf.groupby(['group','phase']).mean()
memerr = memdf.groupby(['group','phase']).sem() * 1.96
memxvals = [[.95,1.95,2.95],[1.05,2.05,3.05]]

for g, group in enumerate(hpc_groups):
    ax[1][0].scatter(memxvals[g],memest.loc[group,'mem_acc'],marker=hpc_groups[group]['marker'],color=hpc_groups[group]['c'],zorder=10)
    ax[1][0].plot(memxvals[g],memest.loc[group,'mem_acc'],color=hpc_groups[group]['c'],zorder=1,linestyle=hpc_groups[group]['line'],marker=hpc_groups[group]['marker'],label=hpc_groups[group]['name'])
    ax[1][0].vlines(memxvals[g],memest.loc[group,'mem_acc'] - memerr.loc[group,'mem_acc'],memest.loc[group,'mem_acc'] + memerr.loc[group,'mem_acc'],
            color=hpc_groups[group]['c'],linestyle=hpc_groups[group]['line'],zorder=9)
ax[1][0].set_xticks(range(1,4))
ax[1][0].set_xticklabels(['Pre-\nconditioning','Conditioning','Extinction'])


ax[1][0].set_title('Recognition memory\n(Day 2)',pad=.1)
ax[1][0].set_ylabel('Hit rate')
ax[1][0].set_xlim(0,ax[1][0].get_xlim()[1])
y2, y1 = ax[1][0].get_ylim()
ax[1][0].arrow(.5,0,0,y1*.03,head_width=.06, head_length=y1*0.02, fc='grey', ec='grey')
ax[1][0].arrow(.5,0,0,y1*-.03,head_width=.06, head_length=y1*0.02, fc='grey', ec='grey')
ax[1][0].text(.5,y1*.055,'more\nCS+',ha='center',va='bottom',fontsize=8)
ax[1][0].text(.5,y1*-.05,'more\nCS-',ha='center',va='top',fontsize=8)

ax[1][0].set_ylim(ax[1][0].get_ylim()[0]*1.15,ax[1][0].get_ylim()[1])
sns.despine(ax=ax[1][0],trim=True,bottom=True)
ax[1][0].hlines(0,.5,ax[1][0].get_xlim()[1],color='k')
ax[1][0].yaxis.set_major_locator(MultipleLocator(.10))  
ax[1][0].yaxis.set_major_formatter(FuncFormatter(no_leading_zeros))


fadf = pd.read_csv('false_alarm_sub_means.csv')
fadf.condition = pd.Categorical(fadf.condition,['CS+','CS-'],ordered=True)
faest = fadf.groupby(['group','condition']).mean()
faerr = fadf.groupby(['group','condition']).sem() * 1.96
faxvals = [[.95,1.95],[1.05,2.05]]

for g, group in enumerate(hpc_groups):
    ax[1][1].scatter(faxvals[g],faest.loc[group,'fa'],marker=hpc_groups[group]['marker'],color=hpc_groups[group]['c'],zorder=10)
    ax[1][1].plot(faxvals[g],faest.loc[group,'fa'],color=hpc_groups[group]['c'],zorder=1,linestyle=hpc_groups[group]['line'],marker=hpc_groups[group]['marker'],label=hpc_groups[group]['name'])
    ax[1][1].vlines(faxvals[g],faest.loc[group,'fa'] - faerr.loc[group,'fa'],faest.loc[group,'fa'] + faerr.loc[group,'fa'],
            color=hpc_groups[group]['c'],linestyle=hpc_groups[group]['line'],zorder=9)
ax[1][1].set_xticks(range(1,3))
ax[1][1].set_xticklabels(['CS+','CS-'])
ax[1][1].set_xlim([.5,2.5])

ax[1][1].set_title('False alarms\n(Day 2)',pad=.1)
ax[1][1].set_ylabel('False alarm rate')

ax[1][1].set_ylim(0,ax[1][1].get_ylim()[1])
sns.despine(ax=ax[1][1],trim=True,bottom=True)

ax[1][1].yaxis.set_major_locator(MultipleLocator(.10))  
ax[1][1].yaxis.set_major_formatter(FuncFormatter(no_leading_zeros))
align_yaxis(ax[1][0],0,ax[1][1],0,0,ax[1][1].get_ylim()[1])

plt.tight_layout()
