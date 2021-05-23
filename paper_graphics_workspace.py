from fg_config import *
from paper_graphics_style import *
import pingouin as pg
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from matplotlib import cm, gridspec
import matplotlib as mpl
# #looking at pfc ers all phases with bars
# paired_pal = [phase_pal[0], [i*1.75 for i in phase_pal[0]], phase_pal[1],[i*1.75 for i in phase_pal[1]], phase_pal[2], [i*1.75 for i in phase_pal[2]]]
# '''hippocampus is key'''
# df = pd.read_csv('all_data_lmm.csv').set_index(['group','phase','condition','subject','stimulus']).sort_index()
# df = df.loc[(['healthy','ptsd'],['acquisition','baseline','extinction']),].reset_index().drop(columns='stimulus')
# df = df.groupby(['group','phase','subject']).mean().reset_index()
# # df['phase_con'] = df.phase + '_' + df.condition

# _rois = {'hc_tail':'Tail (Posterior)', 'hc_head':'Head (Anterior)'}
# fig, ax = plt.subplots(1,len(_rois),figsize=mm2inch(200,100),sharey=True)
# for r, roi in enumerate(_rois):    
#     sns.pointplot(data=df,x='phase',y=f'{roi}_ers',hue='group',palette="Set2",ax=ax[r],order=phase3,dodge=True)#['baseline_CS+','baseline_CS-','acquisition_CS+','acquisition_CS-','extinction_CS+','extinction_CS-'])
#     # sns.stripplot(data=df,x='group',y=f'{roi}_ers',hue='phase_con',dodge=True,color='black',size=2,ax=ax[r])
#     # ax[r].legend(loc='center left',bbox_to_anchor=(1,.5)) if r == len(_rois) - 1 else ax[r].legend_.remove()
#     ax[r].legend_.remove()
#     ax[r].set_ylabel('Encoding-retrieval similarity') if r == 0 else ax[r].set_ylabel('')
#     ax[r].set_title(f'{_rois[roi]}')
#     sns.despine(ax=ax[r],left=True)
#     ax[r].yaxis.set_major_locator(MultipleLocator(.05))   

df = pd.read_csv('pfc_csdif_stats.csv')
df.roi = df.roi.apply(lambda x: x[:-4])
df = df.set_index(['group','phase','roi']).sort_index()

pfc_groups = {'healthy':{'name':'Healthy','stars':['','','***','*','','***']},
                 'ptsd':{'name':'PTSS','stars':['','*','***','','**','']}}

xvals = [(i-.2,i+.2) for i in range(3)]
xvals = [i for j in xvals for i in j]
pal2 = sns.color_palette(['red','dodgerblue'],desat=.8)
pal26 = pal2+pal2+pal2


fig, ax = plt.subplots(2,1,figsize=(mm2inch(86,150)),sharey=True)
for r, group in enumerate(pfc_groups):

    # ax[r].scatter(xvals,df.loc[group,'estimate'].values,c=pal26,s=8)
    # ax[r].plot(xvals[::2],df.loc[group,'estimate'].values[::2],color=pal2[0],linewidth=rcParams['lines.linewidth']*.75,zorder=1)
    # ax[r].plot(xvals[1::2],df.loc[group,'estimate'].values[1::2],color=pal2[1],linewidth=rcParams['lines.linewidth']*.75,zorder=2)

    ax[r].bar(xvals,df.loc[group,'estimate'].values,width=.4,color=pal26)

    ax[r].vlines(xvals,df.loc[group,'asymp.LCL'].values,df.loc[group,'asymp.UCL'].values,color='k',zorder=4,capstyle='round')
    # ax[r].hlines(df.loc[group,'asymp.LCL'].values,[i - .05 for i in xvals],[i + .05 for i in xvals],color=pal26,zorder=5)
    # ax[r].hlines(df.loc[group,'asymp.UCL'].values,[i - .05 for i in xvals],[i + .05 for i in xvals],color=pal26,zorder=6)

    # ax[r].hlines(0,*ax[r].get_xlim(),linestyle=':',color='black',zorder=0)
    ax[r].set_ylim(-.18,.35)
    ax[r].yaxis.set_major_locator(MultipleLocator(.15))   
    # sns.despine(left=True,ax=ax[r])
    ax[r].set_xticks(range(3))
    # ax[r].set_yticks([-.15,0,.15])
    ax[r].set_ylabel('Encoding-retrieval similarity',labelpad=1)
    ax[r].set_xticklabels(['Pre-\nconditioning','Conditioning','Extinction'],ha='center')
    ax[r].set_title(pfc_groups[group]['name'],pad=4)

    # if r == 0:
    #     ax[r].set_xlabel('')

    # else:
    #     # ax[r].set_ylabel('')
    #     # ax[r].set_yticks([])
    ax[r].set_xlim(-.7,ax[r].get_xlim()[1])
    ax[r].arrow(-.5,0,0,.025,head_width=0.06, head_length=0.01, fc='grey', ec='grey')
    ax[r].arrow(-.5,0,0,-.025,head_width=0.06, head_length=0.01, fc='grey', ec='grey')
    ax[r].text(-.5,.04,'more\nCS+',ha='center',va='bottom',fontsize=8)
    ax[r].text(-.5,-.04,'more\nCS-',ha='center',va='top',fontsize=8)

    top = df.loc[group,'asymp.UCL'].values.copy()
    bottom = df.loc[group,'asymp.LCL'].values.copy()
    if group == 'ptsd': top[1] = bottom[1] - .04
    for star, x, y in zip(pfc_groups[group]['stars'],xvals,top): ax[r].text(x,y,star,ha='center',va='bottom',fontsize=10)

legend_elements = [Patch(facecolor=pal2[0],edgecolor=pal2[0],label='dACC'),
                   Patch(facecolor=pal2[1],edgecolor=pal2[1],label='vmPFC')]
# fig.text(.025,.78,'Shown as CS+ - CS-',ha='left',va='bottom',fontsize=8)
fig.legend(handles=legend_elements,ncol=1,frameon=False,loc='upper left',bbox_to_anchor=(0.2,.975))

paired_barplot_annotate_brackets('**', 1, df.loc[('healthy','conditioning'),'asymp.UCL'].values, ax[0].get_ylim(), xtick_spread=.2, dh=.08, barh=.01, fs=10, ax=ax[0])
paired_barplot_annotate_brackets('*', 2, df.loc[('healthy','extinction'),'asymp.UCL'].values, ax[0].get_ylim(), xtick_spread=.2, dh=.08, barh=.01, fs=10, ax=ax[0])
paired_barplot_annotate_brackets('*', 1, df.loc[('ptsd','conditioning'),'asymp.UCL'].values, ax[1].get_ylim(), xtick_spread=.2, dh=.08, barh=.01, fs=10, ax=ax[1])

# ax[0].hlines(0,-.4,ax[0].get_xlim()[1],linestyle=':',color='black',zorder=0)
# ax[1].hlines(0,*ax[1].get_xlim(),linestyle=':',color='black',zorder=0)
plt.tight_layout()

'''colorbars for wholebrain searchlight'''
acq_norm = matplotlib.colors.Normalize(vmin=.05,vmax=.50)
ext_norm = matplotlib.colors.Normalize(vmin=.05,vmax=.3)

fig, ax = plt.subplots(figsize=mm2inch(50,3))
acq = np.array([[.087,.59]])
img = plt.imshow(acq, cmap="Reds")
plt.gca().set_visible(False)
cax = plt.axes([0.1, 0.2, 0.8, 0.6])
plt.colorbar(orientation="horizontal",ticks=[.087,.34,.59], cax=cax)

fig, ax = plt.subplots(figsize=mm2inch(50,3))
ext = np.array([[.1,.3]])
img = plt.imshow(ext, cmap="Blues")
plt.gca().set_visible(False)
cax = plt.axes([0.1, 0.2, 0.8, 0.6])
plt.colorbar(orientation="horizontal",ticks=[.1,.2,.3], cax=cax)

'''hippocampus phase stuff'''
df = pd.read_csv('hpc_phasedif_stats.csv')
df.roi = df.roi.apply(lambda x: x[:-4])
df = df[df.roi.isin(['hc_tail','hc_head'])]
df['roi'] = pd.Categorical(df.roi,categories=['hc_tail','hc_head'],ordered=True)
df = df.set_index(['group','phase','roi']).sort_index()

xvals = [(i-.2,i+.2) for i in range(3)]
xvals = [i for j in xvals for i in j]
pal2 = sns.color_palette(['slateblue','darkorange'],desat=1)
pal26 = pal2+pal2+pal2

fig, ax = plt.subplots(2,1,figsize=(mm2inch(88,150)),sharey=False)
for r, group in enumerate(['healthy','ptsd']):

    # ax[r].scatter(xvals,df.loc[group,'emmean'].values,c=pal26,s=8)
    # ax[r].plot(xvals[::2],df.loc[group,'emmean'].values[::2],color=pal2[0],linewidth=rcParams['lines.linewidth']*.75,zorder=1)
    # ax[r].plot(xvals[1::2],df.loc[group,'emmean'].values[1::2],color=pal2[1],linewidth=rcParams['lines.linewidth']*.75,zorder=2)

    ax[r].bar(xvals,df.loc[group,'emmean'].values,width=.4,color=pal26)
    ax[r].vlines(xvals,df.loc[group,'asymp.LCL'].values,df.loc[group,'asymp.UCL'].values,color='k',zorder=4,capstyle='round')
    # ax[r].hlines(df.loc[group,'asymp.LCL'].values,[i - .05 for i in xvals],[i + .05 for i in xvals],color=pal26,zorder=5)
    # ax[r].hlines(df.loc[group,'asymp.UCL'].values,[i - .05 for i in xvals],[i + .05 for i in xvals],color=pal26,zorder=6)

    # ax[r].hlines(0,*ax[r].get_xlim(),linestyle=':',color='black',zorder=0)
    ax[r].set_ylim(-.05,.14)
    ax[r].yaxis.set_major_locator(MultipleLocator(.05))   
    # sns.despine(left=True,ax=ax[r],trim=True)
    ax[r].set_xlabel('')
    ax[r].set_xticks(range(3))
    ax[r].set_xticklabels(['Pre-\nconditioning','Conditioning','Extinction'],ha='center')
    ax[r].set_ylabel('Encoding-retrieval similarity')
    # if r == 0:
    # else:
        # ax[r].set_ylabel('')
        # ax[r].set_yticks([])

    ax[r].set_title(pfc_groups[group]['name'])

#these are the cross ROI comps
paired_barplot_annotate_brackets('*', 1, df.loc[('healthy','conditioning'),'asymp.UCL'].values, ax[0].get_ylim(),xtick_spread=.2, dh=.025, barh=.01, fs=10, ax=ax[0])
paired_barplot_annotate_brackets('***', 2, df.loc[('healthy','extinction'),'asymp.UCL'].values, ax[0].get_ylim(),xtick_spread=.2, dh=.025, barh=.01, fs=10, ax=ax[0])
paired_barplot_annotate_brackets('*', 1, df.loc[('ptsd','conditioning'),'asymp.UCL'].values, ax[1].get_ylim(),xtick_spread=.2, dh=.025, barh=.01, fs=10, ax=ax[1])

# #within ROI comps
# paired_barplot_annotate_brackets('***',(1-.05,2-.05),df.loc[('healthy','extinction'),'asymp.UCL'].values,ax[0].get_ylim(), dh=.1, barh=.01, fs=10, ax=ax[0])
# paired_barplot_annotate_brackets('*',(1+.05,2+.05),df.loc[('healthy','extinction'),'asymp.UCL'].values,ax[0].get_ylim(), dh=.17, barh=.01, fs=10, ax=ax[0])

# paired_barplot_annotate_brackets('*',(1+.06,2+.05),df.loc[('ptsd','conditioning'),'asymp.UCL'].values,ax[1].get_ylim(), dh=.1, barh=.01, fs=10, ax=ax[1])
# paired_barplot_annotate_brackets('*',(0+.05,1+.04),df.loc[('ptsd','conditioning'),'asymp.UCL'].values,ax[1].get_ylim(), dh=.1, barh=.01, fs=10, ax=ax[1])


legend_elements = [Patch(facecolor=pal2[0],edgecolor=pal2[0],label='Posterior HC'),
                   Patch(facecolor=pal2[1],edgecolor=pal2[1],label='Anterior HC')]
fig.text(.06,.90,'Averaged across CS',ha='left',va='bottom',fontsize=8)
fig.legend(handles=legend_elements,ncol=1,frameon=False,loc='upper left',bbox_to_anchor=(0.025,1))
plt.tight_layout()


'''maybe some scatter lines of ers and pfc diff ers'''
df = pd.read_csv('all_data_lmm.csv').groupby(['group','phase','subject']).mean().reset_index()
df = df[df.phase.isin(['acquisition','extinction'])]

tail = pd.read_csv('hc_tail_pfc_diff_lines.csv').set_index(['group','phase']).sort_index()
head = pd.read_csv('hc_head_pfc_diff_lines.csv').set_index(['group','phase']).sort_index()

pal2 = sns.color_palette(['slategrey','indianred'],desat=1)
fig, ax = plt.subplots(1,2,figsize=mm2inch(150,90),sharey=True)
for g, group in enumerate(['healthy','ptsd']):
    ax[0].plot(tail.loc[(group,'acquisition'),'xvar'].values,tail.loc[(group,'acquisition'),'yvar'].values,color=pal2[g])
    # ax[0].fill_between(tail.loc[(group,'acquisition'),'xvar'].values,tail.loc[(group,'acquisition'),'LCL'].values, tail.loc[(group,'acquisition'),'UCL'].values, color=pal2[g], alpha=.2)

    ax[1].plot(head.loc[(group,'extinction'),'xvar'].values,head.loc[(group,'extinction'),'yvar'].values,color=pal2[g])
    # ax[1].fill_between(head.loc[(group,'extinction'),'xvar'].values,head.loc[(group,'extinction'),'LCL'].values, head.loc[(group,'extinction'),'UCL'].values, color=pal2[g], alpha=.2)


sns.scatterplot(data=df[df.phase=='acquisition'],x='hc_tail_ers',y='pfc_diff_ers',hue='group',ax=ax[0],palette=pal2)
sns.scatterplot(data=df[df.phase=='extinction'],x='hc_head_ers',y='pfc_diff_ers',hue='group',ax=ax[1],palette=pal2)

ax[0].legend_.remove();ax[1].legend_.remove()
ax[0].set_xlabel('Posterior HC ERS');ax[1].set_xlabel('Anterior HC ERS')
ax[0].set_ylabel('vmPFC - dACC ERS');ax[1].set_ylabel('') 
ax[0].set_title('Phase = Conditioning');ax[1].set_title('Phase = Extinction')
ax[0].set_xlim([-.15,.25]);ax[1].set_xlim([-.15,.35])
ax[0].yaxis.set_major_locator(MultipleLocator(.4));ax[1].yaxis.set_major_locator(MultipleLocator(.4))

legend_elements = [Patch(facecolor=pal2[0],edgecolor=pal2[0],label='Healthy'),
                   Patch(facecolor=pal2[1],edgecolor=pal2[1],label='PTSS')]
fig.text(.025,.78,'Averaged across CS',ha='left',va='bottom',fontsize=8)
fig.legend(handles=legend_elements,ncol=1,frameon=False,loc='upper left',bbox_to_anchor=(0.025,.95))
ax[0].spines['left'].set_bounds(ax[0].get_ylim()[0],.5)

plt.tight_layout()

'''a kde of all the vmPFC - dACC ERS data'''
df = pd.read_csv('all_data_lmm.csv')
df = df[df.phase.isin(['acquisition','extinction'])]
# df['dummy'] = ''
df = df.groupby(['group','phase','condition','subject']).mean().reset_index()
pal = sns.color_palette(['red','red','dodgerblue','dodgerblue'],desat=.8)
lss = ['--','-','--','-']
ass = [.5,1,.5,1]
zss = [1,3,2,4]
df['phase_con'] = df.phase + '_' + df.condition


fig = plt.figure(figsize=mm2inch(120,120),constrained_layout=False) 
gs = gridspec.GridSpec(2, 1,height_ratios=[1, 3]) 
ax = plt.subplot(gs[0])

# fig, (ax, ax1) = plt.subplots(2,1,figsize=mm2inch(120,40))
sns.kdeplot(data=df,x='pfc_diff_ers',hue='phase_con',palette=pal,ax=ax,cut=0)
for l, line in enumerate(ax.get_lines()):
    line.set_linestyle(lss[l])
    line.set_alpha(ass[l])
    line.set_zorder(zss[l])
    line.set_c(line.get_c())
ax.legend_.remove()
sns.despine(ax=ax,left=True)
ax.vlines(0,*ax.get_ylim(),linestyle=':',color='black',zorder=0)

ax.arrow(0,.05,.08,0,head_width=0.05, head_length=0.03, fc='grey', ec='grey')
ax.arrow(0,.05,-.08,0,head_width=0.05, head_length=0.03, fc='grey', ec='grey')
ax.text(.12,.05,'more\nvmPFC',ha='left',va='center',fontsize=8)
ax.text(-.12,.05,'more\ndACC',ha='right',va='center',fontsize=8)
ax.set_xlabel('mPFC difference in ERS')


ers = pd.read_csv('ers_slopes_pfc_pred.csv').drop(columns='Unnamed: 0').rename(index={0:'est',1:'LCL',2:'UCL'}).swapaxes('index','columns').reset_index()
ers['stat'] = 'ers'

yvals = [i for i in range(5)]
ax1 = plt.subplot(gs[1])
ax1.scatter(ers.est.values,yvals,s=8,c='k')
ax1.hlines(yvals,ers.LCL.values,ers.UCL.values,color='k',capstyle='round')
ax1.set_xlabel('Predictiveness on\nmPFC difference in ERS')

ax1.set_yticks(range(5))
sns.despine(ax=ax1,left=True)
ax1.set_yticklabels(['Posterior HC','HC Body','Anterior HC','BLA','CeM'],ha='right')

align.xaxes(ax,0,ax1,0)
ers_stars = [('*',-1),('',1),('***',1),('',1),('**',-1)]


for s, pair in enumerate(ers_stars):
    if pair[1] == -1:
        ax1.text(ers.LCL.values[s]-.01,yvals[s]-.18,pair[0],ha='right',va='bottom',fontsize=10)        
    else:
        ax1.text(ers.UCL.values[s]+.01,yvals[s]-.18,pair[0],ha='left',va='bottom',fontsize=10)        

ax1.set_ylim(-.7,ax1.get_ylim()[1])
ax1.vlines(0,*ax1.get_ylim(),linestyle=':',color='black',zorder=0)
ax1.arrow(0,-.45,0.012,0,head_width=0.15, head_length=0.0045, fc='grey', ec='grey')
ax1.arrow(0,-.45,-0.012,0,head_width=0.15, head_length=0.0045, fc='grey', ec='grey')
ax1.text(.018,-.45,'predicts more\nvmPFC',ha='left',va='center',fontsize=8)
ax1.text(-.018,-.45,'predicts more\ndACC',ha='right',va='center',fontsize=8)

ax1.set_ylabel('Subcortical ERS')

legend_elements = [Line2D([0],[0],c=pal[0],linestyle='-',label='Cond.'),
                   Line2D([0],[0],c=pal[2],linestyle='-',label='Ext.'),
                   Line2D([0],[0],c='k',linestyle='-',label='CS+'),
                   Line2D([0],[0],c='k',linestyle='--',label='CS-')]
plt.sca(ax)
legend1 = plt.legend(handles=legend_elements[:2],ncol=1,frameon=False,loc='upper left',bbox_to_anchor=(0,1),numpoints=2,handletextpad=0.1)
ax = plt.gca().add_artist(legend1)
# Create another legend for the second line.
plt.legend(handles=legend_elements[2:], ncol=1,frameon=False,loc='upper right',bbox_to_anchor=(1,1),numpoints=2,handletextpad=0.1)

plt.tight_layout()



'''just a quick plot to show the strength of the estiamtes for pfc diff'''
uni = pd.read_csv('uni_slopes_pfc_pred.csv').drop(columns='Unnamed: 0').rename(index={0:'est',1:'LCL',2:'UCL'}).swapaxes('index','columns')
uni['stat'] = 'uni'

ers = pd.read_csv('ers_slopes_pfc_pred.csv').drop(columns='Unnamed: 0').rename(index={0:'est',1:'LCL',2:'UCL'}).swapaxes('index','columns')
ers['stat'] = 'ers'

df = pd.concat((uni,ers)).reset_index().rename(columns={'index':'roi'})

xvals_top = [i+.1 for i in range(5)]
xvals_bot = [i-.1 for i in range(5)]
pal2 = sns.color_palette(['slategrey','indianred'],desat=1)
ers_stars = [('*',-1),('',1),('***',1),('',1),('**',-1)]
uni_stars = ['***','***','***','***','***']


fig, ax1 = plt.subplots(figsize=mm2inch(120,90))
ax1.scatter(xvals_bot,ers.est.values,s=8,c=pal2[1])
ax1.vlines(xvals_bot,ers.LCL.values,ers.UCL.values,color=pal2[1],capstyle='round')
ax1.tick_params(axis='y', labelcolor=pal2[1])
ax1.set_ylabel('Slope (ERS)',color=pal2[1])
ax2 = ax1.twinx()
ax2.scatter(xvals_top,uni.est.values,s=8,c=pal2[0])
ax2.vlines(xvals_top,uni.LCL.values,uni.UCL.values,color=pal2[0],capstyle='round')
ax2.tick_params(axis='y', labelcolor=pal2[0])
ax2.set_ylabel('Slope (univariate)',rotation=270,color=pal2[0])
sns.despine(left=True,trim=True,ax=ax1)
sns.despine(left=True,trim=True,ax=ax2)
ax2.set_yticks([0,-.0015,-.003])
ax1.set_yticks([-.1,0,.1])
ax1.set_xticks(range(5))
ax1.set_xticklabels(['Posterior HC','HC Body','Anterior HC','BLA','CeM'],ha='center',rotation=45)

ax1.yaxis.set_label_coords(-0.15,.725)
ax2.yaxis.set_label_coords(1.3,.725) 
for s, star in enumerate(uni_stars): ax2.text(xvals_top[s],uni.LCL.values[s]-.0001,star,ha='center',va='top',fontsize=10)
for s, pair in enumerate(ers_stars):
    if pair[1] == -1:
        ax1.text(xvals_bot[s],ers.LCL.values[s]-.01,pair[0],ha='center',va='top',fontsize=10)        
    else:
        ax1.text(xvals_bot[s],ers.UCL.values[s],pair[0],ha='center',va='bottom',fontsize=10)        

ax1.set_xlim(-0.4,5.2)
ax2.set_ylim(-.0035,ax2.get_ylim()[1])
align.yaxes(ax1,0,ax2,0)
ax1.hlines(0,*ax1.get_xlim(),linestyle=':',color='black',zorder=0)

ax1.arrow(5,0,0,.025,head_width=0.1, head_length=0.015, fc='grey', ec='grey')
ax1.arrow(5,0,0,-.025,head_width=0.1, head_length=0.015, fc='grey', ec='grey')
ax1.text(5,.045,'predicts\nvmPFC',ha='center',va='bottom',fontsize=8)
ax1.text(5,-.045,'predicts\ndACC',ha='center',va='top',fontsize=8)

plt.tight_layout()


'''maybe a scatter lines of ev and pfc diff ers'''
df = pd.read_csv('all_data_lmm.csv')
df = df[df.phase.isin(['acquisition','extinction'])]
df = df.groupby(['condition','subject']).mean().reset_index()

ev = pd.read_csv('ev_pfc_diff_lines.csv').set_index(['condition']).sort_index().sort_values(by='xvar')

fig, ax = plt.subplots(figsize=mm2inch(100,90))

ax.plot(ev.loc['CS+','xvar'].values,ev.loc['CS+','yvar'].values,color=cpal[0])
ax.fill_between(ev.loc['CS+','xvar'].values,ev.loc['CS+','LCL'].values, ev.loc['CS+','UCL'].values, color=cpal[0], alpha=.2)

ax.plot(ev.loc['CS-','xvar'].values,ev.loc['CS-','yvar'].values,color=cpal[1])
ax.fill_between(ev.loc['CS-','xvar'].values,ev.loc['CS-','LCL'].values, ev.loc['CS-','UCL'].values, color=cpal[1], alpha=.2)
sns.scatterplot(data=df,x='ev',y='pfc_diff_ers',hue='condition',ax=ax,palette=cpal)

ax.legend_.remove();
ax.set_xlabel('Extinction context evidence')
ax.set_ylabel('vmPFC - dACC ERS')
ax.set_title('Reinstated extinction context\npredicts ERS bias to vmPFC',pad=5)
ax.yaxis.set_major_locator(MultipleLocator(.4))
ax.xaxis.set_major_locator(MultipleLocator(.5))
# plt.tight_layout()

legend_elements = [Patch(facecolor=cpal[0],edgecolor=cpal[0],label='CS+'),
                   Patch(facecolor=cpal[1],edgecolor=cpal[1],label='CS-')]
ax.text(0,-.9,'Averaged across encode phase',ha='left',va='bottom',fontsize=8)
fig.legend(handles=legend_elements,ncol=1,frameon=False,loc='upper left',bbox_to_anchor=(0.04,.95))
ax.spines['left'].set_bounds(ax.get_ylim()[0],.5)

plt.tight_layout()

'''just a small horizontal plot of the context evidence stuff'''
df = pd.read_csv('ev_slopes_pfc_pred.csv')
fig, ax = plt.subplots(figsize=mm2inch(100,30))
ax.scatter(df['ev.trend'],[.15,.25],s=8,c='k')
ax.hlines([.15],df['asymp.LCL'].values[0],df['asymp.UCL'].values[0],color='k',capstyle='round',linestyle='--')
ax.hlines([.25],df['asymp.LCL'].values[1],df['asymp.UCL'].values[1],color='k',capstyle='round')
ax.set_xlabel('Predictiveness on\nmPFC difference in ERS')
ax.set_ylim((0,.35))
ax.vlines(0,*ax.get_ylim(),linestyle=':',color='black',zorder=0)
ax.arrow(0,.05,0.012,0,head_width=0.025, head_length=0.0045, fc='grey', ec='grey')
ax.arrow(0,.05,-0.012,0,head_width=0.025, head_length=0.0045, fc='grey', ec='grey')
ax.text(.018,.05,'predicts more\nvmPFC',ha='left',va='center',fontsize=8)
ax.text(-.018,.05,'predicts more\ndACC',ha='right',va='center',fontsize=8)
ax.set_xlim((-.125,ax.get_xlim()[1]))

legend_elements = [Line2D([0],[0],c='k',linestyle='-',label='CS+'),
                   Line2D([0],[0],c='k',linestyle='--',label='CS-')]
ax.legend(handles=legend_elements, ncol=1,frameon=False,loc='upper left',bbox_to_anchor=(0,1),numpoints=2,handletextpad=0.1)
sns.despine(ax=ax,left=True)
ax.set_yticks([])
ax.text(.215,.22,'***',ha='left',va='bottom')
plt.tight_layout()

'''day 1 behavioral results'''
df = pd.read_csv('day1_behavior.csv')
df = df[df.phase.isin(['fear','late_ext'])]

pal = sns.color_palette(['slategrey',wes_palettes['Chevalier'][0]],desat=None)

fig, (ax1, ax2) = plt.subplots(1,2,figsize=mm2inch(120,60))
sns.barplot(data=df,x='phase',y='scr',hue='group',ax=ax1,palette=pal,saturation=1)
sns.barplot(data=df,x='phase',y='exp',hue='group',ax=ax2,palette=pal,saturation=1)
align_yaxis(ax1,0,ax2,0,-1,.6)
ax1.set_ylabel('Sqrt( SCR )')
ax2.set_ylabel('Expectancy')

ax1.set_yticks([0,.1,.2,.3])
ax2.yaxis.set_major_locator(MultipleLocator(.2))   

for ax in [ax1,ax2]:
    ax.set_xticklabels(['Conditioning','Late\nExtinction'])
    ax.legend_.remove()
    ax.set_xlim(-.77,ax.get_xlim()[1])
    y2, y1 = ax.get_ylim()
    ax.arrow(-.57,0,0,y1*.03,head_width=.06, head_length=y1*0.02, fc='grey', ec='grey')
    ax.arrow(-.57,0,0,y1*-.03,head_width=.06, head_length=y1*0.02, fc='grey', ec='grey')
    ax.text(-.57,y1*.055,'more\nCS+',ha='center',va='bottom',fontsize=8)
    ax.text(-.57,y1*-.05,'more\nCS-',ha='center',va='top',fontsize=8)
    ax.set_xlabel('')
    ax.set_ylim(ax.get_ylim()[0]*1.27,ax.get_ylim()[1])

ax1.set_title('Autonomic arousal',pad=.1)
ax2.set_title('Shock expectancy',pad=.1)

legend_elements = [Patch(facecolor=pal[0],edgecolor=pal[0],label='Healthy'),
                   Patch(facecolor=pal[1],edgecolor=pal[1],label='PTSS')]
fig.legend(handles=legend_elements,ncol=1,frameon=False,loc='upper left',bbox_to_anchor=(0,.975))

ax1.spines['left'].set_bounds(ax1.get_ylim()[0], .32)

plt.tight_layout(pad=.2)


'''trying to show the hippocampal univariate stuff'''
tail = pd.read_csv('hc_tail_slopes_con_phase.csv').rename(columns={'hc_tail_ret_uni.trend':'est'})
tail['roi'] = 'hc_tail'
body = pd.read_csv('hc_body_slopes_con_phase.csv').rename(columns={'hc_body_ret_uni.trend':'est'})
body['roi'] = 'hc_body'
df = pd.concat((tail,body))
yvals = [.1,.125,.15,.175,.3,.325,.35,.375]
fig, ax = plt.subplots(figsize=mm2inch(100,30))
ax.scatter(df['est'],yvals,s=8,c='k')
ax.hlines(yvals,df['asymp.LCL'].values,df['asymp.UCL'].values,color='k',capstyle='round')
ax.hlines([.25],df['asymp.LCL'].values[1],df['asymp.UCL'].values[1],color='k',capstyle='round')
ax.set_xlabel('Predictiveness on\nmPFC difference in ERS')
ax.set_ylim((0,.35))
ax.vlines(0,*ax.get_ylim(),linestyle=':',color='black',zorder=0)
ax.arrow(0,.05,0.012,0,head_width=0.025, head_length=0.0045, fc='grey', ec='grey')
ax.arrow(0,.05,-0.012,0,head_width=0.025, head_length=0.0045, fc='grey', ec='grey')
ax.text(.018,.05,'predicts more\nvmPFC',ha='left',va='center',fontsize=8)
ax.text(-.018,.05,'predicts more\ndACC',ha='right',va='center',fontsize=8)
ax.set_xlim((-.125,ax.get_xlim()[1]))

legend_elements = [Line2D([0],[0],c='k',linestyle='-',label='CS+'),
                   Line2D([0],[0],c='k',linestyle='--',label='CS-')]
ax.legend(handles=legend_elements, ncol=1,frameon=False,loc='upper left',bbox_to_anchor=(0,1),numpoints=2,handletextpad=0.1)
sns.despine(ax=ax,left=True)
ax.set_yticks([])
ax.text(.215,.22,'***',ha='left',va='bottom')
plt.tight_layout()
