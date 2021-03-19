from fg_config import *
from paper_graphics_style import *

df = pd.read_csv('conn_stats_lmm.csv')
df = df[df.phase.isin(['acquisition','extinction'])]
df.phase = df.phase.apply(lambda x: 'conditioning' if x == 'acquisition' else x)
df['phase_con'] = df.condition + '_' + df.phase

pdf = df.groupby(['group','phase','seed','target','subject']).mean().reset_index()
ncdf = df.groupby(['group','seed','target','subject']).mean().reset_index()

animal = df[df.csp_cat == 'animal'][df.seed == 'animal'][df.condition == 'CS+']
tool = df[df.csp_cat == 'tool'][df.seed == 'tool'][df.condition == 'CS+']
cat = pd.concat((animal,tool))
cat = cat[cat.phase.isin(['acquisition','extinction'])][cat.target.isin(['dACC','vmPFC'])]
cat = cat.set_index(['group','phase','seed','csp_cat','target','condition','subject']).sort_index().reset_index()


pfc_pal = sns.color_palette(['#A83465','#22566B'],desat=1)
phase_pal = sns.color_palette(['red','dodgerblue'],desat=.8)
paired_pal = sns.color_palette(['red','lightcoral','dodgerblue','powderblue'],desat=.8)

'''hc head to dACC'''
fig, ax = plt.subplots()
sns.barplot(data=pdf[pdf.seed=='hc_head'][pdf.target=='dACC'],x='group',y='conn',hue='phase',palette=phase_pal,ax=ax)
ax.set_title('aHPC to dACC',fontsize=16)
paired_barplot_annotate_brackets('*', 1, [.4,.4], ax.get_ylim(), xtick_spread=.2, dh=.08, barh=.01, fs=15, ax=ax)
plt.tight_layout()

'''hc head to vmPFC'''
fig, ax = plt.subplots()
sns.barplot(data=pdf[pdf.seed=='hc_head'][pdf.target=='vmPFC'],x='group',y='conn',hue='phase',palette=phase_pal,ax=ax)
ax.set_title('aHPC to vmPFC',fontsize=16)
paired_barplot_annotate_brackets('*', 1, [.2,.2], ax.get_ylim(), xtick_spread=.2, dh=.08, barh=.01, fs=15, ax=ax)
plt.tight_layout()

'''vmPFC to BLA'''
fig, ax = plt.subplots()
sns.barplot(data=df[df.seed=='vmPFC'][df.target=='amyg_bla'],x='group',y='conn',hue='phase_con',palette=paired_pal,ax=ax)
ax.plot([.1,.3],[.33,.33],color='k')
ax.text(.2,.33,'*',fontsize=15,ha='center')
ax.set_title('vmPFC to BLA')
plt.tight_layout()

'''BLA to CEM'''
fig, ax = plt.subplots()
sns.barplot(data=df[df.seed=='amyg_bla'][df.target=='amyg_cem'],x='group',y='conn',hue='phase_con',palette=paired_pal,ax=ax)
ax.plot([-.3,-.1],[.15,.15],color='k')
ax.text(-.2,.17,'<.05 unc.\n~',fontsize=15,ha='center')
ax.set_title('BLA to CEM')
plt.tight_layout()

'''hc head to PFC'''
fig, ax = plt.subplots()
sns.barplot(data=ncdf[ncdf.seed=='hc_head'][ncdf.target.isin(['dACC','vmPFC'])],x='group',y='conn',hue='target',palette=pfc_pal,ax=ax)
ax.set_title('aHPC to mPFC',fontsize=16)
paired_barplot_annotate_brackets('**', 1, [.3,.3], ax.get_ylim(), xtick_spread=.2, dh=.08, barh=.01, fs=15, ax=ax)
plt.tight_layout()

'''hc tail to PFC'''
fig, ax = plt.subplots()
sns.barplot(data=ncdf[ncdf.seed=='hc_tail'][ncdf.target.isin(['dACC','vmPFC'])],x='group',y='conn',hue='target',palette=pfc_pal,ax=ax)
ax.set_title('pHPC to mPFC',fontsize=16)
paired_barplot_annotate_brackets('**', 1, [.25,.25], ax.get_ylim(), xtick_spread=.2, dh=.08, barh=.01, fs=15, ax=ax)
plt.tight_layout()

'''category cortex to vmPFC'''
for group in ['healthy','ptsd']:
    fig, ax = plt.subplots(1,2, sharey=True)
    for t, target in enumerate(['dACC','vmPFC']):
        sns.barplot(data=cat[cat.group == group][cat.target==target], x='seed', y='conn', hue='phase',palette=phase_pal, ax=ax[t])
        ax[t].set_title(f'target = {target}',fontsize=12)
        if group == 'healthy' and target == 'vmPFC':
            paired_barplot_annotate_brackets('*', 0, [.25,.25], ax[t].get_ylim(), xtick_spread=.2, dh=.08, barh=.01, fs=15, ax=ax[t])
        if group == 'healthy' and target == 'dACC':
            paired_barplot_annotate_brackets('~', 0, [.4,.4], ax[t].get_ylim(), xtick_spread=.2, dh=.08, barh=.01, fs=15, ax=ax[t])
    plt.suptitle(group+' CS+',fontsize=15)
    ax[0].legend_.remove()
    plt.tight_layout()    