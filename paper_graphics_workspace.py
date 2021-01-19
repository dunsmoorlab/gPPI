#looking at pfc ers all phases with bars
paired_pal = [phase_pal[0], [i*1.75 for i in phase_pal[0]], phase_pal[1],[i*1.75 for i in phase_pal[1]], phase_pal[2], [i*1.75 for i in phase_pal[2]]]
'''hippocampus is key'''
df = pd.read_csv('all_data_lmm.csv').set_index(['group','phase','condition','subject','stimulus']).sort_index()
df = df.loc[(['healthy','ptsd'],['acquisition','baseline','extinction']),].reset_index().drop(columns='stimulus')
df = df.groupby(['group','phase','condition','subject']).mean().reset_index()
df['phase_con'] = df.phase + '_' + df.condition

_rois = ['dACC','vmPFC']
fig, ax = plt.subplots(1,len(_rois),sharey=True)
for r, roi in enumerate(_rois):    
    sns.barplot(data=df,x='group',y=f'{roi}_ers',hue='phase_con',palette=paired_pal,ax=ax[r],hue_order=['baseline_CS+','baseline_CS-','acquisition_CS+','acquisition_CS-','extinction_CS+','extinction_CS-'])
    # sns.stripplot(data=df,x='group',y=f'{roi}_ers',hue='phase_con',dodge=True,color='black',size=2,ax=ax[r])
    ax[r].legend(loc='center left',bbox_to_anchor=(1,.5)) if r == len(_rois) - 1 else ax[r].legend_.remove()
    ax[r].set_ylabel('Encoding-retrieval similarity') if r == 0 else ax[r].set_ylabel('')
    ax[r].set_title(f'ROI = {roi}')