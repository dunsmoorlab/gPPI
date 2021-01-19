from fg_config import *
phase3 = ['baseline','acquisition','extinction']
subcort = ['hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']

pfc = pd.read_csv('pfc_ers_cleaned.csv').set_index(['roi','condition','phase','group','subject'])
pfc = pfc.loc['dACC'] - pfc.loc['vmPFC']
pfc = pfc.loc['CS+'] - pfc.loc['CS-']
pfc = pfc.reset_index()
pfc = pfc[pfc.phase != 'baseline']

sns.catplot(data=pfc,x='group',y='rsa',hue='phase',palette=phase_pal[1:],kind='bar')

# #PFC ERS
# fig, ax = plt.subplots()
# sns.barplot(data=pfc,x='group',y='rsa',hue='phase',palette=phase_pal,ax=ax)
# ax.set_title('PFC encoding-retrieval similarity\ndorsal-ventral index')
# ax.set_ylabel('dACC - vmPFC (CS+ - CS-)')


#connectivity
'''prep connectivity for lmm in r just focusing on PFC rois'''
conn = pd.read_csv('cleaned_gPPI.csv').set_index(['target','group','phase','condition','subject','seed'])
# conn = conn.loc['vmPFC']
conn = conn.loc['vmPFC'] - conn.loc['amyg_cem']
#conn = conn.loc['CS+'] - conn.loc['CS-']

conn = conn.unstack(level='seed')

conn.columns = conn.columns.droplevel(0)
conn = conn.rename_axis(None, axis=1)
ers = pd.read_csv('ers_cleaned_lmm.csv').drop(columns=['stimulus']).groupby(['group','phase','condition','subject']).mean()
conn = pd.concat((conn,ers),axis=1)

conn = conn.drop('baseline',level='phase')

conn.to_csv('cleaned_gPPI_lmm.csv')
conn = conn.reset_index()
# sns.catplot(data=conn,x='group',y='conn',hue='phase',col='seed',palette=phase_pal,kind='bar')
conn['phase_con'] = conn.phase + '_' + conn.condition

_rois = ['hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']
fig, ax = plt.subplots(1,len(_rois),sharey=True)
for r, roi in enumerate(_rois):    
    sns.barplot(data=conn,x='group',y=f'{roi}',hue='phase_con',palette=paired_pal,ax=ax[r])
    # sns.stripplot(data=df,x='group',y=f'{roi}_ers',hue='phase_con',dodge=True,color='black',size=2,ax=ax[r])
    ax[r].legend(loc='center left',bbox_to_anchor=(1,.5)) if r == len(_rois) - 1 else ax[r].legend_.remove()
    ax[r].set_ylabel('Connectivity') if r == 0 else ax[r].set_ylabel('')
    ax[r].set_title(f'ROI = {roi}')
    plt.suptitle('vmPFC - amyg_cem connectivity')



#correlating this same 



'''correlating univariate to gPPI the same way we did uni-ERS'''
conn = pd.read_csv('cleaned_gPPI.csv').set_index(['condition','group','phase','seed','target','subject']).sort_index()
df = pd.read_csv('all_data.csv').set_index(['group','phase','subject']).sort_index()
conn = conn.loc['CS+'] - conn.loc['CS-']
for seed in subcort:
    for target in ['vmPFC','dACC']:
        for phase in phase3:
            print(seed,target,phase)
            c = skipped_corr(conn.loc[('healthy',phase,seed,target),'conn'],pfc.loc[('healthy',phase),'rsa'], return_dist = True)
            p = skipped_corr(conn.loc[('ptsd',phase,seed,target),'conn'],pfc.loc[('ptsd',phase),'rsa'], return_dist = True)
            diff_p = np.min(((1 - np.mean((c - p) > 0)) * 2,(1 - np.mean((p - c) > 0)) * 2))
            print(f'difference P = {diff_p}')
            print('\n\n')
    input()



'''first mediation'''
df = pd.read_csv('all_data.csv').set_index(['group','phase','subject']).sort_index()
# df['bgroup'] = df.group.apply(lambda x: 1 if x=='healthy' else 0)
c = df.loc[('healthy','extinction')].copy()
pg.mediation_analysis(data=c,x='hc_head-ret-diff_uni',m='amyg_bla-ret-diff_uni',y='vmPFC-diff_ers',n_boot=10000)

'''ctx reinstatement'''
# ev = pd.read_csv('ctx_ev.csv').set_index(['group','phase','condition','subject']).sort_index()
ev = pd.read_csv('ctx_ev.csv').set_index(['condition','group','phase','subject']).sort_index()
ev = ev.loc['CS+'] - ev.loc['CS-']
# ev = pd.read_csv('ctx_ev.csv').groupby(['group','condition','subject']).mean()


thing = 'diff_ers'
for roi in ['vmPFC','dACC']:
# for roi in rois:
    for phase in phase3:
        print(roi,phase)
        c = skipped_corr(ev.loc[('healthy','extinction'),'proba'],df.loc[('healthy',phase),f'{roi}-{thing}'], return_dist = True)
        p = skipped_corr(ev.loc[('ptsd','extinction'),'proba'],df.loc[('ptsd',phase),f'{roi}-{thing}'], return_dist = True)
        diff_p = np.min(((1 - np.mean((c - p) > 0)) * 2,(1 - np.mean((p - c) > 0)) * 2))
        print(f'difference P = {np.round(diff_p,4)}')
        print('\n\n')
    input()



'''getting together a single dataframe with trial-wise data for LMM'''
ev = pd.read_csv('ctx_ev_lmm.csv').set_index(['group','phase','condition','subject','stimulus']).sort_index()
ret_uni = pd.read_csv('ret_uni_lmm.csv').set_index(['group','phase','condition','subject','stimulus']).sort_index()
ers = pd.read_csv('ers_cleaned_lmm.csv').set_index(['group','phase','condition','subject','stimulus']).sort_index()
df = pd.concat([ers,ret_uni,ev],axis=1)
df.to_csv('all_data_lmm.csv')

df = pd.read_csv('all_data_lmm.csv').set_index(['group','phase','condition','subject','stimulus']).sort_index()
df = df.loc[(['healthy','ptsd'],['acquisition','extinction']),].reset_index().drop(columns='stimulus')
for col in ['amyg_bla_ers', 'amyg_cem_ers', 'dACC_ers', 'hc_body_ers','hc_head_ers', 'hc_tail_ers', 'vmPFC_ers', 'pfc_diff_ers','vmPFC_ret_uni', 'dACC_ret_uni', 'hc_tail_ret_uni', 'hc_body_ret_uni','hc_head_ret_uni', 'amyg_bla_ret_uni', 'amyg_cem_ret_uni', 'ev']:
    # df[col] = df.groupby(['subject'])[col].transform(lambda x: zscore(x,ddof=1))
    df[col] = zscore(df[col],ddof=1)
df.to_csv('emo_data_lmm_zscore.csv')


paired_pal = [phase_pal[1],[i*1.75 for i in phase_pal[1]], phase_pal[2], [i*1.75 for i in phase_pal[2]]]
'''hippocampus is key'''
df = pd.read_csv('all_data_lmm.csv').set_index(['group','phase','condition','subject','stimulus']).sort_index()
df = df.loc[(['healthy','ptsd'],['acquisition','extinction']),].reset_index().drop(columns='stimulus')
df = df.groupby(['group','phase','condition','subject']).mean().reset_index()
df['phase_con'] = df.phase + '_' + df.condition

_rois = ['dACC','vmPFC']
fig, ax = plt.subplots(1,len(_rois),sharey=True)
for r, roi in enumerate(_rois):    
    sns.barplot(data=df,x='group',y=f'{roi}_ers',hue='phase_con',palette=paired_pal,ax=ax[r])
    # sns.stripplot(data=df,x='group',y=f'{roi}_ers',hue='phase_con',dodge=True,color='black',size=2,ax=ax[r])
    ax[r].legend(loc='center left',bbox_to_anchor=(1,.5)) if r == len(_rois) - 1 else ax[r].legend_.remove()
    ax[r].set_ylabel('Encoding-retrieval similarity') if r == 0 else ax[r].set_ylabel('')
    ax[r].set_title(f'ROI = {roi}')

'''trying to vis the relationships'''
df = pd.read_csv('all_data_lmm.csv').set_index(['group','phase','condition','subject','stimulus']).sort_index()
df = df.loc[(['healthy','ptsd'],['acquisition','extinction']),].reset_index().drop(columns='stimulus')
df = df.groupby(['group','phase','condition','subject']).mean().reset_index()
df['phase_con'] = df.phase + '_' + df.condition

ev_beta = 9.331e-02
ev_intercept = -7.320e-02


fig, ax = plt.subplots()
sns.scatterplot(data=df,x='ev',y='pfc_diff_ers',hue='phase_con',style='group',palette=paired_pal,ax=ax)
ax.legend(loc='center left',bbox_to_anchor=(1,.5))
# ax.scatter(df.groupby('subject').mean()['ev'],df.groupby('subject').mean()['pfc_diff_ers'])
ax.plot(np.linspace(0,1,100),(np.linspace(0,1,100)*ev_beta)+ev_intercept,color='black')
# ax.set_ylabel('pfc_diff_ers');ax.set_xlabel('ev')


hc_diff_beta = 8.647e-02
hc_intercept = -2.040e-02
fig, ax = plt.subplots()
sns.scatterplot(data=df,x='hc_diff_ers',y='pfc_diff_ers',hue='phase_con',style='group',palette=paired_pal,ax=ax)
ax.legend(loc='center left',bbox_to_anchor=(1,.5))
ax.plot(np.linspace(df.hc_diff_ers.min(),df.hc_diff_ers.max(),100),(np.linspace(df.hc_diff_ers.min(),df.hc_diff_ers.max(),100)*hc_diff_beta)+hc_intercept,color='black')

