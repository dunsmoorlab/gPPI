from fg_config import *
phase3 = ['baseline','acquisition','extinction']
subcort = ['hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']

pfc = pd.read_csv('pfc_ers_cleaned.csv').set_index(['roi','condition','phase','group','subject'])
pfc = pfc.loc['dACC'] - pfc.loc['vmPFC']
pfc = pfc.loc['CS+'] - pfc.loc['CS-']

pfc = pfc.reset_index()
# sns.catplot(data=pfc,x='group',y='rsa',hue='phase',palette=phase_pal,kind='bar')

# #PFC ERS
# fig, ax = plt.subplots()
# sns.barplot(data=pfc,x='group',y='rsa',hue='phase',palette=phase_pal,ax=ax)
# ax.set_title('PFC encoding-retrieval similarity\ndorsal-ventral index')
# ax.set_ylabel('dACC - vmPFC (CS+ - CS-)')


#connectivity
conn = pd.read_csv('cleaned_gPPI.csv').set_index(['target','condition','phase','seed','group','subject'])
conn = conn.loc['amyg_cem'] - conn.loc['vmPFC']
conn = conn.loc['CS+'] - conn.loc['CS-']

conn = conn.reset_index()
# sns.catplot(data=conn,x='group',y='conn',hue='phase',col='seed',palette=phase_pal,kind='bar')

#correlating this same 



'''correlating univariate to gPPI the same way we did uni-ERS'''
conn = pd.read_csv('cleaned_gPPI.csv').set_index(['condition','group','phase','seed','target','subject']).sort_index()
df = pd.read_csv('all_data.csv').set_index(['group','phase','subject']).sort_index()
conn = conn.loc['CS+'] - conn.loc['CS-']
for seed in subcort:
    for target in ['vmPFC','dACC']:
        for phase in phase3:
            print(seed,target,phase)
            c = skipped_corr(conn.loc[('healthy',phase,seed,target),'conn'],df.loc[('healthy',phase),f'{target}-diff_ers'], return_dist = True)
            p = skipped_corr(conn.loc[('ptsd',phase,seed,target),'conn'],df.loc[('ptsd',phase),f'{target}-diff_ers'], return_dist = True)
            diff_p = np.min(((1 - np.mean((c - p) > 0)) * 2,(1 - np.mean((p - c) > 0)) * 2))
            print(f'difference P = {diff_p}')
            print('\n\n')
    input()
