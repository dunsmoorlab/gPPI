from fg_config import *

df = pd.read_csv('all_data_lmm.csv').set_index(['subject','stimulus'])

df['trial_num'] = 0
for sub in all_sub_args:
    df.loc[sub,'trial_num'] = bids_meta(sub).mem_df.reset_index().set_index(['stimulus']).sort_index().dropna(subset=['response'])['index'].values

df = df.reset_index().set_index(['subject','phase','condition','trial_num']).sort_index()
df = df[['vmPFC_ret_uni','dACC_ret_uni','hc_tail_ret_uni','hc_body_ret_uni','hc_head_ret_uni','amyg_bla_ret_uni','amyg_cem_ret_uni']]


vmPFC = df.groupby(['subject','phase','condition']).apply(lambda x: x.corrwith(x.vmPFC_ret_uni,axis=0,method='pearson')).drop(columns='vmPFC_ret_uni').apply(np.arctanh)
dACC = df.groupby(['subject','phase','condition']).apply(lambda x: x.corrwith(x.dACC_ret_uni,axis=0,method='pearson')).drop(columns='dACC_ret_uni').apply(np.arctanh)
hc_tail = df.groupby(['subject','phase','condition']).apply(lambda x: x.corrwith(x.dACC_ret_uni,axis=0,method='pearson')).drop(columns='hc_tail_ret_uni').apply(np.arctanh)
hc_body = df.groupby(['subject','phase','condition']).apply(lambda x: x.corrwith(x.dACC_ret_uni,axis=0,method='pearson')).drop(columns='hc_body_ret_uni').apply(np.arctanh)
hc_head = df.groupby(['subject','phase','condition']).apply(lambda x: x.corrwith(x.dACC_ret_uni,axis=0,method='pearson')).drop(columns='hc_head_ret_uni').apply(np.arctanh)
bla = df.groupby(['subject','phase','condition']).apply(lambda x: x.corrwith(x.dACC_ret_uni,axis=0,method='pearson')).drop(columns='amyg_bla_ret_uni').apply(np.arctanh)
cem = df.groupby(['subject','phase','condition']).apply(lambda x: x.corrwith(x.dACC_ret_uni,axis=0,method='pearson')).drop(columns='amyg_cem_ret_uni').apply(np.arctanh)

stats = pd.DataFrame(index=vmPFC.index)
stats['vmPFC-hc_tail'] = vmPFC['hc_tail_ret_uni']
stats['vmPFC-hc_body'] = vmPFC['hc_body_ret_uni']
stats['vmPFC-hc_head'] = vmPFC['hc_head_ret_uni']
stats['vmPFC-amyg_bla'] = vmPFC['amyg_bla_ret_uni']
stats['vmPFC-amyg_cem'] = vmPFC['amyg_cem_ret_uni']

stats['dACC-hc_tail'] = vmPFC['hc_tail_ret_uni']
stats['dACC-hc_body'] = vmPFC['hc_body_ret_uni']
stats['dACC-hc_head'] = vmPFC['hc_head_ret_uni']
stats['dACC-amyg_bla'] = vmPFC['amyg_bla_ret_uni']
stats['dACC-amyg_cem'] = vmPFC['amyg_cem_ret_uni']

stats = stats.reset_index()
stats['group'] = stats.subject.apply(lgroup)

stats.to_csv('betaseries_lmm.csv',index=False)
