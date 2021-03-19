from fg_config import *

for sub in all_sub_args:
    subj = bids_meta(sub)
    for run in [1,2,3]:
        df = pd.read_csv(f'{subj.subj_dir}/ses-2/func/{subj.fsub}_ses-2_task-memory_run-0{run}_events.tsv',sep='\t')
        df.trial_type = df.trial_type + '_' + df.encode_phase
        df.to_csv(f'{subj.preproc_dir}/func/{subj.fsub}_ses-2_task-memory_run-0{run}_space-MNI152NLin2009cAsym_desc-preproc_denoised_bold_events.tsv',sep='\t',index=False)

from fg_config import *
from scipy.io import loadmat

rois = ['vmPFC','dACC','amyg_cem','amyg_bla','hc_head','hc_body','hc_tail','animal','tool']
conds = ['CS+_acquisition', 'CS+_baseline', 'CS+_extinction', 'CS+_foil', 'CS-_acquisition', 'CS-_baseline', 'CS-_extinction', 'CS-_foil']
df = pd.DataFrame({'conn':0.0},index=pd.MultiIndex.from_product([conds,all_sub_args,rois,rois],names=['condition','subject','seed','target']))

for c, cond in enumerate(conds):
    stats = loadmat(f'conn_project01/results/firstlevel/gPPI_03/resultsROI_Condition00{c+1}.mat')
    for s, sub in enumerate(all_sub_args):
        for si, seed in enumerate(rois):
            for ti, target in enumerate(rois):
                df.loc[(cond,sub,seed,target),'conn'] = stats['Z'][si,ti,s]

df = df.dropna(subset=['conn'])

df['csp_cat'] = ''
df = df.reset_index()
df['group'] = df.subject.apply(lgroup)
df['phase'] = df.condition.apply(lambda x: x[4:])
df['condition'] = df.condition.apply(lambda x: x[:3])
df = df.set_index(['group','phase','condition','seed','target','subject']).sort_index()

df = df.reset_index().set_index('subject')
for sub in all_sub_args:
    subj = bids_meta(sub)
    csp = subj.mem_df[subj.mem_df.trial_type == 'CS+']['stimulus'].values[0].split('/')[0][:-1]
    df.loc[sub,'csp_cat'] = csp
#     csm = subj.mem_df[subj.mem_df.trial_type == 'CS-']['stimulus'].values[0].split('/')[0][:-1]

#     cats[sub] = {csp:'csp',
#                 csm:'csm'}
#     for roi in rois:
#         if roi not in ['animal','tool']:
#             cats[sub][roi] = roi

# df['seed'] = df[['subject','seed']].apply(lambda x: cats[x[0]][x[1]],axis=1)
# df['target'] = df[['subject','target']].apply(lambda x: cats[x[0]][x[1]],axis=1)

df.to_csv('conn_stats_lmm.csv')


pfc = df.reset_index().set_index(['target','group','phase','condition','seed','subject']).sort_index()
pfc = pfc.loc['vmPFC'] - pfc.loc['dACC']
pfc = pfc.dropna(subset=['conn'])
pfc.to_csv('conn_pfc_diff_lmm.csv')


aHPC = df.reset_index()
aHPC = aHPC[aHPC.seed == 'hc_head'][aHPC.target.isin(['vmPFC','dACC'])]