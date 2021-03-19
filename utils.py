#copy the confounds again because they got deleted off of scratch
from fg_config import *
prep_dir = '/Volumes/DunsmoorRed/fc-bids/derivatives/fmriprep'
tmp = '/Volumes/DunsmoorRed/fc-bids/derivatives/tmp'
mkdir(tmp)

for sub in all_sub_args:
    subj = bids_meta(sub)
    subdir = os.path.join(tmp,subj.fsub); mkdir(subdir)
    ses1 = os.path.join(subdir,'ses-1','func'); mkdir(ses1)
    ses2 = os.path.join(subdir,'ses-2','func'); mkdir(ses2)

    for task in tasks:
        _ses = tasks[task]['ses']
        confounds = f'{prep_dir}/{subj.fsub}/ses-{_ses}/func/{subj.fsub}_ses-{_ses}_task-{task}_desc-confounds_regressors.tsv'
        out = f'{subdir}/ses-{_ses}/func/{subj.fsub}_ses-{_ses}_task-{task}_desc-confounds_regressors.tsv'

        os.system(f'cp {confounds} {out}')

'''putting memory into the big data frame with all the data'''
df = pd.read_csv('all_data_lmm.csv').set_index(['subject','stimulus']).sort_index()
df['response'] = 0
df['high_confidence_accuracy'] = ''
df['shock'] = 'CS'

for sub in all_sub_args:
    subj = bids_meta(sub)

    subj.mem_df['shock'] = 'CS'
    subj.mem_df = subj.mem_df.set_index(['encode_phase','stimulus']).sort_index()
    subj.behav['acquisition'] = subj.behav['acquisition'].set_index('stimulus').sort_index()
    subj.mem_df.loc['acquisition','shock'] = subj.behav['acquisition']['shock'].values

    subdf = subj.mem_df.copy().dropna(subset=['response']).reset_index(level='encode_phase').sort_index()

    df.loc[sub,'response'] = subdf.response.values
    df.loc[sub,'high_confidence_accuracy'] = subdf.high_confidence_accuracy.values
    df.loc[sub,'shock'] = subdf.shock.values

df = df.rename(columns={'high_confidence_accuracy':'mem_acc'})
df.to_csv('all_data_lmm.csv')

df = df.reset_index()
'''pfc df'''
pfc = df[df.phase != 'foil'].melt(id_vars=['subject','stimulus','group','phase','condition','response','mem_acc','shock'],value_vars=['dACC_ers','vmPFC_ers'],var_name='roi',value_name='ers')
pfc.roi = pfc.roi.apply(lambda x: x[:-4])
pfc = pfc.set_index(['subject','stimulus','group','phase','condition','response','mem_acc','shock','roi']).sort_index()

pfc_uni = df[df.phase != 'foil'].melt(id_vars=['subject','stimulus','group','phase','condition','response','mem_acc','shock'],value_vars=['dACC_ret_uni','vmPFC_ret_uni'],var_name='roi',value_name='uni')
pfc_uni.roi = pfc_uni.roi.apply(lambda x: x[:-8])
pfc_uni = pfc_uni.set_index(['subject','stimulus','group','phase','condition','response','mem_acc','shock','roi']).sort_index()

pfc = pd.concat((pfc,pfc_uni),axis=1)
pfc.to_csv('pfc_ers_cleaned_lmm.csv')

'''subcortical df'''
hca = df[df.phase != 'foil'].melt(id_vars=['subject','stimulus','group','phase','condition','response','mem_acc','shock'],value_vars=['hc_tail_ers','hc_body_ers','hc_head_ers','amyg_bla_ers','amyg_cem_ers'],var_name='roi',value_name='ers')
hca.roi = hca.roi.apply(lambda x: x[:-4])
hca = hca.set_index(['subject','stimulus','group','phase','condition','response','mem_acc','shock','roi']).sort_index()

hca_uni = df[df.phase != 'foil'].melt(id_vars=['subject','stimulus','group','phase','condition','response','mem_acc','shock'],value_vars=['hc_tail_ret_uni','hc_body_ret_uni','hc_head_ret_uni','amyg_bla_ret_uni','amyg_cem_ret_uni'],var_name='roi',value_name='uni')
hca_uni.roi = hca_uni.roi.apply(lambda x: x[:-8])
hca_uni = hca_uni.set_index(['subject','stimulus','group','phase','condition','response','mem_acc','shock','roi']).sort_index()

hca = pd.concat((hca,hca_uni),axis=1)
hca.to_csv('subcort_ers_cleaned_lmm.csv')

