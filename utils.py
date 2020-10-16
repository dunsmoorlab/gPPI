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