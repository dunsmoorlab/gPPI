from fg_config import *
from subprocess import call
from os.path import join

dest = '/mnt/c/Users/ACH/Desktop/fc-bids'
dest_preproc = join(dest,'preproc')
dest_model = join(dest,'derivatives','model')

for folder in [dest,dest_model]: mkdir(folder)


for sub in all_sub_args:
    print(sub)
    subj = bids_meta(sub)
    # sub_func = join(dest_preproc,subj.fsub,'func');mkdir(sub_func)
    sub_model = join(dest_model,subj.fsub);mkdir(sub_model)
    for run in [1,2]:
        if sub == 107 and run == 2:
            pass
        else:            
            # epi = f'{subj.fsub}_ses-2_task-localizer_run-0{run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'
            # os.system(f'cp -r {join(subj.preproc_dir,"func",epi)} {join(sub_func,epi)}')

            run_model = join(dest_model,subj.fsub,f'localizer_run-0{run}')
            os.system(f'cp -r {join(subj.model_dir,f"localizer_run-0{run}")} {run_model}')
