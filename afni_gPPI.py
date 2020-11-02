from fg_config import *
seeds = ['hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']

def reg_copes(sub):
    subj = bids_meta(sub)

    reg_me = ['cope1','cope2','cope3','cope4','cope5','cope6','cope7','cope8','cope9','cope10','cope11','cope12','cope13','cope14','res4d']

    for run in [1,2,3]:
        _run = f'memory_run-0{run}'
        reg_dir = f'{subj.model_dir}/{_run}/{subj.fsub}_{_run}_gPPI.feat/reg'

        for roi in seeds:
            stats = f'{subj.model_dir}/{_run}/{roi}/source.feat/stats'
            reg_std = f'{subj.model_dir}/{_run}/{roi}/source.feat/reg_std/stats'
            mkdir(reg_std)

            #need to register each cope
            for stat in reg_me:

                reg_cmd = f'flirt -ref {reg_dir}/standard \
                                  -in {stats}/{stat}.nii.gz \
                                  -out {reg_std}/{stat}.nii.gz \
                                  -applyxfm \
                                  -init {reg_dir}/example_func2standard.mat \
                                  -interp trilinear \
                                  -datatype float'
                os.system(reg_cmd)

def smooth_est(sub):
    subj = bids_meta(sub)
    #for each seed region in each run we need to get a spatial noise estimate
    for run in [1,2,3]:
        for seed in seeds:
            reg_std = f'{subj.model_dir}/memory_run-0{run}/{seed}/source.feat/reg_std/stats'
            os.chdir(reg_std)
            est_noise = f'3dFWHMx -mask {std_2009_brain_mask} \
                                  -input {reg_std}/res4d.nii.gz \
                                  -acf >> noise_estimates.txt'
#started at 9:43, need 15 repititions per sub
#
def gPPI_datatables():
    cope_map = {'acquisition':{'CS+':4,
                              'CS-':5},
                 'extinction':{'CS+':7,
                              'CS-':8},
                              }

    out = './gPPI_MVM'
    for seed in seeds:
        # Subj Sess Encode Condition Response InputFile
        df = pd.DataFrame({'InputFile':'','Group':''},
                index=pd.MultiIndex.from_product([all_sub_args,
                                                 ['acquisition','extinction'],
                                                 ['CSp','CSm'],
                                                 [1,2,3]],
                                                 names=['Subj','Encode','Condition','Sess']))
        for sub in all_sub_args:
            if sub < 100:
                df.loc[sub,'Group'] = 'healthy'
            else:
                df.loc[sub,'Group'] = 'ptss'
            subj = bids_meta(sub)
            for phase in ['acquisition','extinction']:
                for c, con in enumerate(cons):
                    for sess in [1,2,3]:
                        df.loc[(sub,phase,consp[c],sess),'InputFile'] = f'{subj.model_dir}/memory_run-0{sess}/{seed}/source.feat/reg_std/stats/cope{cope_map[phase][con]}.nii.gz'
        df = df.reset_index()[['Subj','Group','Encode','Condition','Sess','InputFile']]

        assert df.InputFile.apply(lambda x: os.path.exists(x)).sum() == 576

        df.to_csv(f'{out}/{seed}_dataTable.txt',index=False,sep=' ')

# def collect_smooth():
#     dfs = {}
#     for seed in seeds:
#         dfs[seed] = {}
#         for sub in all_sub_args:
#             subj = bids_meta(sub)
#             dfs[seed][sub] = {}
#             for run in [1,2,3]:
#                 _df[sub][run] = pd.read_csv(f'{subj.model_dir}/memory_run-0{run}/{seed}/source.feat/noise_estimates.txt',sep=' '
#                                 ).loc[0].dropna().values
#         df = pd.DataFrame.from_dict(_df).unstack().reset_index()
#         est = df[0].values.mean(axis=0)
#         est.tofile('sm_events/basic_model_mean_smooth.txt',sep=' ', format='%s')


def run_wrap():
    #for smooth or reg
    for sub in all_sub_args:
        os.system(f"echo singularity run --cleanenv $SCRATCH/bids-apps/neurosft.simg python $HOME/gPPI/wrap_glm_utils.py -s {sub} >> jobs/smooth2_gPPI_job.txt")
    os.system('launch -N 48 -n 48 -J smooth2 -s jobs/smooth2_gPPI_job.txt -m achennings@utexas.edu -p normal -r 10:00:00 -A LewPea_MRI_Analysis')

    #for MVM
    # for seed in ['amyg_cem','amyg_bla','hc_head']:
    # for seed in ['hc_body','hc_tail']:
        # os.system(f'echo "bash -x /home1/05426/ach3377/gPPI/gPPI_MVM/{seed}_MVM.txt" >> gPPI_MVM/{seed}_job.txt')
        # os.system(f'launch -N 1 -n 64 -p largemem512GB -J {seed} -s gPPI_MVM/{seed}_job.txt -m achennings@utexas.edu -r 48:00:00 -A LewPea_MRI_Analysis')

