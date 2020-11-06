from fg_config import *
import nibabel as nib
from nilearn.image import mean_img
seeds = ['hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']
cope_map = {'acquisition':{'CSp':4,
                         'CSm':5},
           'extinction':{'CSp':7,
                         'CSm':8}}
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

def mean_copes(sub):
    print(sub)
    subj = bids_meta(sub)
    gPPI_dir = f'{subj.model_dir}/gPPI_results';mkdir(gPPI_dir)

    for seed in seeds:
        seed_dir = f'{gPPI_dir}/{seed}';mkdir(seed_dir)

        for phase in ['acquisition','extinction']:
            for con in consp:
                cope = mean_img([nib.load(f'{subj.model_dir}/memory_run-0{run}/{seed}/source.feat/reg_std/stats/cope{cope_map[phase][con]}.nii.gz') for run in [1,2,3]])
                nib.save(cope,f'{seed_dir}/{phase}_{con}.nii.gz')

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
            os.system(est_noise)

def gPPI_datatables():
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

def collect_smooth():
    dfs = {}
    for seed in seeds:
        dfs[seed] = {}
        for sub in all_sub_args:
            subj = bids_meta(sub)
            dfs[seed][sub] = {}
            for run in [1,2,3]:
                try:
                    dfs[seed][sub][run] = pd.read_csv(f'{subj.model_dir}/memory_run-0{run}/{seed}/source.feat/reg_std/stats/noise_estimates.txt',sep=' '
                                   ).loc[0].dropna().values
                except:
                    print('whoopsie')
        dfs[seed] = pd.DataFrame.from_dict(dfs[seed]).unstack().reset_index()
        est = dfs[seed][0].values.mean(axis=0)
        est.tofile(f'{SCRATCH}/gPPI_MVM/{seed}_mean_smooth.txt',sep=' ', format='%s')

def paired_ttest(subs=None,name='',seed=None):
    out_parent = '/scratch/05426/ach3377/gPPI_comps'
    out_dir = f'{out_parent}/{name}'
    mkdir(out_dir)
    
    for seed in [seed]:
        seed_dir = f'{out_dir}/{seed}';mkdir(seed_dir)
        setA = ''
        setB = ''
        for s, sub in enumerate(subs):
            subj = bids_meta(sub)
            setA += f'{subj.model_dir}/gPPI_results/{seed}/extinction_CSp.nii.gz '
            setB += f'{subj.model_dir}/gPPI_results/{seed}/acquisition_CSp.nii.gz '

        n_cors = 'export OMP_NUM_THREADS=48'
        cd_cmd = f'cd {seed_dir}'
        clustsim_cmd = f'3dttest++ -setA {setA} \
                                   -setB {setB} \
                                   -AminusB \
                                   -paired \
                                   -Clustsim 48 \
                                   -mask {gm_1mm_thr} \
                                   -prefix {name}_clst-ttest'
        
        script = f'{seed_dir}/ttest_script.txt'
        os.system(f'rm {script}')
    
        for cmd in [n_cors, cd_cmd, clustsim_cmd]:
            os.system(f"echo {cmd} >> {script}")
    
        jobfile = f'/home1/05426/ach3377/gPPI/jobs/{name}_{seed}_gPPI_comp_job.txt'
        os.system(f'rm {jobfile}')

        os.system(f'singularity run --cleanenv \
                        /scratch/05426/ach3377/bids-apps/neurosft.simg \
                        bash -x {script}')
        #not run here, just submiting a job
        # os.system(f'echo singularity run --cleanenv \
        #                 /scratch/05426/ach3377/bids-apps/neurosft.simg \
        #                 bash -x {script} >> {jobfile}')

        # os.system(f'launch -N 1 \
        #                    -n 1 \
        #                    -J 3dttest++ \
        #                    -s {jobfile} \
        #                    -m achennings@utexas.edu \
        #                    -p normal \
        #                    -r 1:00:00 \
        #                    -A LewPea_MRI_Analysis')
paired_ttest(subs=sub_args,name='healthy_CSpE__CSpA',seed='hc_head')
paired_ttest(subs=p_sub_args,name='ptsd_CSpE__CSpA')


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

