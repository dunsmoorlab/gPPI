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
    import time
    seeds = ['hc_tail','hc_body','hc_head','amyg_bla','amyg_cem']
    subj = bids_meta(sub)
    #for each seed region in each run we need to get a spatial noise estimate
    for run in [1,2,3]:
        for seed in seeds:
            reg_std = f'{subj.model_dir}/memory_run-0{run}/{seed}/source.feat/reg_std/stats'
            os.chdir(reg_std)
            time.sleep(1)
            est_noise = f'3dFWHMx -mask {std_2009_brain_mask} \
                                  -input {reg_std}/res4d.nii.gz \
                                  -acf >> noise_estimates.txt'
            os.system(est_noise)
            time.sleep(5)
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
        print(seed)
        dfs[seed] = {}
        for sub in all_sub_args:
            print(sub)
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

def clustsim(seed=None):
    root_dir = '/scratch/05426/ach3377/gPPI_MVM'
    est = np.loadtxt(f'{root_dir}/{seed}_mean_smooth.txt')
    cmd1 = 'export OMP_NUM_THREADS=48'
    cmd2 = f'3dClustSim -OKsmallmask \
                       -mask $SCRATCH/standard/gm_1mm_thr.nii.gz \
                       -acf {est[0]} {est[1]} {est[2]} \
                       > {root_dir}/{seed}_clustsim_output.txt' 
    
    script = f'{root_dir}/{seed}_clustsim_script.txt'
    os.system(f'rm {script}')

    for cmd in [cmd1,cmd2]:
        os.system(f"echo {cmd} >> {script}")

    jobfile = f'/home1/05426/ach3377/gPPI/jobs/{seed}_clustsim_job.txt'
    os.system(f'rm {jobfile}')

    os.system(f'echo singularity run --cleanenv \
                    /scratch/05426/ach3377/bids-apps/neurosft.simg \
                    bash -x {script} >> {jobfile}')
    
    os.system(f'launch -N 1 \
                       -n 1 \
                       -J clustsim \
                       -s {jobfile} \
                       -m achennings@utexas.edu \
                       -p normal \
                       -r 1:00:00 \
                       -A LewPea_MRI_Analysis')

def clusterize(seed,thr1,thr2):
    from nibabel.brikhead import parse_AFNI_header
    here = os.getcwd()
    root = '/Volumes/DunsmoorRed/gPPI_MVM'
    os.chdir('/Volumes/DunsmoorRed/gPPI_MVM')
    header = parse_AFNI_header(f'{seed}+tlrc.HEAD')
    out_dir = f'{root}/{seed}_clusterize';mkdir(out_dir)


    names = header['BRICK_LABS'].replace(':','_'
                                ).replace('  F','_F'
                                ).replace(' F','_F'
                                ).replace(' t', '_t'
                                ).replace(' Z','_Z'
                                ).split('~')

    for i, name in enumerate(names):
        if 'Intercept' in name:
            pass
        elif '_F' in name or '_t' in name or '_Z' in name:
            cmap = f'{out_dir}/{name}_cluster_map.nii.gz';os.system(f'rm {cmap}')
            ctxt = f'{out_dir}/{name}_cluster.txt';os.system(f'rm {ctxt}')
            where = f'{out_dir}/{name}_where.txt';os.system(f'rm {where}')

            if '_F' in name:
                side = '1sided RIGHT_TAIL'
                thr = thr2
            elif '_t' in name or '_Z' in name:
                side = '2sided'
                thr = thr1


            cmd = f"3dClusterize \
                      -inset {seed}+tlrc \
                      -ithr {i} \
                      -mask /Users/ach3377/Desktop/standard/gm_1mm_thr.nii.gz \
                      -{side} p=0.001 \
                      -clust_nvox {thr} \
                      -NN 3 \
                      -pref_map {cmap} \
                      > {ctxt}"
            os.system(cmd)
            

            if os.path.exists(cmap):
                w_cmd = f"whereami -coord_file {ctxt}'[1,2,3]' > {where}"
                os.system(w_cmd)
    os.chdir(here)
# clusterize('hc_tail',75.2,90.5)
# clusterize('hc_body',75.3,89.2)
clusterize('hc_head',76.2,91.1)
clusterize('amyg_bla',77.3,91.2)
clusterize('amyg_cem',75.2,89.8)

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
# paired_ttest(subs=sub_args,name='healthy_CSpE__CSpA',seed='hc_head')
# paired_ttest(subs=p_sub_args,name='ptsd_CSpE__CSpA')


def run_wrap():
    #for smooth or reg
    for sub in all_sub_args:
        os.system(f"echo singularity run --cleanenv $SCRATCH/bids-apps/neurosft.simg python $HOME/gPPI/wrap_glm_utils.py -s {sub} >> jobs/smooth3_gPPI_job.txt")
    os.system('launch -N 48 -n 48 -J smooth3 -s jobs/smooth3_gPPI_job.txt -m achennings@utexas.edu -p normal -r 10:00:00 -A LewPea_MRI_Analysis -d 3152834')

    #for MVM
    # for seed in ['amyg_cem','amyg_bla','hc_head']:
    # for seed in ['hc_body','hc_tail']:
        # os.system(f'echo "bash -x /home1/05426/ach3377/gPPI/gPPI_MVM/{seed}_MVM.txt" >> gPPI_MVM/{seed}_job.txt')
        # os.system(f'launch -N 1 -n 64 -p largemem512GB -J {seed} -s gPPI_MVM/{seed}_job.txt -m achennings@utexas.edu -r 48:00:00 -A LewPea_MRI_Analysis')

