from fg_config import *

def collect_sub_files():

    for sub in all_sub_args:

        subj = bids_meta(sub)

        os.system(f'cp {subj.feat_dir}/localizer_lvl2.gfeat/cope6.feat/stats/cope1.nii.gz {SCRATCH}object_cortex/object_scram/{subj.fsub}.nii.gz')
        os.system(f'cp {subj.feat_dir}/localizer_lvl2.gfeat/cope7.feat/stats/cope1.nii.gz {SCRATCH}object_cortex/a_over_t/{subj.fsub}.nii.gz')
        os.system(f'cp {subj.feat_dir}/localizer_lvl2.gfeat/cope8.feat/stats/cope1.nii.gz {SCRATCH}object_cortex/t_over_a/{subj.fsub}.nii.gz')




def some_stats():

    out_parent = '/scratch/05426/ach3377/object_cortex'
    out_dir = f'{out_parent}/at_ttest'
    mkdir(out_dir)
    
    setA = ''
    for s, sub in enumerate(all_sub_args):
        setA += f'{out_parent}/a_over_t/{bids_meta(sub).fsub}.nii.gz '
    
    n_cors = 'export OMP_NUM_THREADS=24'
    cd_cmd = f'cd {out_dir}'
    clustsim_cmd = f'3dttest++ -setA {setA} \
                               -Clustsim 24 \
                               -mask /scratch/05426/ach3377/object_cortex/object_scram_ttest/object_scram_ClusterMask.nii.gz \
                               -prefix clst-ttest'

    script = f'{out_dir}/ttest_script.txt'
    for cmd in [n_cors, cd_cmd, clustsim_cmd]:
        os.system(f"echo {cmd} >> {script}")
    
    jobfile = f'/home1/05426/ach3377/gPPI/jobs/object_cortex_job2.txt'
    os.system(f'rm {jobfile}')

    #not run here, just submiting a job
    os.system(f'echo singularity run --cleanenv \
                    /scratch/05426/ach3377/bids-apps/neurosft.simg \
                    bash -x {script} >> {jobfile}')

    os.system(f'launch -N 1 \
                       -n 1 \
                       -J 3dttest++ \
                       -s {jobfile} \
                       -m achennings@utexas.edu \
                       -p normal \
                       -r 10:00:00 \
                       -A LewPea_MRI_Analysis')

def cluster_those_stats(folder='/scratch/05426/ach3377/object_cortex/at_ttest',thr=.001,nvox=195,mask='/scratch/05426/ach3377/standard/gm_1mm_thr.nii.gz',tail='two-sided'):

        # here = os.getcwd()
        name = (folder.split('/')[-1]).split('_ttest')[0]
        cd_cmd = f'cd {folder}'
        
        if tail == 'one-sided':
            side = '1sided RIGHT_TAIL'
        elif tail == 'two-sided':
            side = '2sided'

        cmap = f'{name}_ClusterMap.nii.gz';os.system(f'rm {cmap}')
        ceff = f'{name}_ClusterEffEst.nii.gz';os.system(f'rm {ceff}')

        ctxt = f'{name}_cluster.txt'#;os.system(f'rm {ctxt}')
        where = f'{name}_where.txt'#;os.system(f'rm {where}')
        
        cmd = f"3dClusterize -inset clst-ttest+tlrc \
                   -ithr 1 \
                   -idat 0 \
                   -mask {mask} \
                   -NN 3 \
                   -{side} p={thr} \
                   -clust_nvox {nvox} \
                   -pref_map {cmap} \
                   -pref_dat {ceff} >> {ctxt}"

        
        # os.system(cmd)
        
        # if os.path.exists(f'{name}_ClusterMap.nii.gz'):
        w_cmd = f"whereami -coord_file {ctxt}'[13,14,15]' > {where}"
            # os.system(w_cmd)
        
        script = f'{folder}/cluster_script.txt'
        for cmd in [cd_cmd, cmd, w_cmd]:
            os.system(f"echo {cmd} >> {script}")
        
        # os.chdir(here)
    #singularity run --cleanenv /scratch/05426/ach3377/bids-apps/neurosft.simg bash -x /scratch/05426/ach3377/object_cortex/at_ttest/cluster_script.txt