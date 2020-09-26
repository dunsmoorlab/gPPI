from fg_config import *

def apply_mask(mask=None,target=None):

    coor = np.where(mask == 1)
    values = target[coor]
    if values.ndim > 1:
        values = np.transpose(values) #swap axes to get feature X sample
    return values


'''GLM UTILS'''
# encoding_phases = ['baseline','acquisition','extinction']
# memory_phases = ['memory_run-01','memory_run-02','memory_run-03']
# consp = ['CSp','CSm']

pe_map = {'baseline':{'CSp':{'baseline':1,
                          'acquisition':2,
                          'extinction':3},
                   'CSm':{'baseline':4,
                          'acquisition':5,
                          'extinction':6}
                          },
        'acquisition':{'CSp':{'baseline':7,
                          'acquisition':8,
                          'extinction':9},
                   'CSm':{'baseline':10,
                          'acquisition':11,
                          'extinction':12}
                          },
        'extinction':{'CSp':{'baseline':13,
                          'acquisition':14,
                          'extinction':15},
                   'CSm':{'baseline':16,
                          'acquisition':17,
                          'extinction':18}
                          }
        }
def sm_glm_events(sub):
    data_table = pd.DataFrame(columns=['Subj','Sess','Encode','Condition','Response','InputFile'])
    for sub in xcl_sub_args:
        subj = bids_meta(sub)
        events = pd.read_csv(f'sm_events/{subj.fsub}/sm_events.csv').set_index('phase').sort_index()
        events.trial_type = events.trial_type.apply(lambda x: 'CSp' if x == 'CS+' else 'CSm')
        events['PM'] = 1.0

        for mem_phase in memory_phases:

            mem_phase_dir = os.path.join('sm_events',subj.fsub,mem_phase)
            mkdir(mem_phase_dir)

            mem_phase_dat = events.loc[mem_phase].copy(
                                ).set_index(['encode_phase','trial_type']
                                ).sort_values(by='onset')

            for con in consp:

                for encode_phase in encoding_phases:
                    
                    dat = mem_phase_dat.loc[(encode_phase,con)].copy()

                    for sm_resp in encoding_phases:

                        smdat = dat[dat.source_memory == sm_resp].copy()

                        if smdat.shape[0] == 0: 
                            # out_events = pd.DataFrame({'onset':0.0,'duration':0.0,'PM':0.0},index=[0])
                            # os.system(f'echo {sub}\t{mem_phase}\t{encode_phase}\t{con}\t{sm_resp} >> sm_events/missing_evs.txt')
                            pass
                        else:                   
                            # out_events = smdat[['onset','duration','PM']].copy()
                            data_table = data_table.append({'Subj':sub,
                                                            'Sess':mem_phase,
                                                            'Encode':encode_phase,
                                                            'Condition':con,
                                                            'Response':sm_resp,
                                                            'InputFile':f'/scratch/05426/ach3377/fc-bids/derivatives/model/{subj.fsub}/{mem_phase}/source_memory.feat/reg_standard/stats/cope{pe_map[encode_phase][con][sm_resp]}.nii.gz'},
                                                            ignore_index=True)
                        # out_events.to_csv( os.path.join(mem_phase_dir, f'{encode_phase}_{con}_{sm_resp}.txt'),
                        #     sep='\t', float_format='%.8e', index=False, header=False)

                #and the foils
                mem_phase_dat.loc[('foil',con),['onset','duration','PM']].to_csv( os.path.join(mem_phase_dir, f'foil_{con}.txt'),
                            sep='\t', float_format='%.8e', index=False, header=False)
#data_table.to_csv('sm_events/afni_dataTable.txt',index=False,sep='\t')
def copy_sm_events(sub):
    subj = bids_meta(sub)
    for run in [1,2,3]:
        in_folder = f'sm_events/{subj.fsub}/memory_run-0{run}'
        out_folder = os.path.join(subj.model_dir,f'memory_run-0{run}','sm_events')
        os.system(f'cp -r {in_folder} {out_folder}')

    ref = pd.DataFrame({'memory_phase':  np.tile(memory_phases,18),
                        'encoding_phase':np.repeat(encoding_phases,18),
                        'condition':     np.repeat(np.tile(consp,3),9),
                        'source_memory' :np.tile(np.repeat(encoding_phases,3),6)}
                        ).reset_index().rename(columns={'index':'input'})
    ref.input += 1
    ref = ref.set_index(['memory_phase','encoding_phase','condition','source_memory'])

    encode_num = {'baseline':1,'acquisition':7,'extinction':13}
    con_num = {'CSp':0,'CSm':3}
    resp_num = {'baseline':0,'acquisition':1,'extinction':2}
    missing = pd.read_csv('sm_events/missing_evs.txt',sep='\t',header=None
                        ).rename(columns={0:'subject',1:'memory_phase',2:'encoding_phase',3:'condition',4:'source_memory'})
    missing.source_memory = missing.source_memory.apply(lambda x: x[:-1])
    missing['cope'] = 0
    for i in missing.index:missing.loc[i,'cope'] = encode_num[missing.loc[i,'encoding_phase']] + con_num[missing.loc[i,'condition']] + resp_num[missing.loc[i,'source_memory']]


def zero_missing_ev(sub):
    template = 'feats/template_source_memory_lvl2.fsf'

    subj = bids_meta(sub)
    lvl2 = os.path.join(subj.feat_dir,f'{subj.fsub}_source_memory_lvl2.fsf')
    tmp_out = os.path.join(subj.feat_dir,'temp_source_memory_lvl2.fsf')

    sub_missing = missing[missing.subject == sub].copy()
    replacements = {}
    for b in sub_missing.index:
        _input = ref.loc[(sub_missing.loc[b,'memory_phase'],
                        sub_missing.loc[b,'encoding_phase'],
                        sub_missing.loc[b,'condition'],
                        sub_missing.loc[b,'source_memory']),
                        'input']
        _cope = sub_missing.loc[b,'cope']

        search_for = f'set fmri(evg{_input}.{_cope}) 1.0'
        replace_with = f'set fmri(evg{_input}.{_cope}) 0'
        
        grp_search_for = f'set fmri(groupmem.{_input}) 1'
        grp_replace_with = f'set fmri(groupmem.{_input}) 2'
        
        replacements[search_for] = replace_with
        replacements[grp_search_for] = grp_replace_with
        # with open(os.path.join(gPPI_codebase,'feats','%s.fsf'%(template))) as infile: 

    with open(lvl2) as infile:
        with open(tmp_out, 'w') as outfile:
            for line in infile:
                for src, target in replacements.items():
                    line = line.replace(src, target)
                outfile.write(line)

    os.system(f'mv {tmp_out} {lvl2}')

    for cope in range(1,21):
        if sub_missing[sub_missing.cope == cope].shape[0] == 3:
            os.system(f'echo "{sub} {cope}" >> sm_events/missing_evs_group_level.txt')


def group_level_missing_inputs():

    cope_dep = {21:[1,2],
                22:[1,2,3],
                23:[2,5],
                24:[7,8,9],
                25:[14,15],
                26:[13,14,15],
                27:[14,17],
                28:[2,8,14],
                29:[1,8,15],
                30:[1,2,14,15]}

    template = 'feats/template_source_memory_lvl3.fsf'
    out_dir = 'sm_events/group_fsfs'
    out_inputs = 'sm_events/group_inputs'

    missing = pd.read_csv('sm_events/missing_evs_group_level.txt',sep=' ',header=None
                    ).rename(columns={0:'subject',1:'cope'})
    missing['input'] = missing.subject.apply(lambda x: xcl_sub_args.index(x))

    for cope in range(1,31):
        #need this for every cope
        replacements = {'COPEID':f'cope{cope}'}

        if cope <= 20:#technically we don't need it for the foils, but this is the range of single copes
        
            # sub_missing = missing[missing.cope == cope].copy().input.unique()
            sub_missing = missing[missing.cope == cope].copy().subject.unique()

        else:#otherwise we use the entered dependencies to see who needs to Zero'd
            sub_missing = pd.concat([missing[missing.cope == c].copy() for c in cope_dep[cope]]).subject.unique()

        nsub = 34
        for sub in xcl_sub_args:
            subj = bids_meta(sub)
            if sub in sub_missing:
                nsub = nsub -1
            else:
                os.system(f'echo /scratch/05426/ach3377/fc-bids/derivatives/model/{subj.fsub}/all_memory_runs/source_memory.feat/stats/cope{cope}.nii.gz >> sm_events/group_inputs/cope{cope}_inputs.txt')
                pass
        if nsub != 34: print(f'cope{cope}\t{nsub}')

def group_level_autofill_fsf():

    for cope in range(1,31):

        inputs = pd.read_csv(f'sm_events/group_inputs/cope{cope}_inputs.txt',header=None)
        inputs[0] = inputs[0].apply(lambda x: x[:-1])
        template = f'sm_events/group_fsfs/template_cope{cope}.fsf'
        out_feat = f'sm_events/group_fsfs/cope{cope}.fsf'

        #get all the replacements ready        
        replacements = {'source_memory_group_glm/COPEID':f'source_memory_group_glm/cope{cope}'}

        for _j, _input in enumerate(inputs[0]):

            sub_num = '{:03d}'.format(xcl_sub_args[_j])
            search_for   = f'set feat_files({int(_j + 1)}) "/scratch/05426/ach3377/fc-bids/derivatives/model/sub-FC{sub_num}/feats/source_memory_lvl2.gfeat/cope1.feat/stats/COPEID.nii.gz"'
            replace_with = f'set feat_files({int(_j + 1)}) "{_input}"'

            replacements[search_for] = replace_with

        #fill out template
        with open(template) as infile:
            with open(out_feat, 'w') as outfile:
                for line in infile:
                    for src, target in replacements.items():
                        line = line.replace(src, target)
                    outfile.write(line)

        os.system(f'echo feat {gPPI_codebase}{out_feat} >> jobs/source_memory_lvl3_job.txt')

# from nilearn.image import concat_imgs
# from os.path import join
# import nibabel as nib

def cat_mem_runs(sub):
# for sub in xcl_sub_args:
    subj = bids_meta(sub)
    out_dir = join(subj.model_dir,'all_memory_runs');mkdir(out_dir)
    sm_out = out_dir+'/sm_events';mkdir(sm_out) 

    #concatenate the epi data
    mem_data = [nib.load(join(subj.func,f'{subj.fsub}_ses-2_task-memory_run-0{run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz')) for run in [1,2,3]]

    mem_data = concat_imgs(mem_data)
    nib.save(mem_data,join(out_dir,'all_memory_runs.nii.gz'))

    #concatenate the confounds
    confounds = [pd.read_csv(join(subj.model_dir,f'memory_run-0{run}','confounds.txt'),sep='\t',header=None) for run in [1,2,3]] 
    confounds = pd.concat(confounds).reset_index(drop=True)
    
    #run onsets
    confounds = pd.concat((confounds, pd.DataFrame(np.kron(np.eye(3), np.ones((310,1)))) ), axis=1)
    
    confounds.to_csv(join(out_dir,'confounds.txt'),
              sep='\t',float_format='%.8e', index=False, header=False)
    
    #concatenate the regressors
    for con in consp:
        for encode in encoding_phases:
            for response in encoding_phases:
                events = []
                for i, run in enumerate([1,2,3]):
                    
                    evs = pd.read_csv(f'{subj.model_dir}/memory_run-0{run}/sm_events/{encode}_{con}_{response}.txt', sep='\t', header=None)
                    if evs.sum(axis=1)[0] == 0: 
                        pass
                    else:
                        #correct for start time
                        evs[0] += i*620
                        events.append(evs)

                #if regressor is still empty
                if not events:
                    out_events = pd.DataFrame({'onset':0.0,'duration':0.0,'PM':0.0},index=[0])
                    os.system(f'echo {sub}\t{encode}\t{con}\t{response} >> sm_events/missing_evs_concat.txt')

                else:
                    out_events = pd.concat(events).reset_index(drop=True)

                out_events.to_csv( f'{sm_out}/{encode}_{con}_{response}.txt',
                        sep='\t', float_format='%.8e', index=False, header=False)

        #and the foils
        events = []
        for i, run in enumerate([1,2,3]):
    
            evs = pd.read_csv(f'{subj.model_dir}/memory_run-0{run}/sm_events/foil_{con}.txt', sep='\t', header=None)
            if evs.sum(axis=1)[0] == 0: 
                pass
            else:
                #correct for start time
                evs[0] += i*620
                events.append(evs)
        out_events = pd.concat(events).reset_index(drop=True)
        out_events.to_csv( f'{sm_out}/foil_{con}.txt',
                sep='\t', float_format='%.8e', index=False, header=False)

def copy_out_lvl3():
    
    dest = '/scratch/05426/ach3377/source_memory_group_glm/zmaps'
    mkdir(dest)

    for cope in range(1,31):

        in_map = f'/scratch/05426/ach3377/source_memory_group_glm/cope{cope}.gfeat/cope1.feat/stats/zstat1.nii.gz'
        out_map = f'{dest}/cope{cope}_zmap.nii.gz'

        os.system(f'cp {in_map} {out_map}')

def cluster_stats():

    base_dir = '/scratch/05426/ach3377/source_memory_group_glm'
    out_dir = f'{base_dir}/cluster_stats';mkdir(out_dir)
    
    for c in range(1,31):
        cope_dir = f'{base_dir}/cope{c}.gfeat'
        out_cope_dir = f'{out_dir}/cope{c}';mkdir(out_cope_dir)

        in_zmap = f'{cope_dir}/cope1.feat/stats/zstat1.nii.gz'
        zmap = f'{cope_dir}/cope1.feat/thresh_zmap.nii.gz'
        mask = f'{cope_dir}/cope1.feat/mask.nii.gz'
        cope = f'{cope_dir}/cope1.feat/stats/cope1.nii.gz'

        stats = pd.read_csv(f'{cope_dir}/cope1.feat/stats/smoothness',header=None)
        smooth = np.float(stats[0][0].split()[1])
        volume = np.int(stats[0][1].split()[1])

        mask_cmd = f'fslmaths {in_zmap} -mas {mask} {zmap}'

        cluster_cmd = f'$FSLDIR/bin/cluster -i {zmap} \
                        -t 2.575 \
                        --othresh={out_cope_dir}/thresh_zmap \
                        -o cluster_mask_zmap \
                        --connectivity=26 \
                        --mm \
                        --olmax={out_cope_dir}/lmax_zmap_std.txt \
                        --scalarname=Z \
                        -p 0.05 \
                        -d {smooth} \
                        --volume={volume} \
                        -c {cope} \
                        > {out_cope_dir}/cluster_table.txt'

        for cmd in [mask_cmd,cluster_cmd]:
            os.system(cmd)

def switch_to_lvl2_in_afni(x):
    in_str = x.split('/')[-1]
    sub  = x.split('/')[7]
    cope = in_str[4:in_str.find('.nii.gz')]

    return f'/scratch/05426/ach3377/fc-bids/derivatives/model/{sub}/all_memory_runs/source_memory.feat/reg_standard/stats/cope{cope}.nii.gz'

def hack_lvl2():
    subj = bids_meta(4)

    dest = '/scratch/05426/ach3377/hack_lvl2';mkdir(dest)

    in_lvl2 = f'{subj.feat_dir}/source_memory_lvl2.gfeat'
    out_lvl2 = f'{dest}/lvl2.gfeat'
    os.system(f'cp -r {in_lvl2} {out_lvl2}')

    for run in [1,2,3]:
        in_run = f'{subj.model_dir}/memory_run-0{run}/source_memory.feat'
        out_run = f'{dest}/run{run}.feat'
        os.system(f'cp -r {in_run} {out_run}')


def afni_commands():
    RE = {  0:{'name':'Intercept',
                'dof':[1,1579],
                'thr':6.6509},

            1:{'name':'Encode',
                'dof':[2,1579],
                'thr':4.6186},

            2:{'name':'Condition',
                'dof':[1,1579],
                'thr':6.6509},

            3:{'name':'Response',
                'dof':[2,1579],
                'thr':4.6186},

            4:{'name':'Encode_Condition',
                'dof':[2,1579],
                'thr':4.6186},

            5:{'name':'Encode_Response',
                'dof':[4,1579],
                'thr':3.3329},

            6:{'name':'Condition_Response',
                'dof':[2,1579],
                'thr':4.6186},

            7:{'name':'Encode_Condition_Response',
               'dof':[4,1579],
                'thr':3.3329},
            }
    for effect in RE:
        os.system(f"3dcalc -a SM_LME+tlrc'[{effect}]' \
                -expr 'fift_t2p(a,{RE[effect]['dof'][0]},{RE[effect]['dof'][1]})' \
                -datum float \
                -prefix {RE[effect]['name']}_pmap.nii.gz")

        os.system(f"fslmaths {RE[effect]['name']}_pmap.nii.gz -mul -1 -add 1 -mas ../standard/MNI152NLin2009cAsym_T1_1mm_brain_mask.nii.gz {RE[effect]['name']}_1-pmap.nii.gz")

def afni_fwhmx(sub):
    inputs = pd.read_csv('sm_events/afni_dataTable.txt',sep=' ').InputFile

    subj = bids_meta(sub)
    for run in [1,2,3]:
        feat_dir = f'{subj.model_dir}/memory_run-0{run}/source_memory.feat'
        resid = f'{feat_dir}/stats/res4d.nii.gz'
        resid_std = f'{feat_dir}/reg_standard/stats/res4d.nii.gz'
        
        os.chdir(feat_dir)
        
        reg_cmd = f'flirt -ref {feat_dir}/reg/standard \
                          -in {resid} \
                          -out {resid_std} \
                          -applyxfm \
                          -init {feat_dir}/reg/example_func2standard.mat \
                          -interp trilinear \
                          -datatype float'
        
        est_noise = f'3dFWHMx -mask {std_2009_brain_mask} \
                              -input {resid_std} \
                              -acf >> noise_estimates.txt'

        for cmd in [reg_cmd, est_noise]:
            os.system(cmd)
# for sub in xcl_sub_args:
#     os.system(f"echo singularity run --cleanenv $SCRATCH/bids-apps/neurosft.simg python $HOME/gPPI/wrap_glm_utils.py -s {sub} >> jobs/afni_fwhm_job.txt")
# os.system('launch -N 48 -n 48 -J smooth -s jobs/afni_fwhm_job.txt -m achennings@utexas.edu -p normal -r 2:00:00 -A LewPea_MRI_Analysis')

def collect_fwhm():
    _df = {}
    for sub in all_sub_args:#xcl_sub_args:
        subj = bids_meta(sub)
        _df[sub] = {}
        for run in [1,2,3]:
            _df[sub][run] = pd.read_csv(f'{subj.model_dir}/memory_run-0{run}/basic_model.feat/noise_estimates.txt',sep=' '
                            ).loc[0].dropna().values
    df = pd.DataFrame.from_dict(_df).unstack().reset_index()
    est = df[0].values.mean(axis=0)
    est.tofile('sm_events/basic_model_mean_smooth.txt',sep=' ', format='%s')

def clustsim(file='sm_events/basic_model_mean_smooth.txt'):
    est = np.loadtxt(file)
    cmd1 = 'export OMP_NUM_THREADS=48'
    cmd2 = f'3dClustSim -OKsmallmask \
                       -mask $SCRATCH/standard/gm_1mm_thr.nii.gz \
                       -acf {est[0]} {est[1]} {est[2]} \
                       > sm_events/basic_model_clustsim_output.txt' 
    
    script = 'sm_events/clustsim_script.txt'
    os.system(f'rm {script}')

    for cmd in [cmd1,cmd2]:
        os.system(f"echo {cmd} >> {script}")

    jobfile = f'/home1/05426/ach3377/gPPI/jobs/clustsim_job.txt'
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
                       -r 5:00:00 \
                       -A LewPea_MRI_Analysis')

def clusterize(folder,file,thr1,thr2):
    # os.chdir('../../Desktop/3dLME_results')
    # header = parse_AFNI_header('SM_LME+tlrc.head')
    os.chdir(f'../../Desktop/{folder}')
    header = parse_AFNI_header(f'{file}.HEAD')


    names = header['BRICK_LABS'].replace(':','_'
                                ).replace(' F','_F'
                                ).replace(' t', '_t'
                                ).split('~')

    for i, name in enumerate(names):
        if 'Intercept' in name:
            pass
        elif '_F' in name or '_t' in name:
            cmap = f'{name}_cluster_map.nii.gz';os.system(f'rm {cmap}')
            ctxt = f'{name}_cluster.txt';os.system(f'rm {ctxt}')
            where = f'{name}_where.txt';os.system(f'rm {where}')

            if '_F' in name:
                side = '1sided RIGHT_TAIL'
                thr = thr1
            elif '_t' in name:
                side = '2sided'
                thr = thr2

            cmd = f"3dClusterize \
                      -inset {file} \
                      -ithr {i} \
                      -mask ../standard/gm_1mm_thr.nii.gz \
                      -{side} p=0.005 \
                      -clust_nvox {thr} \
                      -NN 3 \
                      -pref_map {cmap} \
                      > {ctxt}"
            os.system(cmd)
            

            if os.path.exists(cmap):
                w_cmd = f"whereami -coord_file {ctxt}'[1,2,3]' > {where}"
                os.system(w_cmd)

def brik_break(file):
    from nibabel.brikhead import parse_AFNI_header
    import os
    import nibabel as nib

    header = parse_AFNI_header(f'{file}.HEAD')

    tmp = '.nii'
    os.system(f'3dAFNItoNIFTI {file} -prefix {tmp}')
    img = nib.load(f'{file.split("+")[0]}.nii')

    names = header['BRICK_LABS'].split('~')

    assert len(names) == img.shape[-1]

    for i, name in enumerate(names):
        nib.save(img.slicer[:,:,:,:,i],f'{name}.nii.gz')

    os.system(f'rm {file.split("+")[0]}.nii')

def extract_clusters(folder):
    from nilearn.image import get_data, new_img_like
    
    os.chdir(folder)
    effects = [i.split('_where')[0] for i in os.listdir() if 'where' in i]

    for eff in effects:
        out_dir = f'{eff}_cluster_masks';mkdir(out_dir)

        cmap = f'{eff}_cluster_map.nii.gz'
        cluster_data = get_data(cmap)
        for clust in np.unique(cluster_data):
            if clust > 0:
                mask_data = np.zeros(cluster_data.shape)
                mask_data[np.where(cluster_data == clust)] = 1

                out_mask = new_img_like(nib.load(cmap),mask_data.astype(int),copy_header=True)
                nib.save(out_mask,f'{out_dir}/cluster_{clust}_mask.nii.gz')

def extract_pe(effects=['Response','Condition','Condition_Response'],):
    from nilearn.image import get_data
    #create output df
    phases = ['baseline','acquisition','extinction']
    runs = [1,2,3]

    for eff in ['Condition_Response']:
        mask_dir = f'{SCRATCH}3dLME_results/{eff}_cluster_masks'
        clusters = {i.split('_')[1].split('-')[1]:{'num':i.split('_')[1].split('-')[0],'mask':get_data(f'{mask_dir}/{i}')} for i in os.listdir(mask_dir)}
        
        df = pd.DataFrame(columns=['pe'],index=pd.MultiIndex.from_product([xcl_sub_args,phases,consp,phases,runs,list(clusters.keys())],
                            names=['subject','encode_phase','condition','source_memory','run','roi']))

        #use the inputs from the LME data table
        inputs = pd.read_csv('sm_events/afni_dataTable.txt',sep=' '
            ).set_index(['Subj','Sess','Encode','Condition','Response'])

        for sub in xcl_sub_args:
            print(sub)
            subj = bids_meta(sub)
            for encode in phases:
                print(encode)
                for con in consp:
                    for response in phases:
                        for run in runs:
                            #check if it needs to be run or if missing
                            try:
                                file = inputs.loc[(sub,f'memory_run-0{run}',encode,con,response),'InputFile']
                                pe_img = get_data(file)
                                for roi in clusters:
                                    df.loc[(sub,encode,con,response,run,roi),'pe'] = apply_mask(mask=clusters[roi]['mask'],target=pe_img).mean()
                            except:
                                pass
        df.to_csv(f'sm_events/{eff}_extracted_pe.csv')

def basic_model_reg_smooth(sub):

    subj = bids_meta(sub)
    for run in [1,2,3]:
        feat_dir = f'{subj.model_dir}/memory_run-0{run}/basic_model.feat'
        resid = f'{feat_dir}/stats/res4d.nii.gz'
        resid_std = f'{feat_dir}/reg_standard/stats/res4d.nii.gz'
        
        os.chdir(feat_dir)
        
        featreg_cmd = f'featregapply {feat_dir}' 

        res_reg_cmd = f'flirt -ref {feat_dir}/reg/standard \
                          -in {resid} \
                          -out {resid_std} \
                          -applyxfm \
                          -init {feat_dir}/reg/example_func2standard.mat \
                          -interp trilinear \
                          -datatype float'
        
        est_noise = f'3dFWHMx -mask {std_2009_brain_mask} \
                              -input {resid_std} \
                              -acf >> noise_estimates.txt'

        for cmd in [featreg_cmd, res_reg_cmd, est_noise]:
            os.system(cmd)
    # for sub in all_sub_args:
    #    os.system(f"echo singularity run --cleanenv $SCRATCH/bids-apps/neurosft.simg python $HOME/gPPI/wrap_glm_utils.py -s {sub} >> jobs/basic_model_reg_smooth_job.txt")
    # os.system('launch -N 48 -n 48 -J smooth -s jobs/basic_model_reg_smooth_job.txt -m achennings@utexas.edu -p normal -r 2:00:00 -A LewPea_MRI_Analysis')


basic_pe_map = {   'baseline':{'CSp':1,
                               'CSm':2},
                'acquisition':{'CSp':4,
                               'CSm':5},
                 'extinction':{'CSp':7,
                               'CSm':8},
                       'foil':{'CSp':10,
                               'CSm':11}}
def basic_model_dataTable():
    data_table = pd.DataFrame(columns=['Subj','Sess','Group','Encode','Condition','InputFile'])
    phases = ['baseline','acquisition','extinction','foil']
    memory_phases = ['memory_run-01','memory_run-02','memory_run-03']
    
    for sub in all_sub_args:
        subj = bids_meta(sub)
        for encode_phase in phases:
            for con in consp:
                for mem_phase in memory_phases:
                    data_table = data_table.append({
                        'Subj':sub,
                        'Sess':mem_phase,
                        'Group':'healthy' if sub < 100 else 'ptss',
                        'Encode':encode_phase,
                        'Condition':con,
                        'InputFile':f'/scratch/05426/ach3377/fc-bids/derivatives/model/{subj.fsub}/{mem_phase}/basic_model.feat/reg_standard/stats/cope{basic_pe_map[encode_phase][con]}.nii.gz'},
                        ignore_index=True)
    data_table.to_csv('sm_events/memory_dataTable.txt',index=False,sep=' ')
    data_table = data_table[data_table.Encode.isin(phases[:-1])]
    data_table.to_csv('sm_events/memory_no_foils_dataTable.txt',index=False,sep=' ')

