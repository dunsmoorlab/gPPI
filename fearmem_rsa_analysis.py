from index_rsm import get_square

#get the source memory data into the recognition memory df
sm_convert = {1:'baseline',2:'acquisition',3:'extinction'}
sm_df = pd.read_csv('../fearmem/cleaned_full_sm.csv').set_index('subject')
sm_df.source_memory = sm_df.source_memory.apply(lambda x: sm_convert[x])
for sub in xcl_sub_args:
    subj = bids_meta(sub)
    sm_dat = sm_df.loc[sub].copy().set_index('stimulus')
    mem_dat = pd.read_csv(f'rsa_results/{subj.fsub}/reordered_mem_labels.csv')#.drop(columns='index')
    mem_dat['source_memory'] = ''
    for i in mem_dat.index:
        if mem_dat.loc[i,'encode_phase'] == 'foil': pass
        else: mem_dat.loc[i,'source_memory'] = sm_dat.loc[mem_dat.loc[i,'stimulus'],'source_memory']
    # mem_dat.to_csv(f'rsa_results/{subj.fsub}/reordered_mem_labels.csv',index=False)
    mkdir(f'sm_events/{subj.fsub}')
    mem_dat.to_csv(f'sm_events/{subj.fsub}/sm_events.csv',index=False)


#now do the rsa
c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)

gmats = {'healthy':c.mem_mats,
         'ptsd':p.mem_mats}

seeds = ['hc_tail', 'hc_body', 'hc_head', 'amyg_bla', 'amyg_cem']
rois = ['rACC','sgACC'] + seeds
# rois = bn_rois

sm_con = ['baseline','acquisition']
df = pd.DataFrame({'rsa':0.0},index=pd.MultiIndex.from_product([sm_con,cons,rois,xcl_sub_args],names=['response_phase','condition','roi','subject']))
for roi in rois:
    for con in cons:
        for group in groups:
            for s in subs:
                sub = subjects[group][s]
                if sub in xcl_sub_args:
                    subj = bids_meta(sub)
                    mem_dat = pd.read_csv(f'rsa_results/{subj.fsub}/reordered_mem_labels.csv')
                    baseline_correct = mem_dat[mem_dat.trial_type == con][mem_dat.encode_phase == 'extinction'][mem_dat.source_memory == 'extinction'].index
                    baseline_acquisition = mem_dat[mem_dat.trial_type == con][mem_dat.encode_phase == 'extinction'][mem_dat.source_memory == 'acquisition'].index
                    acquisition_correct = mem_dat[mem_dat.trial_type == con][mem_dat.encode_phase == 'acquisition'][mem_dat.source_memory == 'acquisition'].index
                    

                    b_correct_acq_correct = gmats[group][sub][roi][baseline_correct][:,acquisition_correct].mean()
                    b_acq_acq_correct = gmats[group][sub][roi][baseline_acquisition][:,acquisition_correct].mean()

                    df.loc[('baseline',con,roi,sub),'rsa'] = b_correct_acq_correct
                    df.loc[('acquisition',con,roi,sub),'rsa'] = b_acq_acq_correct

df = df.reset_index()
# df['group'] = df.subject.apply(lgroup)
sns.barplot(data=df[df.roi=='rACC'],x='condition',y='rsa',hue='response_phase')#,kind='bar',col='group')
sns.swarmplot(data=df[df.roi=='rACC'],x='condition',y='rsa',hue='response_phase',dodge=True,color='black')#,kind='bar',col='group')

sns.catplot(data=df,x='condition',y='rsa',hue='response_phase',col='roi',kind='bar')
                    # square = get_square(gmats[group], roi, sub, s, mem_slices, con, 'baseline')

#baseline to acq similarity regardless of source memory
df = pd.DataFrame({'rsa':0.0},index=pd.MultiIndex.from_product([cons,rois,all_sub_args],names=['condition','roi','subject']))
for roi in rois:
    for con in cons:
        for group in groups:
            for s in subs:
                sub = subjects[group][s]
                if sub in all_sub_args:
                    subj = bids_meta(sub)
                    subj = bids_meta(sub)
                    mem_dat = pd.read_csv(f'rsa_results/{subj.fsub}/reordered_mem_labels.csv')
                    
                    baseline = mem_dat[mem_dat.trial_type == con][mem_dat.encode_phase == 'extinction'].index
                    acquisition = mem_dat[mem_dat.trial_type == con][mem_dat.encode_phase == 'acquisition'].index
                    
                    ba = gmats[group][sub][roi][baseline][:,acquisition].mean()
                    df.loc[(con,roi,sub),'rsa'] = ba
df = df.reset_index()
df.roi = df.roi.apply(pfc_rename)
sns.catplot(data=df,x='condition',y='rsa',col='roi',kind='strip')



#     for el_phase in el:
#         for con in cons:
#             for roi in rois: print(roi)
#                 for sub in subs:                                        
#                     cae.loc[(mem,el_phase,con,roi,sub_args[sub]),'rsa'] = get_square(c.mats,roi,sub_args[sub],sub,split_slices,con,'late_acquisition',mem,con,'%s_extinction'%(el_phase),mem)
#                     pae.loc[(mem,el_phase,con,roi,p_sub_args[sub]),'rsa'] = get_square(p.mats,roi,p_sub_args[sub],sub,split_slices,con,'late_acquisition',mem,con,'%s_extinction'%(el_phase),mem)
# ae = pd.concat((cae,pae)).reset_index()
# ae['group'] = ae.subject.apply(lgroup)
encoding_phases = ['baseline','acquisition','extinction']
memory_phases = ['memory_run-01','memory_run-02','memory_run-03']
consp = ['CSp','CSm']

def sm_glm_events(sub):
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
                        out_events = pd.DataFrame({'onset':0.0,'duration':0.0,'PM':0.0},index=[0])
                        os.system(f'echo {sub}\t{mem_phase}\t{encode_phase}\t{con}\t{sm_resp} >> sm_events/missing_evs.txt')
                    else:                   
                        out_events = smdat[['onset','duration','PM']].copy()
                    
                    out_events.to_csv( os.path.join(mem_phase_dir, f'{encode_phase}_{con}_{sm_resp}.txt'),
                        sep='\t', float_format='%.8e', index=False, header=False)

            #and the foils
            mem_phase_dat.loc[('foil',con),['onset','duration','PM']].to_csv( os.path.join(mem_phase_dir, f'foil_{con}.txt'),
                        sep='\t', float_format='%.8e', index=False, header=False)

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

        replacements[search_for] = replace_with

        # with open(os.path.join(gPPI_codebase,'feats','%s.fsf'%(template))) as infile: 
    with open(template) as infile:
        with open('feats/rep_test.fsf', 'w') as outfile:
            for line in infile:
                for src, target in replacements.items():
                    line = line.replace(src, target)
                outfile.write(line)
