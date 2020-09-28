from index_rsm import get_square
from roi_rsa import *
from fg_config import *

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
rois = ['rACC','sgACC','thalamus_clst','RSP_clst','dACC_clst','lOFC_clst'] + seeds
# rois = bn_rois

'''BASELINE TO ACQUISITION_CORRECT, SPLIT BY SM_RESPONSE'''
sm_con = ['baseline','acquisition']
df = pd.DataFrame({'rsa':0.0},index=pd.MultiIndex.from_product([sm_con,cons,rois,xcl_sub_args],names=['response_phase','condition','roi','subject']))
for roi in rois:
    for con in cons:
        for group in groups:
            for s in subs:
                sub = subjects[group][s]
                if sub in xcl_sub_args:
                    subj = bids_meta(sub)
                    mem_dat = pd.read_csv(f'sm_events/{subj.fsub}/sm_events.csv')
                    acquisition_correct = mem_dat[mem_dat.trial_type == con][mem_dat.encode_phase == 'acquisition'][mem_dat.source_memory == 'acquisition'].index
                    baseline_correct = mem_dat[mem_dat.trial_type == con][mem_dat.encode_phase == 'baseline'][mem_dat.source_memory == 'baseline'].index
                    baseline_acquisition = mem_dat[mem_dat.trial_type == con][mem_dat.encode_phase == 'baseline'][mem_dat.source_memory == 'acquisition'].index
                    

                    b_correct_b_acq = gmats[group][sub][roi][baseline_correct][:,acquisition_correct].mean()
                    b_acq_acq_correct = gmats[group][sub][roi][baseline_acquisition][:,acquisition_correct].mean()

                    df.loc[('baseline',con,roi,sub),'rsa'] = b_correct_b_acq
                    df.loc[('acquisition',con,roi,sub),'rsa'] = b_acq_acq_correct

df = df.reset_index()
# df['group'] = df.subject.apply(lgroup)
spal = list((wes_palettes['Darjeeling1'][-1],wes_palettes['Darjeeling1'][0],wes_palettes['Darjeeling1'][1],))
# sns.swarmplot(data=df[df.roi=='rACC'],x='condition',y='rsa',hue='response_phase',dodge=True,color='black')#,kind='bar',col='group')
# sns.catplot(data=df,x='condition',y='rsa',hue='response_phase',col='roi',kind='bar')
                    # square = get_square(gmats[group], roi, sub, s, mem_slices, con, 'baseline')
for roi in rois:
    print(roi)
    print(pg.wilcoxon(df.rsa[df.response_phase == 'baseline'][df.condition == 'CS+'][df.roi==roi],
                        df.rsa[df.response_phase == 'acquisition'][df.condition == 'CS+'][df.roi==roi]))
    print('\n')
'''(CS+B_A vs. CS+B_B) to CS+A_A vmPFC finding'''
fig, ax = plt.subplots(figsize=(8,6))
sns.barplot(data=df[df.roi=='sgACC'],x='condition',y='rsa',hue='response_phase',palette=[spal[0],spal[1]],ax=ax)#,kind='bar',col='group')
ax.set_ylabel('Similarity to acquisition "correct"')
ax.set_title('vmPFC')

'''CS+B_A vs. CS-B_A bla & head'''
fig, ax = plt.subplots(figsize=(8,6))
sns.barplot(data=df[df.roi.isin(['amyg_bla','hc_head'])], x='roi', y='rsa', hue='condition', palette=[spal[1],'white'], edgecolor=spal[1])
ax.set_title('CS+B_A vs. CS-B_A')

'''BASELINE TO ACQ. SIMILARITY REGARDLESS OF SOURCE MEMORY'''
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


'''B_B VS. B_A ERS'''
df = pd.concat((c.df,p.df)
            ).set_index(['subject','stimulus']
            ).sort_index(
            ).loc[xcl_sub_args]
sm = pd.read_csv('../fearmem/cleaned_full_sm.csv')
stims = sm.stimulus.unique()
sm = sm.set_index(['subject','stimulus']).sort_index()

df['source_memory'] = 0
for sub in xcl_sub_args:
    for stim in stims:
        df.loc[(sub, stim), 'source_memory'] = sm.loc[(sub, stim), 'source_memory']

phase_convert = {1:'baseline',2:'acquisition',3:'extinction'}
spal = list((wes_palettes['Darjeeling1'][-1],wes_palettes['Darjeeling1'][0],wes_palettes['Darjeeling1'][1],))

df.source_memory = df.source_memory.apply(lambda x: phase_convert[x])
df = df.groupby(['source_memory','trial_type','encode_phase','roi','subject']).mean()
df = df.reset_index()

#filter out rois here
pfc = df[df.roi.isin(rois)]
pfc.roi = pd.Categorical(pfc.roi,categories=rois)
sns.catplot(data=pfc[pfc.encode_phase=='extinction'], x='trial_type', y='rsa', 
            hue='source_memory', col='roi', palette=spal, hue_order=['baseline','acquisition','extinction'], kind='bar')
statdf = pfc.set_index(['roi','encode_phase','trial_type','source_memory','subject']).sort_index()
for roi in rois:
    s = pg.wilcoxon(statdf.loc[(roi,'extinction','CS+','extinction'),'rsa'],
            statdf.loc[(roi,'extinction','CS+','acquisition'),'rsa'])
    print(roi,'\n',s)

fig, ax = plt.subplots(figsize=(8,6))
sns.barplot(data=statdf.loc['sgACC','baseline',slice('CS+','CS-'),'acquisition'].reset_index()
     ,x='trial_type', y='rsa', palette=[spal[1],'white'], edgecolor=spal[1])
ax.set_title('CS+B_A vs. CS-B_A reinstatement')
ax.set_xlabel('condition')



sns.boxplot(data=statdf.loc['hc_tail','baseline',slice('CS+','CS-'),'acquisition'].reset_index()
            ,y='rsa',hue='trial_type')
pg.plot_paired(data=statdf.loc['hc_tail','baseline',slice('CS+','CS-'),'acquisition'].reset_index()
                ,dv='rsa', within='trial_type', subject='subject')