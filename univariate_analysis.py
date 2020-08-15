from fg_config import *
from nilearn.image import get_data
def apply_mask(mask=None,target=None):

    coor = np.where(mask == 1)
    values = target[coor]
    if values.ndim > 1:
        values = np.transpose(values) #swap axes to get feature X sample
    return values

# ROIS = ['sgACC','rACC','lh_hpc','rh_hpc','lh_amyg','rh_amyg']
ROIS = ['rACC','sgACC']
day2_copes = {'acq_csp_csm':6,
              'ext_csp_csm':9,
              'acq_ext':13,
              'ext_acq':14}

day1_copes = {'csp':1,
              'csm':2}
              # 'csp_csm':3,
              # 'csm_csp':4}

mem_phases = ['memory_run-01','memory_run-02','memory_run-03']
# encode_phases = ['baseline','acquisition','extinction']
encode_phases = ['acquisition']
models = ['reg','dn']
all_phases = encode_phases+mem_phases

encode_df = pd.DataFrame(columns=['beta'], index=pd.MultiIndex.from_product([ROIS,day1_copes,all_sub_args,encode_phases,models], names=['roi','cope','subject','phase','model']))
# mem_df = pd.DataFrame(columns=['conn'], index=pd.MultiIndex.from_product([ROIS,day2_copes,ROIS,all_sub_args,mem_phases], names=['seed','cope','target','subject','phase']))

for sub in all_sub_args:
    subj = bids_meta(sub)
    print(sub)
    
    for roi in ROIS:
        print(roi)        
        #load target mask here and extract for all seeds and phases
        try:
            mask_dat = get_data(os.path.join(subj.masks,'%s.nii.gz'%(roi)))
        except ValueError:
            mask_dat = get_data(os.path.join(subj.masks,'%s_mask.nii.gz'%(roi)))

        for phase in encode_phases:#all_phases

            for model in models:

                if model == 'reg':

                    feat_dir = os.path.join(subj.model_dir,phase,'%s_%s_reg_gPPI.feat'%(subj.fsub,phase))
                else:
                    feat_dir = os.path.join(subj.model_dir,phase,'%s_%s_gPPI.feat'%(subj.fsub,phase))

                if phase in mem_phases:
                    copes = day2_copes
                    df = mem_df
                elif phase in encode_phases:
                    copes = day1_copes
                    df = encode_df
                
                for cope in copes:

                    cope_data = get_data(os.path.join(feat_dir,'stats','cope%s.nii.gz'%(copes[cope])))

                    extracted = apply_mask(mask=mask_dat,target=cope_data)

                    df.loc[(roi,cope,sub,phase,model),'beta'] = extracted.mean()

# mem_df = mem_df.reset_index()
# mem_df.conn = mem_df.conn.astype(float)
# mem_df = mem_df.groupby(['seed','cope','target','subject']).mean().reset_index()
# mem_df.to_csv('extracted_mem_gPPI.csv',index=False)

encode_df = encode_df.reset_index()
encode_df.beta = encode_df.beta.astype(float)
encode_df.to_csv('univariate_gPPI_comp.csv',index=False)

df = pd.read_csv('univariate_gPPI_comp.csv')
df.roi = df.roi.apply(lambda x: 'vmPFC' if x == 'sgACC' else 'dACC')
df.model = df.model.apply(lambda x: 'denoised before' if x == 'dn' else 'confounds in model')
df.cope = df.cope.apply(lambda x: 'CS+' if x == 'csp' else 'CS-')
df = df.drop(columns='phase')
df = df.set_index(['cope','roi','model','subject'])


pg.ttest(df.loc[('CS+','dACC','confounds in model'),'beta'],df.loc[('CS-','dACC','confounds in model'),'beta'],paired=True)
pg.ttest(df.loc[('CS+','dACC','denoised before'),'beta'],df.loc[('CS-','dACC','denoised before'),'beta'],paired=True)

pg.ttest(df.loc[('CS+','vmPFC','confounds in model'),'beta'],df.loc[('CS-','vmPFC','confounds in model'),'beta'],paired=True)
pg.ttest(df.loc[('CS+','vmPFC','denoised before'),'beta'],df.loc[('CS-','vmPFC','denoised before'),'beta'],paired=True)


sns.catplot(data=df.reset_index(),x='roi',y='beta',hue='cope',col='model',kind='bar',sharey=False)

df = (df.loc['CS+'] - df.loc['CS-'])




#########actual extraction###########
seeds = ['hc_tail','hc_body','hc_head','amyg_bla','amyg_cem'] + bn_rois

# day2_copes = {'acq_csp_csm':6,
#               'ext_csp_csm':9,
#               'acq_ext':13,
#               'ext_acq':14}
copes = {'CS+B':1,
         'CS-B':2,
         'CS+A':3,
         'CS-A':4,
         'CS+E':5,
         'CS-E':6,
         'CS+F':7,
         'CS-F':8}

mem_phases = ['memory_run-01','memory_run-02','memory_run-03']

df = pd.DataFrame(columns=['beta'], index=pd.MultiIndex.from_product([seeds,copes,all_sub_args,mem_phases], names=['roi','cope','subject','phase']))

for sub in all_sub_args:
    subj = bids_meta(sub)
    print(sub)
    
    for roi in seeds:
        print(roi)        
        #load target mask here and extract for all seeds and phases
        try:
            mask_dat = get_data(os.path.join(subj.masks,'%s.nii.gz'%(roi)))
        except ValueError:
            mask_dat = get_data(os.path.join(subj.masks,'%s_mask.nii.gz'%(roi)))

        for phase in mem_phases:#all_phases

            feat_dir = os.path.join(subj.model_dir,phase,'%s_%s_gPPI.feat'%(subj.fsub,phase))

            for cope in copes:

                cope_data = get_data(os.path.join(feat_dir,'stats','cope%s.nii.gz'%(copes[cope])))

                extracted = apply_mask(mask=mask_dat,target=cope_data)

                df.loc[(roi,cope,sub,phase),'beta'] = extracted.mean()


##########analyze it############
df = pd.read_csv('univariate_gPPI.csv')
df = df.groupby(['cope','subject','roi']).mean()

beta = (df.loc['CS+E'] - df.loc['CS-E']).rename(columns={'beta':'ext_csp_csm'})
beta['acq_csp_csm'] = (df.loc['CS+A'] - df.loc['CS-E'])
beta['bsl_csp_csm'] = (df.loc['CS+B'] - df.loc['CS-B'])
beta['ext_acq'] = (df.loc['CS+E'] - df.loc['CS+A'])
beta['foil_csp_csm'] = (df.loc['CS+F'] - df.loc['CS-F'])
beta = beta.reset_index().melt(id_vars=['roi','subject'],
                        value_vars=['ext_csp_csm', 'acq_csp_csm','bsl_csp_csm','ext_acq','foil_csp_csm'],
                        value_name='beta',var_name='cope')
beta['group'] = beta.subject.apply(lgroup)
beta = beta.set_index('subject').drop([20,120]).reset_index()

# sns.catplot(data=beta[beta.cope != 'ext_acq'],x='roi',y='beta',hue='cope',row='group',kind='bar')
# sns.catplot(data=beta[beta.cope == 'ext_acq'],x='roi',y='beta',hue='group',kind='bar')
# sns.catplot(data=beta[beta.cope == 'foil_csp_csm'],x='roi',y='beta',hue='group',kind='bar')



copes = {'ext_acq':     ['CS+E','CS+A'],
         'ext_csp_csm': ['CS+E','CS-E'],
         'acq_csp_csm': ['CS+A','CS-A']}

stats = stats = pd.DataFrame(columns={'w':0.0,'p':0.0,'p_fdr':0.0,'p_mask':0.0,'cles':0.0},
                     index=pd.MultiIndex.from_product([groups,copes,bn_rois],
                     names=['group','cope','roi'])
                     ).sort_index()

df = df.reset_index()
df['group'] = df.subject.apply(lgroup)
df = df.set_index(['cope','group','roi','subject'])
for group in groups:
    for cope in copes:
        con1, con2 = copes[cope][0], copes[cope][1]
        for roi in bn_rois:
            wres = pg.wilcoxon(df.loc[(con1,group,roi),'beta'],df.loc[(con2,group,roi),'beta'])
            stats.loc[(group,cope,roi),['w','p','cles']] = wres[['W-val','p-val','CLES']].values

        stats.loc[(group,cope),'p_fdr'] = pg.multicomp(list(stats.loc[(group,cope),'p'].values),method='fdr_bh')[1]


#sub cortical stats
stats = stats = pd.DataFrame(columns={'w':0.0,'p':0.0,'p_fdr':0.0,'p_mask':0.0,'cles':0.0},
                     index=pd.MultiIndex.from_product([groups,copes,seeds],
                     names=['group','cope','roi'])
                     ).sort_index()

df = df.reset_index()
df['group'] = df.subject.apply(lgroup)
df = df.set_index(['cope','group','roi','subject'])
for group in groups:
    for cope in copes:
        con1, con2 = copes[cope][0], copes[cope][1]
        for roi in seeds:
            wres = pg.wilcoxon(df.loc[(con1,group,roi),'beta'],df.loc[(con2,group,roi),'beta'])
            stats.loc[(group,cope,roi),['w','p','cles']] = wres[['W-val','p-val','CLES']].values

        stats.loc[(group,cope),'p_fdr'] = pg.multicomp(list(stats.loc[(group,cope),'p'].values),method='fdr_bh')[1]
