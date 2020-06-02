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
df = df.drop(columns='phase')
sns.barplot(data=df,x='roi',y='beta',hue='model')


df = df.set_index(['cope','roi','model','subject'])
df = (df.loc['csp'] - df.loc['csm'])
