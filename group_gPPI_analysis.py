from fg_config import *
from nilearn.image import get_data
def cluster():
    cope_dir = os.path.join(SCRATCH,'group_gPPI_out')
    ROIS = ['mOFC','dACC','lh_hpc','rh_hpc','lh_amyg','rh_amyg']
    groups = ['healthy','ptsd']
    cons = ['ext_csp_acq_csp','acq_csp_ext_csp']

    df = pd.DataFrame(columns=['p','peak'],index=pd.MultiIndex.from_product([ROIS,cons,ROIS,groups],names=['seed','con','target','group']))
    for seed in ROIS:
        seed_dir = os.path.join(cope_dir,seed)
        
        for con in cons:
            con_dir = os.path.join(seed_dir,con)
            
            for roi in ROIS:
                mask = os.path.join(group_masks,'%s_group_mask.nii.gz'%(roi))
            
                for group in groups:
                    zimg = os.path.join(con_dir,'%s_0.nii.gz'%(group))
                    otxt = os.path.join(con_dir,'%s_%s_smooth.txt'%(roi,group))
                    # os.system('smoothest --zstat=%s --mask=%s >> %s'%(zimg,mask,otxt))

                    zmas = os.path.join(con_dir,'%s_0_%s.nii.gz'%(group,roi))
                    # os.system('fslmaths %s -mas %s %s'%(zimg, mask, zmas))

                    smooth = pd.read_csv(otxt,sep='\t',header=None)
                    dlh = smooth[0][0][4:]
                    vol = smooth[0][1][7:]

                    ocls = os.path.join(con_dir, '%s_%s_cluster.txt'%(roi,group))
                    os.system('rm %s'%(ocls))
                    os.system('cluster --zstat=%s --zthresh=1.96 -p 0.05 --dlh=%s --volume=%s >> %s'%(zmas,dlh,vol,ocls))

                    clust = pd.read_csv(ocls,sep='\t',header=0)
                    if clust.shape[0] == 0: 
                        pass
                    elif clust.shape[0] > 0:
                        print(seed,con,roi,group,'\n\n',clust,'\n\n')
                        df.loc[(seed,con,roi,group),'p'] = clust.P[0]
                        df.loc[(seed,con,roi,group),'peak'] = clust.MAX[0]
                        if clust.shape[0] > 1:
                            print('OOPSIE')

def apply_mask(mask=None,target=None):

    coor = np.where(mask == 1)
    values = target[coor]
    if values.ndim > 1:
        values = np.transpose(values) #swap axes to get feature X sample
    return values

def mask_and_extract():
    ROIS = ['sgACC','rACC','lh_hpc','rh_hpc','lh_amyg','rh_amyg']
    day2_copes = {'acq_csp_csm':6,
                  'ext_csp_csm':9,
                  'acq_ext':13,
                  'ext_acq':14}
    
    day1_copes = {'csp':1,
                  'csm':2,
                  'csp_csm':3}
    
    mem_phases = ['memory_run-01','memory_run-02','memory_run-03']
    encode_phases = ['baseline','acquisition','extinction']
    all_phases = encode_phases+mem_phases


    mem_df = pd.DataFrame(columns=['conn'], index=pd.MultiIndex.from_product([ROIS,day2_copes,ROIS,all_sub_args,mem_phases], names=['seed','cope','target','subject','phase']))
    encode_df = pd.DataFrame(columns=['conn'], index=pd.MultiIndex.from_product([ROIS,day1_copes,ROIS,all_sub_args,encode_phases], names=['seed','cope','target','subject','phase']))

    for sub in all_sub_args:
        subj = bids_meta(sub)
        print(sub)
        
        for target in ROIS:
            print(target)        
            
            #load target mask here and extract for all seeds and phases
            if 'h' in target:
                in_mask = os.path.join(subj.masks,'%s.nii.gz'%(target))
            else:
                in_mask = os.path.join(subj.masks,'%s_mask.nii.gz'%(target))
            mask_dat = get_data(in_mask)

            for seed in ROIS:

                for phase in all_phases:

                    seed_dir = os.path.join(subj.model_dir,phase,seed,'source.feat')

                    if phase in mem_phases:
                        copes = day2_copes
                        df = mem_df
                    elif phase in encode_phases:
                        copes = day1_copes
                        df = encode_df
                    
                    for cope in copes:

                        cope_data = get_data(os.path.join(seed_dir,'stats','cope%s.nii.gz'%(copes[cope])))

                        extracted = apply_mask(mask=mask_dat,target=cope_data)

                        df.loc[(seed,cope,target,sub,phase),'conn'] = extracted.mean()
    
    mem_df = mem_df.reset_index()
    mem_df.conn = mem_df.conn.astype(float)
    mem_df = mem_df.groupby(['seed','cope','target','subject']).mean().reset_index()
    mem_df.to_csv('extracted_mem_gPPI.csv',index=False)

    encode_df = encode_df.reset_index()
    encode_df.conn = encode_df.conn.astype(float)
    encode_df.to_csv('extracted_encode_gPPI.csv',index=False)

def group_gPPI_clean():
    import os
    from fg_config import SCRATCH, mkdir
    import nibabel as nib
    from nilearn.image import index_img

    mem_copes = {'baseline_csp':1,
             'baseline_csm':2,
             'baseline_csp_csm':3,
             'acquisition_csp':4,
             'acquisition_csm':5,
             'acquisition_csp_csm':6,
             'extinction_csp':7,
             'extinction_csm':8,
             'extinction_csp_csm':9,
             'foil_csp':10,
             'foil_csm':11,
             'foil_csp_csm':12,
             'acq_csp_ext_csp':13,
             'ext_csp_acq_csp':14}

    encode_copes = {'csp':1,
                    'csm':2,
                    'csp_csm':3}

    encode_phases = ['baseline','acquisition','extinction']


    stats = {'healthy_0':1,
             '0_healthy':2,
             'ptsd_0':3,
             '0_ptsd':4,
             'healthy_ptsd':5,
             'ptsd_healthy':6,
             'healhty_ptsd_0':7,
             '0_healthy_ptsd':8}

    # for roi in ['sgACC','rACC','lh_amyg','rh_amyg','lh_hpc','rh_hpc']:
    for roi in ['rh_hpc']:
        wd = os.path.join(SCRATCH,'group_gPPI',roi)
        mem_od = os.path.join(SCRATCH,'group_gPPI_out',roi,'retrieval');mkdir(mem_od)
        encode_od = os.path.join(SCRATCH,'group_gPPI_out',roi,'encoding');mkdir(encode_od)

        for cope in mem_copes:
            out = os.path.join(mem_od,cope);mkdir(out)
            for stat in stats:
                infile = os.path.join(wd,'cope%s++++++.gfeat'%(mem_copes[cope]),'cope1.feat','stats','zstat%s.nii.gz'%(stats[stat]))
                outfile = os.path.join(out,'%s.nii.gz'%(stat))
                os.system('cp %s %s'%(infile, outfile))

        for phase in encode_phases:
            for cope in encode_copes:
                out = os.path.join(encode_od,phase,cope);mkdir(out)
                for stat in stats:
                    infile = os.path.join(wd,phase+'.gfeat','cope%s.feat'%(encode_copes[cope]),'stats','zstat%s.nii.gz'%(stats[stat]))
                    outfile = os.path.join(out,'%s.nii.gz'%(stat))
                    os.system('cp %s %s'%(infile, outfile))


def pconvert(p):
    if p < .001:
        return '***'
    elif p < .01:
        return '**'
    elif p < .05:
        return '*'
    elif p < .1:
        return '~'
    else:
        return ''



def graph_gPPI():
    from fg_config import lgroup

    import matplotlib.pyplot as plt
    import pingouin as pg
    import seaborn as sns

    ROIS = ['rACC','sgACC','lh_hpc','rh_hpc','lh_amyg','rh_amyg']
    # COPES = ['acq_ext','ext_acq']
    COPES = ['ext_acq']
    groups = ['healthy','ptsd']

    df = pd.read_csv('extracted_mem_gPPI.csv')
    df = df.groupby(['seed','cope','target','subject']).mean().reset_index()
    df['group'] = df.subject.apply(lgroup) 
    df = df.set_index(['cope','group','seed','target']).sort_index()

    stats = pd.DataFrame(columns=['t','p'],index=pd.MultiIndex.from_product([groups,ROIS,ROIS],names=['group','seed','target']))
    gstats = pd.DataFrame(columns=['t','p'],index=pd.MultiIndex.from_product([ROIS,ROIS],names=['seed','target']))
    
    for seed in ROIS:
        for group in groups:
            for target in ROIS:
                tres = pg.ttest(df.loc[('ext_acq',group,seed,target),'conn'].values,0,tail='two-sided')
                stats.loc[(group,seed,target)][['t','p']] = tres.loc['T-test'][['T','p-val']]

        # gres = pg.ttest(df.loc[('ext_acq','healthy',seed,target),'conn'].values,df.loc[('ext_acq','ptsd',seed,target),'conn'].values,paired=False)
        # gstats.loc[(seed,target)][['t','p']] = gres.loc['T-test'][['T','p-val']]

    
    # mask = np.zeros([len(ROIS),len(ROIS)])
    # mask[np.diag_indices_from(mask)] = True 

    # fig, (gax, gcbar) = plt.subplots(2,2,gridspec_kw={'height_ratios':(.9,.05),'hspace':.5})
    # for j, cope in enumerate(COPES):
        
    #     gt = gstats.loc[(cope),'t'].unstack(level=-1).astype(float).loc[ROIS][ROIS]
    #     gp = gstats.loc[(cope),'p'].apply(pconvert).unstack(level=-1).astype(str).loc[ROIS][ROIS]

    #     sns.heatmap(gt,mask=mask,ax=gax[j],square=True,
    #                 annot=gp,fmt='',cmap='PRGn',center=0,vmin=-3,vmax=3,
    #                 cbar_ax=gcbar[j],cbar_kws={'orientation':'horizontal'})
    #     gax[j].set_title(cope + '_group_comp')

    fig1, (ax1, cbar1) = plt.subplots(2,2,gridspec_kw={'height_ratios':(.9,.05),'hspace':.5})
    fig2, (ax2, cbar2) = plt.subplots(2,2,gridspec_kw={'height_ratios':(.9,.05),'hspace':.5})

    for i, group in enumerate(groups):
        t = stats.loc[(group),'t'].unstack(level=-1).astype(float).loc[ROIS][ROIS]
        p = stats.loc[(group),'p'].apply(pconvert).unstack(level=-1).astype(str).loc[ROIS][ROIS]

        # sns.heatmap(t,mask=mask,ax=ax[i],square=True,
        #             annot=p,fmt='',cmap='PRGn',center=0,vmin=-3,vmax=3,
        #             cbar_ax=cbar[i],cbar_kws={'orientation':'horizontal'})
        # # ax[i].set_title(group + '_' + cope)

        pfc_targ_t = t.loc[('rh_hpc', 'lh_hpc', 'rh_amyg', 'lh_amyg'),['rACC','sgACC']].T
        pfc_targ_p = p.loc[('rh_hpc', 'lh_hpc', 'rh_amyg', 'lh_amyg'),['rACC','sgACC']].T

        sns.heatmap(pfc_targ_t,ax=ax1[i],annot=pfc_targ_p,square=True,
                    fmt='',cmap='PRGn',center=0,vmin=-3,vmax=3,
                    cbar_ax=cbar1[i],cbar_kws={'orientation':'horizontal'})
        ax1[i].set_title(group + ' ext vs. acq')

        
        pfc_seed_t = t.loc[('rACC','sgACC'),['rh_hpc', 'lh_hpc', 'rh_amyg', 'lh_amyg']].T
        pfc_seed_p = p.loc[('rACC','sgACC'),['rh_hpc', 'lh_hpc', 'rh_amyg', 'lh_amyg']].T
        
        sns.heatmap(pfc_seed_t,ax=ax2[i],annot=pfc_seed_p,square=True,
                    fmt='',cmap='PRGn',center=0,vmin=-3,vmax=3,
                    cbar_ax=cbar2[i],cbar_kws={'orientation':'horizontal'})
        ax2[i].set_title(group + ' ext vs. acq')


    ####encoding!####
    edf = pd.read_csv('extracted_encode_gPPI.csv')
    edf = edf.set_index(['cope','phase','seed','target','subject'])
    edf = (edf.loc['csp','extinction'] - edf.loc['csp','acquisition']).reset_index()
    edf['group'] = edf.subject.apply(lgroup)
    edf = edf.set_index(['group','seed','target'])

    estats = pd.DataFrame(columns=['t','p'],index=pd.MultiIndex.from_product([groups,ROIS,ROIS],names=['group','seed','target']))
    
    for seed in ROIS:
        for group in groups:
            for target in ROIS:
                etres = pg.ttest(edf.loc[(group,seed,target),'conn'].values,0,tail='two-sided')
                estats.loc[(group,seed,target)][['t','p']] = etres.loc['T-test'][['T','p-val']]
    
    fig3, (ax3, cbar3) = plt.subplots(2,2,gridspec_kw={'height_ratios':(.9,.05),'hspace':.5})
    fig4, (ax4, cbar4) = plt.subplots(2,2,gridspec_kw={'height_ratios':(.9,.05),'hspace':.5})

    for i, group in enumerate(groups):
        t = estats.loc[(group),'t'].unstack(level=-1).astype(float).loc[ROIS][ROIS]
        p = estats.loc[(group),'p'].apply(pconvert).unstack(level=-1).astype(str).loc[ROIS][ROIS]

        # sns.heatmap(t,mask=mask,ax=ax[i],square=True,
        #             annot=p,fmt='',cmap='PRGn',center=0,vmin=-3,vmax=3,
        #             cbar_ax=cbar[i],cbar_kws={'orientation':'horizontal'})
        # # ax[i].set_title(group + '_' + cope)

        pfc_targ_t = t.loc[('rh_hpc', 'lh_hpc', 'rh_amyg', 'lh_amyg'),['rACC','sgACC']].T
        pfc_targ_p = p.loc[('rh_hpc', 'lh_hpc', 'rh_amyg', 'lh_amyg'),['rACC','sgACC']].T

        sns.heatmap(pfc_targ_t,ax=ax3[i],annot=pfc_targ_p,square=True,
                    fmt='',cmap='PRGn',center=0,vmin=-3,vmax=3,
                    cbar_ax=cbar3[i],cbar_kws={'orientation':'horizontal'})
        ax3[i].set_title(group + ' ext vs. acq')

        
        pfc_seed_t = t.loc[('rACC','sgACC'),['rh_hpc', 'lh_hpc', 'rh_amyg', 'lh_amyg']].T
        pfc_seed_p = p.loc[('rACC','sgACC'),['rh_hpc', 'lh_hpc', 'rh_amyg', 'lh_amyg']].T
        
        sns.heatmap(pfc_seed_t,ax=ax4[i],annot=pfc_seed_p,square=True,
                    fmt='',cmap='PRGn',center=0,vmin=-3,vmax=3,
                    cbar_ax=cbar4[i],cbar_kws={'orientation':'horizontal'})
        ax4[i].set_title(group + ' ext vs. acq')




