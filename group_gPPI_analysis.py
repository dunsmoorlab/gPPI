from fg_config import *
from nilearn.image import get_data

def apply_mask(mask=None,target=None):

    coor = np.where(mask == 1)
    values = target[coor]
    if values.ndim > 1:
        values = np.transpose(values) #swap axes to get feature X sample
    return values

def mask_and_extract():
    ROIS = ['mOFC','dACC','lh_hpc','rh_hpc','lh_amyg','rh_amyg']
    COPES = ['acq_ext','ext_acq']
    PHASES = ['memory_run-01','memory_run-02','memory_run-03']
    
    df = pd.DataFrame(columns=['conn'], index=pd.MultiIndex.from_product([ROIS,COPES,ROIS,all_sub_args,PHASES], names=['seed','cope','target','subject','phase']))

    cons = {'acq_ext':13,'ext_acq':14}

    for sub in all_sub_args:
        subj = bids_meta(sub)
        print(sub)
        for phase in PHASES:
            print(phase)
            phase_dir = os.path.join(subj.model_dir,phase)

            for seed in ROIS:
                print(seed)
                seed_dir = os.path.join(phase_dir,seed,'source.feat')

                for target in ROIS:

                    if 'h' in target:
                        in_mask = os.path.join(subj.masks,'%s.nii.gz'%(target))
                    else:
                        in_mask = os.path.join(subj.masks,'%s_mask.nii.gz'%(target))

                    phase_mask = os.path.join(seed_dir,'%s_mask.nii.gz'%(target))
                    fsl_mask = os.path.join(seed_dir,'mask.nii.gz')
                    
                    # os.system('fslmaths %s -mas %s %s'%(in_mask,fsl_mask,phase_mask))

                    mask_dat = get_data(phase_mask)

                    for cope in COPES:

                        cope_data = get_data(os.path.join(seed_dir,'stats','zstat%s.nii.gz'%(cons[cope])))

                        extracted = apply_mask(mask=mask_dat,target=cope_data)

                        df.loc[(seed,cope,target,sub,phase),'conn'] = extracted.mean()

df = df.groupby(['seed','cope','target','subject']).mean()


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


