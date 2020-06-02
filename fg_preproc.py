from fg_config import *
import nibabel as nib
from nilearn.image import mean_img, get_data, threshold_img, new_img_like, clean_img
from subprocess import Popen

class fmriprep_preproc():

    def __init__(self,sub):

        self.subj = bids_meta(sub)
        self.space = 'MNI152NLin2009cAsym'
        # self.space = 'T1w'
    
    def be_t1(self): #brain extract T1w

        # t1 = os.path.join(self.subj.prep_dir,'anat','%s_desc-preproc_T1w.nii.gz'%(self.subj.fsub))
        # t1_mask = os.path.join(self.subj.prep_dir,'anat','%s_desc-brain_mask.nii.gz'%(self.subj.fsub))
        t1 = os.path.join(self.subj.prep_dir,'anat','%s_space-%s_desc-preproc_T1w.nii.gz'%(self.subj.fsub,self.space))
        t1_mask = os.path.join(self.subj.prep_dir,'anat','%s_space-%s_desc-brain_mask.nii.gz'%(self.subj.fsub,self.space))

        os.system('cp %s %s'%(t1,self.subj.t1))
        os.system('cp %s %s'%(t1_mask,self.subj.t1_mask))
        os.system('fslmaths %s -mas %s %s'%(self.subj.t1,self.subj.t1_mask,self.subj.t1_brain))

    def refvol(self): #create reference bold, bold mask, and mask the ref

        refs = [nib.load( os.path.join(self.subj.prep_dir,folder[0],file) )
                        for folder in os.walk(self.subj.prep_dir)
                        for file in folder[2]
                        if '%s_boldref.nii.gz'%(self.space) in file]
        ref = mean_img(refs)
        nib.save(ref,self.subj.refvol)
        
        masks = [nib.load( os.path.join(self.subj.prep_dir,folder[0],file) )
                        for folder in os.walk(self.subj.prep_dir)
                        for file in folder[2]
                        if '%s_desc-brain_mask.nii.gz'%(self.space) in file]
        mask = mean_img(masks)
        mask = threshold_img(mask,.1)
        nib.save(mask,self.subj.refvol_mask)

        os.system('fslmaths %s -bin %s'%(self.subj.refvol_mask,self.subj.refvol_mask))

        os.system('fslmaths %s -mas %s %s'%(self.subj.refvol,self.subj.refvol_mask,self.subj.refvol_brain))
    
    def be_bold(self): #mask all of the preproc bold files with the new mask

        for folder in os.walk(self.subj.prep_dir):
            for file in folder[2]:
                if '%s_desc-preproc_bold.nii.gz'%(self.space) in file:
                    os.system('fslmaths %s -mas %s %s'%(os.path.join(self.subj.prep_dir,folder[0],file), self.subj.refvol_mask, os.path.join(self.subj.func,file)))

    def bbreg(self):
        # reg = 'export SUBJECTS_DIR=%s; bbregister --s %s --mov %s --bold --reg %s'%(fs_dir,self.subj.fsub,self.subj.refvol_brain,self.subj.fs_regmat)
        # v2v = 'export SUBJECTS_DIR=%s; mri_vol2vol --subject %s --targ %s --mov %s --reg %s --nearest --inv --o %s'%(fs_dir,self.subj.fsub,os.path.join(self.subj.fs_dir,'mri','aparc+aseg.mgz'),self.subj.refvol_brain,self.subj.fs_regmat,self.subj.faa)
        # mask = 'fslmaths %s -mas %s %s'%(self.subj.faa,self.subj.refvol_mask,self.subj.faa)
        # for cmd in [reg,v2v,mask]: os.system(cmd)

        os.system('fslmaths %s/ses-1/func/%s_ses-1_task-extinction_space-%s_desc-aparcaseg_dseg.nii.gz -mas %s %s'%(self.subj.prep_dir,self.subj.fsub,self.space,self.subj.refvol_mask,self.subj.faa))

    def fs_mask(self):

        rois = {
                 'mOFC':{'l':1014,
                         'r':2014},

                 'amyg':{'l':18,
                         'r':54},

                 'hpc':{'l':17,
                        'r':53},

                 'dACC':{'l':1002,
                         'r':2002},

                 'ins':{'l':1035,
                        'r':2035}
                 }

        for roi in rois:
            for hemi in rois[roi]:
                #make the hemisphere specific masks
                out = os.path.join(self.subj.masks,'%sh_%s.nii.gz'%(hemi,roi))
                thr = rois[roi][hemi]
                cmd = 'fslmaths %s -thr %s -uthr %s -bin %s'%(self.subj.faa,thr,thr,out)
                os.system(cmd)
            #make the joint mask
            l = os.path.join(self.subj.masks,'lh_%s.nii.gz'%(roi))
            r = os.path.join(self.subj.masks,'rh_%s.nii.gz'%(roi))
            out = os.path.join(self.subj.masks,'%s_mask.nii.gz'%(roi))
            cmd = 'fslmaths %s -add %s -bin %s'%(l,r,out)
            os.system(cmd)

    def brainnetome(self):
        pass


    def hca_mask(self):

        fsn2t1w = os.path.join(self.subj.prep_dir,'anat','%s_from-fsnative_to-T1w_mode-image_xfm.txt'%(self.subj.fsub))
        t1w2std = os.path.join(self.subj.prep_dir,'anat','%s_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5'%(self.subj.fsub)) 

        t1w_in_t1w_space = os.path.join(self.subj.prep_dir,'anat','%s_desc-preproc_T1w.nii.gz'%(self.subj.fsub))
            
        hc_rois = {'hc_head':232,
                   'hc_body':231,
                   'hc_tail':226}

        amyg_rois = {'amyg_bla':[7001,7003],
                     'amyg_cem':[7004,7010]}


        for hemi in ['l','r']:  

            in_label  = os.path.join(self.subj.fs_dir,'mri',hemi+'h.hippoAmygLabels-T1.v21.HBT.FSvoxelSpace.mgz')
            working_label = os.path.join(self.subj.masks,hemi+'h_hca_mask.nii.gz')
            
            t1_convert  = 'antsApplyTransforms -i %s \
                                               -t %s \
                                               -r %s \
                                               -o %s \
                                               -n MultiLabel'%(in_label,fsn2t1w,t1w_in_t1w_space,working_label)

            std_convert = 'antsApplyTransforms -i %s \
                                               -t %s \
                                               -r %s \
                                               -o %s \
                                               -n MultiLabel'%(working_label,t1w2std,std_2009_brain,working_label)

            resamp = 'ResampleImage 3 %s %s 65x77x65 1'%(working_label,working_label)

            for cmd in [t1_convert,std_convert,resamp]: os.system(cmd)


            #do the whole hippocamppus
            hemi_hpc = 'fslmaths %s -thr 226 -uthr 232 -bin %s'%(working_label,os.path.join(self.subj.masks,hemi+'h_hpc.nii.gz'))
            os.system(hemi_hpc)

            #and then head/body/tail
            for roi in hc_rois:
                thr = str(hc_rois[roi])
                out = os.path.join(self.subj.masks,hemi+'h_%s.nii.gz'%(roi))
                hc_cmd = ['fslmaths', working_label,
                              '-thr', thr,
                             '-uthr', thr,
                              '-bin', out]
                Popen(hc_cmd).wait()

            #do the whole amygdala
            hemi_amyg = 'fslmaths %s -thr 7001 -uthr 7010 -bin %s'%(working_label,os.path.join(self.subj.masks,hemi+'h_amyg.nii.gz'))
            os.system(hemi_amyg)

            #and then bla/cem
            for roi in amyg_rois:
                lthr = str(amyg_rois[roi][0])
                uthr = str(amyg_rois[roi][1])
                out = os.path.join(self.subj.masks,hemi+'h_%s.nii.gz'%(roi))
                amyg_cmd = ['fslmaths', working_label,
                                    '-thr', lthr,
                                   '-uthr', uthr,
                                    '-bin', out ]
                Popen(amyg_cmd).wait()


        for roi in ['hpc','amyg','hc_head','hc_body','hc_tail','amyg_bla','amyg_cem']:
            cmd = ['fslmaths', os.path.join(self.subj.masks,'rh_'+roi+'.nii.gz'),
                       '-add', os.path.join(self.subj.masks,'lh_'+roi+'.nii.gz'),
                       '-bin', os.path.join(self.subj.masks,roi+'.nii.gz')]
            Popen(cmd).wait()

    def group_mask(self):
        
        #masks = ['sgACC','rSMA','rACG'] 
        # masks = ['rACC']
        #masks = ['A32sg','A32p','A24cd','A24rv','A14m','A11m','A13','A10m','A9m','A8m','A6m',]
        masks = ['rACC','sgACC']        

        for roi in masks:
            in_mask = os.path.join(group_masks,'%s_group_mask.nii.gz'%(roi))

            out_mask = os.path.join(self.subj.masks,'%s_mask.nii.gz'%(roi))
            
            os.system('flirt -in %s -ref %s -applyxfm -init %s -out %s -interp nearestneighbour'%(
                    in_mask, self.subj.refvol_brain, self.subj.std2ref, out_mask))

    def fsl_reg(self):

        os.system('flirt -in %s -ref %s -dof 12 -omat %s'%(self.subj.refvol_brain, std_2009_brain, self.subj.ref2std))
        os.system('convert_xfm -omat %s -inverse %s'%(self.subj.std2ref,self.subj.ref2std))

        os.system('flirt -in %s -ref %s -dof 12 -omat %s'%(self.subj.refvol_brain, std_2009_brain_3mm, self.subj.ref2std3))
        os.system('convert_xfm -omat %s -inverse %s'%(self.subj.std32ref,self.subj.ref2std3))

    def apply_reg_betas(self):
        for task in tasks:
            inbeta = os.path.join(self.subj.beta,'%s_beta.nii.gz'%(task))
            outbeta = os.path.join(self.subj.beta,'%s_beta_std.nii.gz'%(task))
            os.system('flirt -in %s -ref %s -applyxfm -init %s -out %s'%(inbeta,std_2009_brain_3mm,self.subj.ref2std3,outbeta))

    def apply_reg_rsa(self):
        for phase in ['baseline','acquisition','extinction']:
            for con in ['CSp','CSm']:
                infile = os.path.join(self.subj.weights,'%s_%s.nii.gz'%(phase,con))
                outfile = os.path.join(self.subj.weights,'%s_%s_std.nii.gz'%(phase,con))
                os.system('flirt -in %s -ref %s -applyxfm -init %s -out %s'%(infile,std_2009_brain_3mm,self.subj.ref2std3,outfile))
    
    def denoise(self):

        for task in tasks:
            # if 'memory' in task:
            if tasks[task]['ses'] == 1 and self.subj.num in [105,106]:
                tr = 2.23
            else:
                tr = 2
            
            if task == 'localizer_run-02' and self.subj.num == 107:
                pass
            else:

                inbold = os.path.join(self.subj.func,'%s_ses-%s_task-%s_space-%s_desc-preproc_bold.nii.gz'%(self.subj.fsub,tasks[task]['ses'],task,self.space))
                outbold = os.path.join(self.subj.func,'%s_ses-%s_task-%s_space-%s_desc-preproc_denoised_bold.nii.gz'%(self.subj.fsub,tasks[task]['ses'],task,self.space))
                confounds = os.path.join(self.subj.model_dir,task,'confounds.txt')

                tmp = clean_img(nib.load(inbold), detrend=False, standardize=True, confounds=confounds,
                            low_pass=None, high_pass=None, t_r=tr, ensure_finite=True, mask_img=None)

                nib.save(tmp,outbold)

        # standard      = '/work/IRC/ls5/opt/apps/fsl-5.0.10/data/standard/MNI152_T1_1mm_brain'
        # standard_head = '/work/IRC/ls5/opt/apps/fsl-5.0.10/data/standard/MNI152_T1_1mm'
        # standard_mask = '/work/IRC/ls5/opt/apps/fsl-5.0.10/data/standard/MNI152_T1_1mm_brain_mask_dil'
        
        # os.system('epi_reg --epi=%s --t1=%s --t1brain=%s --out=%s'%(self.subj.refvol_brain,self.subj.t1,self.subj.t1_brain,self.subj.ref2t1[:-4]))

        # os.system('flirt -in %s -ref %s -omat %s \
        #                 -cost corratio -dof 12 -searchrx -90 90 \
        #                 -searchry -90 90 -searchrz -90 90 -interp trilinear'%(self.subj.t1_brain,standard,self.subj.t12std))

        # os.system('fnirt --in=%s --aff=%s --cout=%s \
        #                  --config=T1_2_MNI152_2mm \
        #                  --ref=%s --refmask=%s --warpres=10,10,10'%(self.subj.t1,self.subj.t12std,self.subj.t12std_warp,standard_head,standard_mask))

        # os.system('applywarp --ref=%s --in=%s \
        #                      --out=%s --warp=%s \
        #                      --premat=%s --interp=nn'%(standard,self.subj.faa,self.subj.saa,self.subj.t12std_warp,self.subj.ref2t1))
def brainnetome_group():
    from nilearn.image import resample_to_img

    rois = {
             'A32sg':[187,188],
             'A32p':[179,180],
             'A24cd':[183,184],
             'A24rv':[177,178],
             'A14m':[41,42],
             'A11m':[47,48],
             'A13':[49,50],
             'A10m':[13,14],
             'A9m':[11,12],
             'A8m':[1,2],
             'A6m':[9,10]
             }
    gmask = os.path.join(SCRATCH,'fc-bids','derivatives','group_masks')
    atlas = os.path.join(gmask,'BNA-maxprob-thr0-1mm.nii.gz')
    gprob = os.path.join(gmask,'std_gm_thr.nii.gz')
    
    std = nib.load(std_2009_brain)


    for roi in rois:

        out_mask = os.path.join(gmask,'%s_group_mask.nii.gz'%(roi))
        mask_cmd = 'fslmaths %s -thr %s -uthr %s -bin %s'%(atlas,rois[roi][0],rois[roi][1],out_mask)
        thr_cmd = 'fslmaths %s -mas %s -bin %s'%(gprob,out_mask,out_mask)

        os.system(mask_cmd)
        resamp = resample_to_img(nib.load(out_mask),std,interpolation='nearest')
        nib.save(resamp,out_mask)
        os.system(thr_cmd)


def std_space_check():
    from nilearn.image import get_data
    pdir = 'D:\\fc-bids\\derivatives\\fmriprep'
    for sub in all_sub_args:
        print(sub)
        subj = bids_meta(sub)
        sub_dir = os.path.join(pdir,subj.fsub)
        aparcs = glob('%s/ses-*/func/%s_ses-*_task-**_space-MNI152NLin2009cAsym_desc-aparcaseg_dseg.nii.gz'%(sub_dir,subj.fsub))
        for a in aparcs:
            assert np.array_equal(get_data(aparcs[0]),get_data(a))

def group_std_masks():
    from nilearn.image import get_data, new_img_like, resample_img
    import nibabel as nib
    
    saa = {}
    
    thr = {'mOFC':     [1014,2014],
            'dACC':     [1002,2002],
            'lh_amyg':  [18],
            'rh_amyg':  [54],
            'lh_hpc':   [17],
            'rh_hpc':   [53]}

    masks = {}
    for roi in thr: masks[roi] = {}


    for sub in all_sub_args:
        subj = bids_meta(sub)
        print(sub)
        # os.system('flirt -in %s -out %s -ref %s -applyxfm -init %s -interp nearestneighbour'%(subj.faa,subj.saa,std_2009_brain,subj.ref2std))

        saa[sub] = get_data(subj.faa)

        for roi in thr:
            if len(thr[roi]) == 2:
                l = np.zeros(saa[sub].shape)
                l[np.where(saa[sub] == thr[roi][0])] = 1

                r = np.zeros(saa[sub].shape)
                r[np.where(saa[sub] == thr[roi][1])] = 1

                m = l+r
            
            else:
                m = np.zeros(saa[sub].shape)
                m[np.where(saa[sub] == thr[roi][0])] = 1

            masks[roi][sub] = m

    for roi in masks:
        masks[roi] = np.array([masks[roi][i] for i in masks[roi]])

        masks[roi] = np.sum(masks[roi],axis=0) / 48
        
        masks[roi] = new_img_like(subj.faa,masks[roi],copy_header=False)
        # masks[roi].header['pixdim'] = [1.,3., 3., 3., 0., 0., 0., 0.,]
        # masks[roi] = new_img_like(std_2009_brain,masks[roi],copy_header=False)
        std = nib.load(std_2009_brain)
        masks[roi] = resample_img(masks[roi],target_affine=std.affine,target_shape=std.shape,interpolation='nearest')
        
        # nib.save(masks[roi],os.path.join(group_masks,'%s_3mm_raw.nii.gz'%(roi)))
        # nib.save(masks[roi],os.path.join(group_masks,'%s_1mm_raw.nii.gz'%(roi)))
        nib.save(masks[roi],os.path.join(group_masks,'%s_1mm.nii.gz'%(roi)))
        
        os.system('fslmaths %s -bin %s'%(os.path.join(group_masks,'%s_1mm.nii.gz'%(roi)), os.path.join(group_masks,'%s_group_mask.nii.gz'%(roi))))

    for hemi in ['l','r']:
        amyg = get_data(os.path.join(group_masks,'%sh_amyg_1mm.nii.gz'%(hemi)))
        hpc  = get_data(os.path.join(group_masks,'%sh_hpc_1mm.nii.gz'%(hemi)))

        a_bin = np.where(amyg > hpc, 1, 0)
        h_bin = np.where(hpc > amyg, 1, 0)

        a_out = new_img_like(std_2009_brain,a_bin,copy_header=False)
        h_out = new_img_like(std_2009_brain,h_bin,copy_header=False)

        nib.save(a_out,os.path.join(group_masks,'%sh_amyg_group_mask.nii.gz'%(hemi)))
        nib.save(h_out,os.path.join(group_masks,'%sh_hpc_group_mask.nii.gz'%(hemi)))


def copy_events_confounds():
    from bids_model import bids_events
    dest = os.path.join(SCRATCH,'preproc')
    for sub in all_sub_args:
        subj = bids_meta(sub)
        sub_dest = os.path.join(dest,subj.fsub)
        mkdir(sub_dest)
        events = os.path.join(sub_dest,'events')
        confounds = os.path.join(sub_dest,'confounds')
        mkdir(events)
        mkdir(confounds)
        for task in tasks:
            c = os.path.join(subj.model_dir,task,'confounds.txt')
            c_out = os.path.join(confounds,'%s_task-%s_confounds.txt'%(subj.fsub,task))
            os.system('cp %s %s'%(c,c_out))

            e = bids_events(sub).phase_events(task)
            e.to_csv(os.path.join(events, '%s_task-%s_events.tsv'%(subj.fsub,task)), sep='\t', index=False)

def backup_betas():
    #from fg_config import *
    dest = os.path.join(SCRATCH,'fc-bids','derivatives','beta_backup')
    for sub in all_sub_args[1:]:
        subj = bids_meta(sub)
        indir = subj.beta
        outdir = os.path.join(dest,subj.fsub)
        mkdir(outdir)
        os.system('cp -r %s %s'%(indir, outdir))


    source = os.path.join(SCRATCH,'fc-bids','derivatives','beta_backup')
    for sub in all_sub_args:
        subj = bids_meta(sub)
        indir = os.path.join(source,subj.fsub,'lss_betas')
        outdir = subj.preproc_dir
        os.system('rm -r %s'%(subj.beta))
        os.system('cp -r %s %s'%(indir,outdir))
