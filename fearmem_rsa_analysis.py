from index_rsm import get_square

c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)

seeds = ['hc_tail', 'hc_body', 'hc_head', 'amyg_bla', 'amyg_cem']
rois = ['rACC'] + seeds

sm_con = ['baseline','acquisition']
cdf = pd.DataFrame(index=pd.MultiIndex.from_product([sm_con,cons,rois,sub_args],names=['memory_phase','el','condition','roi','subject']))
pdf = pd.DataFrame(index=pd.MultiIndex.from_product([sm_con,cons,rois,p_sub_args],names=['memory_phase','el','condition','roi','subject']))

    for el_phase in el:
        for con in cons:
            for roi in rois:
                print(roi)
                for sub in subs:                                        
                    cae.loc[(mem,el_phase,con,roi,sub_args[sub]),'rsa'] = get_square(c.mats,roi,sub_args[sub],sub,split_slices,con,'late_acquisition',mem,con,'%s_extinction'%(el_phase),mem)
                    pae.loc[(mem,el_phase,con,roi,p_sub_args[sub]),'rsa'] = get_square(p.mats,roi,p_sub_args[sub],sub,split_slices,con,'late_acquisition',mem,con,'%s_extinction'%(el_phase),mem)
ae = pd.concat((cae,pae)).reset_index()
ae['group'] = ae.subject.apply(lgroup)

for sub in smt_sub_args:
    csub = sub_args[]