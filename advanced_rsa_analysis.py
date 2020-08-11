###############fear to extinction similarity, day1 and day2
c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)


cae = pd.DataFrame(index=pd.MultiIndex.from_product([memcon,cons,rois,sub_args],names=['memory_phase','condition','roi','subject']))
pae = pd.DataFrame(index=pd.MultiIndex.from_product([memcon,cons,rois,p_sub_args],names=['memory_phase','condition','roi','subject']))
for mem in memcon:
    for con in cons:
        for roi in rois:
            for sub in subs:
                cae.loc[(mem,con,roi,sub_args[sub]),'rsa'] = c.mats[roi][sub,slice3[con]['acquisition'][mem],slice3[con]['extinction'][mem]].mean()
                pae.loc[(mem,con,roi,p_sub_args[sub]),'rsa'] = p.mats[roi][sub,slice3[con]['acquisition'][mem],slice3[con]['extinction'][mem]].mean()
ae = pd.concat((cae,pae)).reset_index()
ae['group'] = ae.subject.apply(lgroup)

bn = ae.set_index('roi')
bn = bn.loc[bn_rois].reset_index()
bn = bn[bn.condition == 'CS+']
bn = bn.set_index('subject').drop([20,120]).reset_index(
        ).set_index(['group','roi','memory_phase','subject']).sort_index()

sns.catplot(data=bn.reset_index(),x='roi',y='rsa',hue='memory_phase',row='group',kind='bar')


stats = pd.DataFrame(columns=['w','p','pc','cles'],index=pd.MultiIndex.from_product([['healthy','ptsd'],rois],names=['group','roi']))


for group in ['healthy','ptsd']:
    for roi in rois:
        wres = pg.wilcoxon(bn.loc[(group,roi,'encoding'),'rsa'],bn.loc[(group,roi,'retrieval'),'rsa'],tail='two-sided').values
        stats.loc[(group,roi),'w']    = wres[0,0]
        stats.loc[(group,roi),'p']    = wres[0,2]
        stats.loc[(group,roi),'cles'] = wres[0,-1]
    stats.loc[(group,),'pc'] = pg.multicomp(list(stats.loc[(group),'p'].values),method='fdr_bh')[1]
# stats.p = 1 - stats.ps
# stats.p = stats.p.apply(lambda x: 0 if x < .95 else x)
stats['p_mask'] = stats.pc.apply(lambda x: 0 if x >.05 else 1)
stats['cles_disp'] = stats.cles * stats.p_mask

