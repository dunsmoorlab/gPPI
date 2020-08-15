from index_rsm import get_square 

cmaps = {'baseline':'Greys',
         'acquisition':'Purples',
         'early_extinction':sns.light_palette('seagreen',as_cmap=True),
         'extinction':'Greens'}

###############fear to extinction similarity, day1 and day2###############
#replicating labar 2020

c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)

el = ['early','late']
seeds = ['hc_tail', 'hc_body', 'hc_head', 'amyg_bla', 'amyg_cem']
rois = bn_rois + seeds

cae = pd.DataFrame(index=pd.MultiIndex.from_product([memcon,el,cons,rois,sub_args],names=['memory_phase','el','condition','roi','subject']))
pae = pd.DataFrame(index=pd.MultiIndex.from_product([memcon,el,cons,rois,p_sub_args],names=['memory_phase','el','condition','roi','subject']))
for mem in memcon:
    for el_phase in el:
        for con in cons:
            for roi in rois:
                print(roi)
                for sub in subs:                                        
                    cae.loc[(mem,el_phase,con,roi,sub_args[sub]),'rsa'] = get_square(c.mats,roi,sub_args[sub],sub,split_slices,con,'late_acquisition',mem,con,'%s_extinction'%(el_phase),mem)
                    pae.loc[(mem,el_phase,con,roi,p_sub_args[sub]),'rsa'] = get_square(p.mats,roi,p_sub_args[sub],sub,split_slices,con,'late_acquisition',mem,con,'%s_extinction'%(el_phase),mem)
ae = pd.concat((cae,pae)).reset_index()
ae['group'] = ae.subject.apply(lgroup)


ae = ae.set_index('subject').drop([20,120]).reset_index()
        # ).set_index(['group','roi','memory_phase','subject']).sort_index()

sns.catplot(data=ae[ae.condition=='CS+'][ae.memory_phase=='encoding'],x='roi',y='rsa',hue='el',row='group',kind='bar')
sns.catplot(data=ae[ae.condition=='CS+'][ae.memory_phase=='retrieval'],x='roi',y='rsa',hue='el',row='group',kind='bar')

# sns.catplot(data=ae[ae.condition=='CS-'][ae.memory_phase=='encoding'],x='roi',y='rsa',hue='el',row='group',kind='bar')
# sns.catplot(data=ae[ae.condition=='CS-'][ae.memory_phase=='retrieval'],x='roi',y='rsa',hue='el',row='group',kind='bar')


ae = ae.set_index(['group','memory_phase','roi','condition','el','subject'])
stats = pd.DataFrame(columns=['w','p','cles','p_fdr'],
                         index=pd.MultiIndex.from_product([groups,memcon,cons,rois],
                         names=['group','memory_phase','condition','roi']))
for group in groups:
    for memory_phase in memcon:
        for con in cons:
            for roi in rois:
                wres = pg.wilcoxon(ae.loc[(group,memory_phase,roi,con,'early'),'rsa'], ae.loc[(group,memory_phase,roi,con,'late'),'rsa'], tail='greater')
                stats.loc[(group,memory_phase,con,roi),['w','p','cles']] = wres[['W-val','p-val','CLES']].values
            stats.loc[(group,memory_phase,con),'p_fdr'] = pg.multicomp(list(stats.loc[(group,memory_phase,con),'p'].values),method='fdr_bh')[1]
stats['p_mask'] = stats.p_fdr.apply(lambda x: 0 if x >.05 else 1)
stats['cles_disp'] = stats.p_mask * stats.cles
for group in groups:
    for memory_phase in memcon:
        for con in cons:
            disp = stats.loc[group,memory_phase,con].copy()
            if disp.cles_disp.max() == 0: 
                pass
            else:
                bnsurf(data=disp,val='cles_disp',min_val=.5,max_val=stats.cles_disp.max(),cmap='Oranges',out='replication/%s_%s_%s'%(group,memory_phase,con))

#doing the CS+ vs. CS- comparison
ae = ae.reset_index().set_index(['el','group','memory_phase','roi','condition','subject'])
ae = (ae.loc['early'] - ae.loc['late'])

stats = pd.DataFrame(columns=['w','p','cles','p_fdr'],
                         index=pd.MultiIndex.from_product([groups,memcon,rois],
                         names=['group','memory_phase','roi']))

for group in groups:
    for memory_phase in memcon:
        for roi in rois:
            wres = pg.wilcoxon(ae.loc[(group,memory_phase,roi,'CS+'),'rsa'], ae.loc[(group,memory_phase,roi,'CS-'),'rsa'], tail='greater')
            stats.loc[(group,memory_phase,roi),['w','p','cles']] = wres[['W-val','p-val','CLES']].values
        stats.loc[(group,memory_phase,bn_rois),'p_fdr'] = pg.multicomp(list(stats.loc[(group,memory_phase,bn_rois),'p'].values),method='fdr_bh')[1]
stats['p_mask'] = stats.p_fdr.apply(lambda x: 0 if x >.05 or pd.isna(x) else 1)
stats['cles_disp'] = stats.p_mask * stats.cles



#just one roi for vis. purposes
a = ae[ae.group == 'healthy'][ae.roi == 'A9m'][ae.condition == 'CS+']
a = a.rename(columns={'el':'Extinction half'})
fig, ax = plt.subplots(figsize=(8,5))
sns.barplot(data=a,x='memory_phase',y='rsa',hue='Extinction half',ax=ax)
ax.set_ylabel('Similarity to late acquisition')
ax.set_xlabel('Memory phase')
ax.set_title('ROI = A9m')

#########################Early extinction to late extinction vs. late acquisition#######################
cae = pd.DataFrame(index=pd.MultiIndex.from_product([memcon,['acquisition','extinction'],cons,rois,sub_args],names=['memory_phase','phase','condition','roi','subject']))
pae = pd.DataFrame(index=pd.MultiIndex.from_product([memcon,['acquisition','extinction'],cons,rois,p_sub_args],names=['memory_phase','phase','condition','roi','subject']))
for mem in memcon:
    for phase in ['acquisition','extinction']:
        for con in cons:
            for roi in rois:
                print(roi)
                for sub in subs:                                        
                    cae.loc[(mem,phase,con,roi,sub_args[sub]),'rsa'] = get_square(c.mats,roi,sub_args[sub],sub,split_slices,con,'early_extinction',mem,con,'late_%s'%(phase),mem)
                    pae.loc[(mem,phase,con,roi,p_sub_args[sub]),'rsa'] = get_square(p.mats,roi,p_sub_args[sub],sub,split_slices,con,'early_extinction',mem,con,'late_%s'%(phase),mem)
ae = pd.concat((cae,pae)).reset_index()
ae['group'] = ae.subject.apply(lgroup)

ae = ae.set_index('subject').drop([20,120]).reset_index()
sns.catplot(data=ae[ae.condition=='CS+'][ae.memory_phase=='encoding'],x='roi',y='rsa',hue='phase',row='group',kind='bar')
sns.catplot(data=ae[ae.condition=='CS+'][ae.memory_phase=='retrieval'],x='roi',y='rsa',hue='phase',row='group',kind='bar')










#############WITHIN PHASE SIMILARITY########################
phase2 = ['acquisition','extinction']
c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)

cdf = pd.DataFrame(index=pd.MultiIndex.from_product([phase2,memcon,cons,bn_rois,sub_args],names=['phase','memory_phase','condition','roi','subject']))
pdf = pd.DataFrame(index=pd.MultiIndex.from_product([phase2,memcon,cons,bn_rois,p_sub_args],names=['phase','memory_phase','condition','roi','subject']))

for phase in phase2:
    for memory_phase in memcon:
        for con in cons:
            for roi in rois:
                for sub in subs:                                                
                    cdf.loc[(phase,memory_phase,con,roi,sub_args[sub]),'rsa'] = get_square(c.mats,roi,sub_args[sub],sub,slice3,con,phase,memory_phase,con,phase,memory_phase)
                    pdf.loc[(phase,memory_phase,con,roi,p_sub_args[sub]),'rsa'] = get_square(p.mats,roi,p_sub_args[sub],sub,slice3,con,phase,memory_phase,con,phase,memory_phase)
df = pd.concat((cdf,pdf))
# df = df.drop('baseline')
# diff = (df.loc['extinction'] - df.loc['acquisition'])

df = df.reset_index()
df = df.set_index('subject').drop([20,120]).reset_index()
df['group'] = df.subject.apply(lgroup)
df = df.set_index(['group','phase','memory_phase','condition','roi','subject'])

contrasts = ['encoding','retrieval','ER']
stats = pd.DataFrame(columns=['w','p','cles','p_fdr'],
                     index=pd.MultiIndex.from_product([groups,phase2,memcon,bn_rois],
                     names=['group','phase','memory_phase','roi']))

for group in groups:
    # for con in cons:
    for phase in phase2:
        for memory_phase in memcon:
            for roi in rois:
                wres = pg.wilcoxon(df.loc[(group,phase,memory_phase,'CS+',roi),'rsa'], df.loc[(group,phase,memory_phase,'CS-',roi),'rsa'],tail='greater')
                stats.loc[(group,phase,memory_phase,roi),['w','p','cles']] = wres[['W-val','p-val','CLES']].values
            stats.loc[(group,phase,memory_phase,bn_rois),'p_fdr'] = pg.multicomp(list(stats.loc[(group,phase,memory_phase,bn_rois),'p'].values),method='fdr_bh')[1]
stats.p = stats.p.apply(lambda x: x[0] if type(x) == list else x)
stats.cles = stats.cles.apply(lambda x: x[0] if type(x) == list else x)
stats['p_mask'] = stats.p_fdr.apply(lambda x: 0 if x >.05 or pd.isna(x) else 1)

# stats.loc[('healthy','CS+','retrieval','A14m'),'p_mask'] = 1

stats['cles_disp'] = stats.cles * stats.p_mask

stats = stats.sort_index()
cmaps = {'acquisition':'Purples',
         'extinction':'Greens'}
for group in groups:
    # for con in cons:
    for phase in phase2:
        MAX = stats.loc[(slice('healthy','ptsd'),phase),'cles_disp'].max()
        for memory_phase in memcon:
            disp = stats.loc[group,phase,memory_phase].copy()
            if disp.cles_disp.max() == 0: 
                pass
            else:
                bnsurf(data=disp,val='cles_disp',min_val=.5,max_val=MAX,cmap=cmaps[phase],out='within_phase/%s_%s_%s'%(group,phase,memory_phase))

#difference between encoding and retreival

diff = df.reset_index().set_index(['phase','group','memory_phase','condition','roi','subject'])
diff = (diff.loc['extinction'] - diff.loc['acquisition'])
sns.catplot(data=diff.loc['healthy'].reset_index(),x='roi',y='rsa',hue='memory_phase',row='condition',palette=cpal,kind='bar')
sns.catplot(data=diff.loc['ptsd'].reset_index(),x='roi',y='rsa',hue='memory_phase',row='condition',palette=cpal,kind='bar')

diff_stats = pd.DataFrame(columns=['w','p','cles','p_fdr'],
                     index=pd.MultiIndex.from_product([groups,cons,rois],
                     names=['group','condition','roi']))
for group in groups:
    for con in cons:
        for roi in rois:
            wres = pg.wilcoxon(diff.loc[(group,'retrieval',con,roi),'rsa'], diff.loc[(group,'encoding',con,roi),'rsa'],tail='two-sided')
            diff_stats.loc[(group,con,roi),['w','p','cles']] = wres[['W-val','p-val','CLES']].values
        diff_stats.loc[(group,con,bn_rois),'p_fdr'] = pg.multicomp(list(diff_stats.loc[(group,con,bn_rois),'p'].values),method='fdr_bh')[1]
diff_stats.p = diff_stats.p.apply(lambda x: x[0] if type(x) == list else x)
diff_stats.cles = diff_stats.cles.apply(lambda x: x[0] if type(x) == list else x)
diff_stats['p_mask'] = diff_stats.p_fdr.apply(lambda x: 0 if x >.05 or pd.isna(x) else 1)
 
#####################CS+E early to CS+ late encoding vs. retreival#####################
cdf = pd.DataFrame(index=pd.MultiIndex.from_product([memcon,cons,rois,sub_args],names=['memory_phase','condition','roi','subject']))
pdf = pd.DataFrame(index=pd.MultiIndex.from_product([memcon,cons,rois,p_sub_args],names=['memory_phase','condition','roi','subject']))

for memory_phase in memcon:
    for con in cons:
        for roi in rois:
            for sub in subs:                                                
                cdf.loc[(memory_phase,con,roi,sub_args[sub]),'rsa'] = get_square(c.mats,roi,sub_args[sub],sub,slice3,con,'extinction',memory_phase,con,'acquisition',memory_phase)
                pdf.loc[(memory_phase,con,roi,p_sub_args[sub]),'rsa'] = get_square(p.mats,roi,p_sub_args[sub],sub,slice3,con,'extinction',memory_phase,con,'acquisition',memory_phase)
df = pd.concat((cdf,pdf)).reset_index()
df = df.set_index('subject').drop([20,120]).reset_index()
df['group'] = df.subject.apply(lgroup)
df = df.set_index(['condition','memory_phase','group','roi','subject'])
# diff = (df.loc['retrieval'] - df.loc['encoding'])
diff = (df.loc['CS+'] - df.loc['CS-'])
# diff = df.reset_index()
# sns.catplot(data=diff[diff.group == 'healthy'],x='roi',y='rsa',hue='condition',row='memory_phase',palette=cpal,kind='bar')
df = df.reset_index().set_index(['group','condition','memory_phase','roi','subject'])


stats = pd.DataFrame(columns=['w','p','cles','p_fdr'],
                     index=pd.MultiIndex.from_product([groups,memcon,rois],
                     names=['group','memory_phase','roi']))

for group in groups:
    # for memory_phase in memcon:
    for memory_phase in memcon:
        for roi in rois:
            wres = pg.wilcoxon(df.loc[(group,'CS+',memory_phase,roi),'rsa'], df.loc[(group,'CS-',memory_phase,roi),'rsa'],tail='two-sided')
            stats.loc[(group,memory_phase,roi),['w','p','cles']] = wres[['W-val','p-val','CLES']].values
        stats.loc[(group,memory_phase,bn_rois),'p_fdr'] = pg.multicomp(list(stats.loc[(group,memory_phase,bn_rois),'p'].values),method='fdr_bh')[1]
stats.p = stats.p.apply(lambda x: x[0] if type(x) == list else x)
stats.cles = stats.cles.apply(lambda x: x[0] if type(x) == list else x)
stats['p_mask'] = stats.p_fdr.apply(lambda x: 0 if x >.05 or pd.isna(x) else 1)
stats['cles_disp'] = stats.cles * stats.p_mask


for group in groups:
    for memory_phase in memcon:
        MAX = stats.loc[(slice('healthy','ptsd'),memory_phase),'cles_disp'].max()
        disp = stats.loc[group,memory_phase].copy()
        if disp.cles_disp.max() == 0: 
            pass
        else:
            bnsurf(data=disp,val='cles_disp',min_val=.5,max_val=MAX,cmap='Oranges',out='replication/Ext_Acq_CS_comp_%s_%s'%(group,memory_phase))


####################at retrieval, where are there differences in item level reinstatement between CS+E and CS+A#####################
c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)

cdf = c.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = p.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
mdf = pd.concat((cdf,pdf)).reset_index()
# mdf = (mdf.loc['CS+'] - mdf.loc['CS-']).reset_index()

mdf = mdf.set_index('subject').drop([20,120]).reset_index()
mdf['group'] = mdf.subject.apply(lgroup)
mdf = mdf.set_index(['group','roi','encode_phase','trial_type','subject']).sort_index()

stats = pd.DataFrame(columns=['w','p','cles','p_corr'],index=pd.MultiIndex.from_product([groups,bn_rois],names=['group','roi']))
for group in ['healthy','ptsd']:
    for roi in bn_rois:
        wres = pg.wilcoxon(mdf.loc[(group,roi,'extinction','CS+'),'rsa'],mdf.loc[(group,roi,'acquisition','CS+'),'rsa'],tail='two-sided')
        stats.loc[(group,roi),['w','p','cles']] = wres[['W-val','p-val','CLES']].values
    stats.loc[(group),'p_corr'] = pg.multicomp(list(stats.loc[(group),'p'].values),method='fdr_bh')[1]
stats['p_mask'] = stats.p.apply(lambda x: 0 if x >.05 else 1)
stats['cles_disp'] = stats.cles * stats.p_mask

for group in groups:
    disp = stats.loc[group].copy()
    if disp.cles_disp.max() == 0: 
        pass
    else:
        bnsurf(data=disp,val='cles_disp',min_val=.25,max_val=.75,mid_val=.5,cmap='PRGn',out='advanced_rsa/CSpE_CSpA_%s_%s'%(group,phase))

#################memory logreg -again?
R = np.random.RandomState(42)

def boot_roi_logreg(bdf,n_boot=10000):
    logreg = LogisticRegression(solver='lbfgs')

    boot_res = np.zeros(n_boot)
    
    subs = bdf.reset_index()['subject'].unique()
    n_subs = subs.shape[0]

    for i in range(n_boot):
        y = np.zeros(n_subs * 24)
    
        while len(np.unique(y)) == 1:
            _samp = R.choice(subs,n_subs)
            X = bdf.loc[_samp,'rsa'].values.reshape(-1,1)
            y = bdf.loc[_samp,'high_confidence_accuracy'].values
        
        logreg.fit(X,y)
        boot_res[i] = logreg.coef_

    avg = boot_res.mean()
    pval = 1 - np.mean(boot_res > 0)
    CI = [np.percentile(boot_res,5),np.percentile(boot_res,100)] 

    return avg, pval, CI

df = pd.concat((c.df,p.df)).set_index('subject').drop([18,20,120]).reset_index()
df['group'] = df.subject.apply(lgroup)
df = df.set_index(['roi','group','encode_phase','trial_type','subject']).sort_index()
df = df.loc[bn_rois]

stats = pd.DataFrame(columns=['avg','p','CI'],index=pd.MultiIndex.from_product([['healthy','ptsd'],phase3,cons,bn_rois],names=['group','encode_phase','condition','roi']))

for group in groups:
    print(group)
    for phase in phase3:
        print(phase)
        for con in cons:
            print(con)
            for roi in bn_rois: 
                stats.loc[(group,phase,con,roi),['avg','p','CI']] = boot_roi_logreg(df.loc[(roi,group,phase,con)].copy())

self.betas = betas


##########subcortical main effect emotional ERS#############
c = group_roi_rsa(group='control',ext_split=False,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=False,fs=True,hemi=False)

cdf = c.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = p.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
mdf = pd.concat((cdf,pdf)).reset_index()
# mdf = (mdf.loc['CS+'] - mdf.loc['CS-']).reset_index()

mdf = mdf.set_index('subject').drop([20,120]).reset_index()
mdf['group'] = mdf.subject.apply(lgroup)
mdf = mdf.set_index(['group','roi','encode_phase','trial_type','subject']).sort_index()
