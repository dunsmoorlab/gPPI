from roi_rsa import *
c = group_roi_rsa(group='control',ext_split=True,fs=True,hemi=False)
p = group_roi_rsa(group='ptsd',ext_split=True,fs=True,hemi=False)


cons = ['CS+','CS-']
rois = ['mOFC', 'dACC', 'amyg', 'hpc', 'ins']
phases = {'baseline':24,'acquisition':24,'early_extinction':8,'extinction':16}
subs = range(24)

from graphing_functions import *
##############Item level##############################
cdf = c.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
cdf = (cdf.loc['CS+'] - cdf.loc['CS-']).reset_index()
cdf['group'] = cdf.subject.apply(lgroup)

pdf = p.df.groupby(['trial_type','encode_phase','roi','subject']).mean()
pdf = (pdf.loc['CS+'] - pdf.loc['CS-']).reset_index()

pdf['group'] = pdf.subject.apply(lgroup)

pfc = pd.concat((cdf,pdf))
pfc.roi = pfc.roi.apply(lambda x: 'vmPFC' if x == 'mOFC' else x)

pfc = pfc.set_index(['group','roi','encode_phase'])
cscomp('control',pfc,['vmPFC','dACC'],phases=phases.keys())
cscomp('ptsd',pfc,['vmPFC','dACC'],phases=phases.keys())

##############Set level##############################

slices={'CS+':{
                 'baseline':{'encoding':slice(0,24),
                             'retrieval':slice(144,168)},
              
              'acquisition':{'encoding':slice(24,48),
                             'retrieval':slice(168,192)},

         'early_extinction':{'encoding':slice(48,56),
                            'retrieval':slice(192,200)},
               
               'extinction':{'encoding':slice(56,72),
                             'retrieval':slice(200,216)}},
        'CS-':{
                 'baseline':{'encoding':slice(72,96),
                            'retrieval':slice(216,240)},
              
              'acquisition':{'encoding':slice(96,120),
                            'retrieval':slice(240,264)},

         'early_extinction':{'encoding':slice(120,128),
                            'retrieval':slice(264,272)},
               
               'extinction':{'encoding':slice(128,144),
                            'retrieval':slice(272,288)}}}



csl = pd.DataFrame(index=pd.MultiIndex.from_product([cons,rois,phases,subs],names=['condition','roi','encode_phase','subject']))
psl = pd.DataFrame(index=pd.MultiIndex.from_product([cons,rois,phases,subs],names=['condition','roi','encode_phase','subject']))
for con in cons:
    for roi in rois:
        for phase in phases:
            for sub in subs:
                csl.loc[(con,roi,phase,sub),'rsa'] = c.mats[roi][sub,slices[con][phase]['encoding'],slices[con][phase]['retrieval']][np.eye(phases[phase]) == 1].mean()
                psl.loc[(con,roi,phase,sub),'rsa'] = p.mats[roi][sub,slices[con][phase]['encoding'],slices[con][phase]['retrieval']][np.eye(phases[phase]) == 1].mean()

csl = (csl.loc['CS+'] - csl.loc['CS-']).reset_index()
psl = (psl.loc['CS+'] - psl.loc['CS-']).reset_index()
csl.subject = csl.subject.astype(int)
psl.subject = psl.subject.astype(int) + 101
sl = pd.concat([csl,psl])
sl['group'] = sl.subject.apply(lgroup)

sl = sl.set_index(['group','roi','encode_phase'])
cscomp('control',sl,['mOFC','dACC'],phases=phases.keys())
cscomp('ptsd',sl,['mOFC','dACC'],phases=phases.keys())
