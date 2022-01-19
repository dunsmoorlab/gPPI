from fg_config import *

fa = {sub: bids_meta(sub).mem_df for sub in all_sub_args}
fa = pd.concat(fa.values())
fa = fa[fa.memory_condition == 'New'].reset_index(drop=True)
fa['fa'] = fa.high_confidence_accuracy.apply(lambda x: 1 if x == 'FA' else 0)
fa['group'] = fa.subject.apply(lgroup)
fa = fa[['group','subject','trial_type','stimulus','fa']]
fa = fa.rename(columns={'trial_type':'condition'})
fa = fa.groupby(['subject','condition','group']).mean()
fa.to_csv('false_alarm_sub_means.csv')
# sns.pointplot(data=fa.reset_index(),x='group',y='fa',hue='condition')

'''just getting high confidence hit rate'''
hits = {sub: bids_meta(sub).mem_df for sub in all_sub_args}
hits = pd.concat(hits.values())
hits = hits[hits.memory_condition == 'Old'].reset_index(drop=True)
hits['hit'] = hits.high_confidence_accuracy.apply(lambda x: 1 if x == 'H' else 0)
hits['group'] = hits.subject.apply(lgroup)
hits = hits[['subject','group','encode_phase','trial_type','stimulus','hit']]
hits = hits.rename(columns={'trial_type':'condition','encode_phase':'phase'})
hits = hits.groupby(['subject','condition','group','phase']).mean()
hits['fa'] = fa.fa
hits['cr'] = hits.hit - corrected_recognition_sub_meanshits.fa
hits.to_csv('.csv')

'''doing hit rate by phase/group and confidence'''
df = {sub: bids_meta(sub).mem_df for sub in all_sub_args}
df = pd.concat(df.values())
df = df[df.memory_condition == 'Old'].reset_index(drop=True)
df.low_confidence_accuracy = df.low_confidence_accuracy.apply(lambda x: 1 if x == 'H' else 0) 
df.high_confidence_accuracy = df.high_confidence_accuracy.apply(lambda x: 1 if x == 'H' else 0)
df = df.rename(columns={'trial_type':'condition','encode_phase':'phase'})
df = df[['subject','phase','condition','low_confidence_accuracy','high_confidence_accuracy']]
df = df.groupby(['subject','phase','condition']).mean().reset_index()
df['group'] = df.subject.apply(lgroup)
acc = df.drop(columns='subject').groupby(['group','phase','condition']).mean()
acc = acc.rename(columns={'low_confidence_accuracy':'low_hit','high_confidence_accuracy':'high_hit'})
err = df.drop(columns='subject').groupby(['group','phase','condition']).std(ddof=1)
err = err.rename(columns={'low_confidence_accuracy':'low_std','high_confidence_accuracy':'high_std'})
df = pd.concat((acc,err),axis=1)
df = df.reindex(axis='index', level=1, labels=['baseline','acquisition','extinction'])

def f(x,y):
    return f'{round(x,2)} ({round(y,2)})'
df['low'] = df.apply(lambda x: f(x.low_hit, x.low_std), axis=1)
df['high'] = df.apply(lambda x: f(x.high_hit, x.high_std), axis=1)
df = df[['low','high']]
df.to_csv('mem_data_for_paper.csv')


'''random rsa graphing'''
# df = pd.read_csv('cb_response_rsa.csv').groupby(['condition','group','subject','phase','roi']).mean()#.reset_index()
# df = (df.loc['CS+'] - df.loc['CS-']).reset_index()
# df = df[df.phase.isin(['conditioning','extinction'])]
# rois = ['vmPFC','dACC']
# sns.catplot(data=df[subdf.roi.isin(rois)],x='phase',y='rand_rsa',hue='roi',row='group',kind='bar')

df = pd.read_csv('cb_response_rsa.csv')
df['item vs. off diagonal'] = df.ers - df.off_ers
df = df[df.phase.isin(['conditioning','extinction'])]
df = df.groupby(['group','subject','roi']).mean().reset_index()
df.group = df.group.apply(lambda x: 'Healthy' if x == 'healthy' else 'PTSS')
o = df.groupby(['group','roi']).mean().sort_values('item vs. off diagonal',ascending=False).reset_index(level=-1).loc['Healthy','roi']
g = sns.catplot(data=df,x='roi',y='item vs. off diagonal',row='group',kind='bar',order=o.values,palette='crest')
g = g.set_xticklabels(['fusiform','dACC','BLA','Precuneus','vmPFC','HC Body','Ant. Insula','CeM','pHC','aHC']