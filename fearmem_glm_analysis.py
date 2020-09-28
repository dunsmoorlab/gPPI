from fg_config import *
def convert_to_arr(x):
    for j in ['\n','[',']']:
        x = x.replace(j,'')
    x = [i for i in x.split(' ') if i != '']
    return np.array(x).astype(float)
effects = ['Condition','Response','Condition_Response'] 
phases = ['baseline','acquisition','extinction']

eff = 'Response'

df = pd.read_csv(f'sm_events/{eff}_extracted_pe.csv')
df.condition = df.condition.apply(lambda x: 'CS+' if x == 'CSp' else 'CS-')
df = df.groupby(['roi','encode_phase','condition','source_memory','subject']).mean().sort_index()

resp_eff = df.groupby(['roi','source_memory','subject']).mean()
q = sns.catplot(data=resp_eff.reset_index(),x='source_memory',y='pe',col='roi',order=phases,kind='bar',palette=spal)


roi = 'cuneus'

roidf = df.loc[roi].copy()



fig, ax = plt.subplots(1,3,figsize=(12,5),sharey=True)
for i, phase in enumerate(phases):
    dat = roidf.loc[phase].copy().reset_index()
    sns.barplot(data=dat,x='condition',y='pe',hue='source_memory',ax=ax[i],
        palette=spal)
    if i != 0:ax[i].set_ylabel('')
    else:ax[i].set_ylabel('Beta')
    ax[i].legend_.remove()
# legend_elements = [Patch(facecolor=spal[0],edgecolor=None,label='Baseline'),
#                    Patch(facecolor=spal[1],edgecolor=None,label='Acquisition'),
#                    Patch(facecolor=spal[2],edgecolor=None,label='Extinction')]
# ax[0].legend(handles=legend_elements,loc='upper right')
# ax[0].legend_.set_title('Source memory\n      response')
# ax[0].set_xlabel('Baseline')
# ax[1].set_xlabel('Acquisition')
# ax[2].set_xlabel('Extinction')
