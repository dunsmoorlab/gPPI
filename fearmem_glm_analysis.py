from fg_config import *
effects = ['Condition','Response','Condition_Response'] 
phases = ['baseline','acquisition','extinction']

df = pd.read_csv('sm_events/extracted_pe.csv'
    ).dropna(subset=['pe']
    )
def convert_to_arr(x):
    for j in ['\n','[',']']:
        x = x.replace(j,'')
    x = [i for i in x.split(' ') if i != '']
    return np.array(x).astype(float)

df.pe = df.pe.apply(convert_to_arr)
df = pd.concat((df,df.pe.apply(pd.Series)),axis=1).drop(columns='pe')
df.condition = df.condition.apply(lambda x: 'CS+' if x == 'CSp' else 'CS-')
df = df.groupby(['effect','subject','encode_phase','condition','source_memory']).mean().drop(columns='run').sort_index()
df.to_csv('sm_events/cleaned_pe_estimates.csv')
'''now actually analyze'''

indf = pd.read_csv('sm_events/cleaned_pe_estimates.csv').set_index(['effect','subject','encode_phase','condition','source_memory'])

eff = 'Response'
cluster = '0'

df = indf.loc[eff,cluster].reset_index().rename(columns={cluster:'pe'})
fig, ax = plt.subplots(1,3,figsize=(12,5),sharey=True)
for i, phase in enumerate(phases):
    dat = df[df.encode_phase == phase].copy()
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
