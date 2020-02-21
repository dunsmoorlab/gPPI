import numpy as np
import nibabel as nib
from scipy.stats import pearsonr, zscore
from scipy.signal import convolve
from brainiak_utils import _double_gamma_hrf, convolve_hrf
from sklearn.linear_model import Ridge, RidgeCV
from scipy.fftpack import dct
from scipy.linalg import toeplitz

data = nib.load('../../nistats-bids/sub-FC001/ses-1/func/sub-FC001_ses-1_task-acquisition_bold.nii.gz')
data = data.get_data()

bold = data[35,10,14,:]
# bold = bold - bold.mean()
x = range(bold.shape[0])
plt.plot(x,bold)

#####
dt = 2 #microscan time?
N = 259 #number of timepoints
TR = 2 #length of TR in seconds
NT = TR/dt#NUMBER OF BINS PER TR - if this is set to 1 then I think it models just the midpoint of the TR, which is what we want for slice time corrected data
k = np.arange(0,N*NT,NT,dtype=int) #holdover from matlab code
###in python this code does:
hrf = np.array(_double_gamma_hrf(temporal_resolution=NT/TR))
# plt.subplots(); plt.plot(hrf)
def dct_mat(N_,K_):
    n = np.array((range(0,N_))).T
    C_ = np.zeros((n.shape[0],K_))
    C_[:,0] = np.ones(n.shape[0])/np.sqrt(N_)
    for q in range(1,K_):
        C_[:,q] = np.sqrt(2/N_)* np.cos( np.pi*(2*n)* (q) /(2*N_)) #what is this???/?????
    return C_
xb = dct_mat(int(N*NT),N)

# def dft_mat(N_,K_):
#     n = np.array((range(0,N_))).T
#     C_ = np.zeros((n.shape[0],K_))
#     for k in range(K_):
#         y =  np.cos(np.pi*k*(2*n+1)/(2*N))
#         if k == 0: C_[:,k] = y * np.sqrt(1/(4*N_))
#         else:      C_[:,k] = y * np.sqrt(1/(2*N_))
#     return C_
# xb = dft_mat(int(N*NT),N)

Hxb = np.zeros((N,N))
for i in range(N):
    Hx = convolve(xb[:,i],hrf)
    Hxb[:,i] = Hx[k]
xb = xb[:,:]

####ethan help code
reg = Ridge(alpha=.001,solver='lsqr',fit_intercept=True,normalize=False,max_iter=1000)
reg.fit(Hxb,bold)
neuronal = np.matmul(xb,reg.coef_)
plt.subplots();plt.plot(x,neuronal)
plt.legend(['Neuronal signal'])
# recon = convolve(neuronal,hrf)[:bold.shape[0]] + reg.intercept_
recon = convolve(neuronal , hrf)[:bold.shape[0]] + reg.intercept_

plt.subplots()
plt.plot(x,bold)
# plt.plot(x,recon)
plt.plot(x,recon)
plt.title(pearsonr(bold,recon))
plt.legend(['Original BOLD','Reconstructed BOLD'])

alphas = [.001,.1,1,10,1000]
for a in alphas:
    reg = Ridge(alpha=a,solver='lsqr',fit_intercept=True,normalize=False,max_iter=1000)
    reg.fit(Hxb,bold)
    neuronal = np.matmul(xb,reg.coef_)
    plt.plot(x,neuronal)
plt.legend(alphas)


labels = pd.read_csv('../../nistats-bids/sub-FC001/ses-1/func/sub-FC001_ses-1_task-acquisition_events.tsv',sep='\t')
labels = labels[labels.trial_type == 'CS+'][['onset','duration']].round()
design = np.zeros(N)
for i in labels.index:
    onset = int(labels.loc[i,'onset'])
    duration = int(labels.loc[i,'duration'])
    design[onset:onset+duration] = 1

lazy = convolve(design,hrf)[:bold.shape[0]] * (bold - bold.mean())
better = convolve((neuronal-neuronal.mean())*design,hrf)[:bold.shape[0]]
plt.subplots();plt.plot(x,design)
plt.subplots()
plt.plot(x,lazy)
plt.plot(x,better)
plt.legend(['BOLD * HRF(Task)','HRF(Neuronal * task)'])