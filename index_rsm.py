from fg_config import *
ref = pd.read_csv('na_reference.csv').set_index(['subject','trial_type','encode_phase']).sort_index()
con_convert = {'CS+':'CSp_trial',
               'CS-':'CSm_trial'}
phase_convert = {'baseline':'baseline',
                'acquisition':'acquisition',
                'extinction':'extinction',
                'early_acquisition':'acquisition',
                'late_acquisition':'acquisition',
                'early_extinction':'extinction',
                'late_extinction':'extinction'}


#pythonic indexing
ref.CSp_trial = ref.CSp_trial - 1
ref.CSm_trial = ref.CSm_trial - 1

def get_square(mats,roi,subject,sindex,slices,condition1,phase1,memcon1,condition2,phase2,memcon2):
    
    square = mats[roi][sindex, 

            slices[condition1][phase1][memcon1],

            slices[condition2][phase2][memcon2]]

    if memcon1 == 'retrieval':

        dat = ref.loc[subject,condition1,phase_convert[phase1]].copy()

        if dat.response.isna().sum() == 0: pass
        else:
            to_drop = dat[con_convert[condition1]][dat.response.isna()].astype(int).values
            square = np.delete(square,to_drop,axis=0)

    if memcon2 == 'retrieval':

        dat = ref.loc[subject,condition2,phase_convert[phase2]].copy()

        if dat.response.isna().sum() == 0: pass
        else:
            to_drop = dat[con_convert[condition2]][dat.response.isna()].astype(int).values
            square = np.delete(square,to_drop,axis=1)

    if memcon1 == memcon2 and phase1 == phase2 and condition1 == condition2:
        #if we are along the diagonal, just get the lower triangle
        return square[np.tril_indices(square.shape[0],k=-1)].mean()
    
    else:
        return square.mean()