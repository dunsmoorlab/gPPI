
%You  need to add spm5 to the path (I think)
addpath /space/raid/fmri/spm5

%First you need to tell Matlab where to look for the ppi stuff.
%There's a copy of spm_PEB.m in here but I don't think it matters if it
%uses the spm5 version or this one.  I checked and they seem mostly identical

addpath ~/MatlabCode/fslppi_JM

% The path to the feat directory
 featdir='/space/raid6/data/engel/rad/Analysis/Intact_Scrambled/2003/2003_S1.feat';

%the path to the seed voxel.  You'll need to make this file on your
% own, probably using fslmeants with your ROI mask for the seed.

voi='/space/raid6/data/engel/rad/Analysis/Intact_Scrambled/2003/2003_S1.feat/seed.txt';

%the call works as follows: featdir, voi, contrast that picks out the
% regressor you're interetested in and the last number is a flag where
% 1 plots the output and 0 doesn't.  I've set it at 1 so you can see
% what the code makes


%I need to run it twice, once for each regressor, since I'm interested
%in the [1 -1] PPI interaction, I must run [1 0] and [0 1] separately.

[A1 B1]=fsl_ppi(featdir, voi, [1 0], 1);

[A2 B2]=fsl_ppi(featdir, voi, [0 1], 1);

%Nothing has been saved out, but you'll want to save out the
% regressors that will be fed into your PPI feat

%Your feat model will contain everything you had in it originally, but
% now you'll add the seed (you already saved that) and your 2 PPI
% regressors


%These are already convolved, so use 1 column format with *NO* HRF
%convolution in FSL.  
ppi_01=A1.ppi;
ppi_02=A2.ppi;

%you might want to use a better path, this will just stick the
%regressor in the feat directory

save('ppi_01.txt', 'ppi_01', '-ascii')
save('ppi_02.txt', 'ppi_02', '-ascii')
