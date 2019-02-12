i=1;
freq=1000;
TR=2;

temp=sprintf('C:\\Users\\hatzi\\Desktop\\fMRI data\\ds116\\sub-0%d\\func\\task001_run01\\task001_run01.feat\\filtered_func_data.nii.gz',i);
tempi=load_untouch_nii(temp);
tempimg=reshape(tempi.img,[],size(tempi.img,4));

temp=sprintf('C:\\Users\\hatzi\\Desktop\\fMRI data\\data eeg fmri\\sub0%d\\model\\model001\\onsets\\task001_run001\\cond003.txt',i);
a1=importdata(temp);

temp=sprintf('C:\\Users\\hatzi\\Desktop\\fMRI data\\data eeg fmri\\sub0%d\\model\\model001\\onsets\\task001_run001\\cond002.txt',i);
a2=importdata(temp);

temp=sprintf('C:\\Users\\hatzi\\Desktop\\fMRI data\\data eeg fmri\\sub0%d\\EEG\\task001_run001\\EEG_rereferenced.mat',i);
load(temp);

sizet=size(tempimg,2);
toepl_hrf=toeplitz([hrf(1) zeros(1,sizet-1)],[hrf zeros(1,sizet-1)]);
Nt=TR*freq;
%%Downsampling matrx

d=zeros(1,size(tempimg,2)*Nt);
d(1)=1;
for i=1:size(tempimg,2)/Nt
    downs(i,:)=circshift(d,[1,Nt*(i-1)]);
end
    
%%downsampl+hrf matrix
X=toepl_hrf'*downs;


