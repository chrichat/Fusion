for i=1:14
    if (i~=2)
        Fz_o(i,:)=squeeze(mean(TEEG(7,:,(2:2:100),i),3));
        Fz_n(i,:)=squeeze(mean(TEEG(7,:,(1:3:100),i),3));
        Pz_o(i,:)=squeeze(mean(TEEG(25,:,(2:2:100),i),3));
        Pz_n(i,:)=squeeze(mean(TEEG(25,:,(1:2:100),i),3));
     end
end


s=std(Fz_o);
s2=std(Fz_n);


mfz_o=mean(Fz_o);
mfz_n=mean(Fz_n);
y=1:700;
figure;plot(y,mfz_o)
% errorbar(y,mfz_o,s)
hold on
plot(mfz_n,'r');
% errorbar(y,mfz_n,s2)
hline=refline(0,0);
hline.Color = [0.1 0.1 0.1];

ss=std(Pz_o);
ss2=std(Pz_n);
figure
mpz_o=mean(Pz_o);
mpz_n=mean(Pz_n);
figure;plot(circshift(mpz_o/max(abs(mpz_o)),[0,30]))
% errorbar(y,mpz_o,ss)
% hold on
% plot(mpz_n/max(mpz_n),'r');
% errorbar(y,mpz_n,ss2)
hold on
plot(circshift(solra.factors.B(:,7)/max(abs(solra.factors.B(:,7))),[30,0]),'r')
hline=refline(0,0);
hline.Color = [0.1 0.1 0.1];
legend('Oddball ERP extracted from Mean Analysis', 'Oddball ERP extracted from Soft Coupled CPDs')

for i=1:14
    temp=sprintf('C:\\Users\\hatzi\\Desktop\\fMRI data\\data eeg fmri\\sub0%d\\EEG\\task001_run001\\EEG_rereferenced.mat',i);
    load(temp); 
    temp=sprintf('C:\\Users\\hatzi\\Desktop\\fMRI data\\data eeg fmri\\sub0%d\\model\\model001\\onsets\\task001_run001\\cond003.txt',i);
    a1=importdata(temp);
    sizea1=size(a1,1);
    for j=1:sizea1
        oddmatrix(j,:,:)=data_reref(:,round(a1(j,1)*1000):round(a1(j,1)*1000+700));
    end
    temp=sprintf('C:\\Users\\hatzi\\Desktop\\fMRI data\\data eeg fmri\\sub0%d\\model\\model001\\onsets\\task001_run001\\cond002.txt',i);
    b1=importdata(temp);
    sizeb1=size(b1,1);
    for j=1:sizeb1
        normatrix(j,:,:)=data_reref(:,round(b1(j,1)*1000):round(b1(j,1)*1000+700));
    end      
    temp=sprintf('C:\\Users\\hatzi\\Desktop\\fMRI data\\data eeg fmri\\sub0%d\\EEG\\task001_run002\\EEG_rereferenced.mat',i);
    load(temp);
    temp=sprintf('C:\\Users\\hatzi\\Desktop\\fMRI data\\data eeg fmri\\sub0%d\\model\\model001\\onsets\\task001_run002\\cond003.txt',i);
    a2=importdata(temp);
    sizea2=size(a2,1);
    for j=1:sizea2
        oddmatrix(sizea1+j,:,:)=data_reref(:,round(a2(j,1)*1000):round(a2(j,1)*1000+700));
    end
    temp=sprintf('C:\\Users\\hatzi\\Desktop\\fMRI data\\data eeg fmri\\sub0%d\\model\\model001\\onsets\\task001_run002\\cond002.txt',i);
    b2=importdata(temp);
    sizeb2=size(b2,1);
    for j=1:sizeb2
        normatrix(sizeb1+j,:,:)=data_reref(:,round(b2(j,1)*1000):round(b2(j,1)*1000+700));
    end      
    temp=sprintf('C:\\Users\\hatzi\\Desktop\\fMRI data\\data eeg fmri\\sub0%d\\EEG\\task001_run003\\EEG_rereferenced.mat',i);
    load(temp);
    temp=sprintf('C:\\Users\\hatzi\\Desktop\\fMRI data\\data eeg fmri\\sub0%d\\model\\model001\\onsets\\task001_run003\\cond003.txt',i);
    a3=importdata(temp);
    sizea3=size(a3,1);
    for j=1:sizea3
        oddmatrix(sizea1+sizea2+j,:,:)=data_reref(:,round(a3(j,1)*1000):round(a3(j,1)*1000+700));
    end
    temp=sprintf('C:\\Users\\hatzi\\Desktop\\fMRI data\\data eeg fmri\\sub0%d\\model\\model001\\onsets\\task001_run003\\cond002.txt',i);
    b3=importdata(temp);
    sizeb3=size(b3,1);
    for j=1:sizeb3
        normatrix(sizeb1+sizeb2+j,:,:)=data_reref(:,round(b3(j,1)*1000):round(b3(j,1)*1000+700));
    end  
    aFz_o(i,:)=squeeze(mean(oddmatrix(:,7,:),1));
    aFz_n(i,:)=squeeze(mean(normatrix(:,7,:),1));
    aPz_o(i,:)=squeeze(mean(oddmatrix(:,25,:),1));
    aPz_n(i,:)=squeeze(mean(normatrix(:,25,:),1));
end

s=std(aFz_o);
s2=std(aFz_n);


amfz_o=mean(aFz_o);
amfz_n=mean(aFz_n);
figure;plot(amfz_o)
% errorbar(y,mfz_o,s)
hold on
plot(amfz_n,'r');
% errorbar(y,mfz_n,s2)
hline=refline(0,0);
hline.Color = [0.1 0.1 0.1];
axis tight

ss=std(aPz_o);
ss2=std(aPz_n);

ampz_o=mean(aPz_o);
ampz_n=mean(aPz_n);
figure;plot(ampz_o)
% errorbar(y,mpz_o,ss)
hold on
plot(ampz_n,'r');
% errorbar(y,mpz_n,ss2)
hline=refline(0,0);
hline.Color = [0.1 0.1 0.1];
axis tight

plotme(aPz_o)

%%%TOPOPLOTS%%
for i=1:size(aPz_o,1)
    [x maxo]=max(mpz_o);
    ap300(i,:)=(mean(squeeze(oddmatrix(:,1:34,maxo)),1));
end
topoplotme(ap300)

for i=1:size(Pz_o,1)
    [x maxo]=max(mpz_o);
    p300(i,:)=(mean(squeeze(TEEG(1:34,maxo,(2:2:102),i)),2));
end
topoplotme(p300)