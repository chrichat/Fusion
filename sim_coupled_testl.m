R1=6;
R2=6;
B=1; %For exactly the same subject factors it can be manually set or random analogously with next line
Nt=4; %Downsampling factor
% B=randn(1,1); %Keep the two matrixes analogous 


[tfmri,teeg]=Simulate_time(1,0); % The one is for the same HRF everywhere so the asumption of the model will be correct 
%I do not think that you need something extra for the time beeing, with 0
%then we have a different HRF per source. The second 0 is if we do not want ERPs
%(trial -mode in your case is not needed either).
[sfmri,seeg]=Simulate_space();

load('data/signal_hrf.mat')
hrf=sig_hrf(:,1)';

A1=abs(randn(6,6));
A2=A1;

i=6;
snr1=[0.5 0.833 1 1.4 2 4];
j=4;
snr2=[2 5 10 12 14 306];
teeg(958:end,:)=0;
Xeeg=seeg([1,2,3,4,5,6],:)'*kr(A1,teeg(:,[1,2,3,4,5,6])+randn(size(teeg))*1e-5)'; % Source 5 is not selected - blind for EEG - unshare
Xeeg=Xeeg + mean(std(Xeeg))*sqrt(snr1(i))*randn(size(Xeeg)); % 
% Xeeg2=awgn(Xeeg,1,'measured');
ww=jader(Xeeg,6);
topoplotme(pinv(ww)')
temp=abs((ww*Xeeg));
temp=reshape(temp,size(temp,1),[],size(A1,2));
plotme(abs(mean(temp,3)))
cor_ica_eegtemp(i,:)=max(abs(corr(teeg,mean(temp,3)')),[],2);
cor_ica_eegspace(i,:)=max(abs(corr(seeg',pinv(ww))),[],2);
Teeg=reshape(Xeeg,size(Xeeg,1),size(teeg,1),[]);

Xfmri=sfmri([1,2,3,4,5,6],:)'*kr(A2,tfmri([1,2,3,4,5,6],:)')'; % Source 5 is not selected - blind for EEG - unshare
% Xfmri=awgn(Xfmri,0.1,'measured');
Xfmri=Xfmri + mean(std(Xfmri))*sqrt(snr2(j))*randn(size(Xfmri)); % 
Tfmri=reshape(Xfmri,size(Xfmri,1),size(tfmri,2),[]);
ww=jader(Xfmri',6);
plotme(ww*Xfmri',70,70)
temp=reshape(pinv(ww)',size(ww,1),[],size(A1,2));
plotme(abs(mean(temp,3)))
cor_ica_fmrispac(i,:)=max(abs(corr(tfmri',mean(temp,3)')),[],2);
cor_ica_fmritemp(i,:)=max(abs(corr(sfmri',(ww*Xfmri')')),[],2);

% [U1,output] = cpd_gevd(Teeg,R1);
% [U2,output] = cpd_gevd(Tfmri,R2);
Ueeg{1,2}=teeg;
Ueeg{1,1}=seeg';
Ufmri{1,2}=tfmri';
Ufmri{1,1}=sfmri';

Teeg=Teeg./(norm(Teeg(:)));
Tfmri=Tfmri./(norm(Tfmri(:)));

[U1,~] = cpd_gevd(Teeg,R1);
[U2,~] = cpd_gevd(Tfmri,R2);

%Toeplitz matrix created from hrf 
sizet=size(Teeg,2)/Nt;
toepl_hrf=toeplitz([hrf(1) zeros(1,sizet-1)],[hrf zeros(1,sizet-1)]);

%%Downsampling matrx
d=eye(size(Teeg,2));
downs=d(1:Nt:end,:);

%%downsampl+hrf matrix
X=toepl_hrf'*downs;

%%FIX permutatin and
%%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
tempcorr=abs(corr(X*U1{1,2},U2{1,2}));
[scorr,ind]=max(tempcorr,[],2);
[~,id_un] = unique(ind,'rows','stable'); % Find unique values
dup_id=setdiff(1:size(ind,1),id_un); %Find replicated voxels and then replace them
% Check if some of the correlations are with the same map
while size(dup_id)~=0
    inds=find(ind==ind(dup_id(1)));
    [~,tempid]=min(scorr(inds));
    tempcorr(inds(tempid),ind(inds(tempid)))=0;
    [scorr,ind]=max(tempcorr,[],2);
    [~,id_un] = unique(ind,'rows','stable'); % Find unique values
    dup_id=setdiff(1:size(ind,1),id_un); %Find replicated voxels and then replace them
end
U1{1,1}(:,ind)=U1{1,1};
U1{1,2}(:,ind)=U1{1,2};
U1{1,3}(:,ind)=U1{1,3};
%%%

model1=struct;
model1.variables.a=U1{1};
model1.variables.b=U1{2};
model1.variables.c=U1{3};
model1.factors.A='a';
model1.factors.B='b';
model1.factors.C='c';
model1.factorizations.T.data=Teeg;
model1.factorizations.T.cpd={'A','B','C'};
[sol1,output]=sdf_nls(model1,'Display',100,'CGMaxIter',500,'TolFun', eps^2, 'TolX', eps, 'MaxIter', 300);
U11{1,1}=sol1.factors.A;
U11{1,2}=sol1.factors.B;
[ corune,corune_t,corune_s] = mean_corr( Ueeg,U11 );

model2=struct;
model2.variables.d=U2{1};
model2.variables.e=U2{2};
model2.variables.f=U2{3};
model2.factors.D='d';
model2.factors.E='e';
model2.factors.F='f';
model2.factorizations.T.data=Tfmri;
model2.factorizations.T.cpd={'D','E','F'};
[sol2,output]=sdf_nls(model2,'Display',100,'CGMaxIter',500,'TolFun', eps^2, 'TolX', eps, 'MaxIter', 300);
U22{1,1}=sol2.factors.D;
U22{1,2}=sol2.factors.E;
[ corunf,corunf_t,corunf_s] = mean_corr( Ufmri,U22 );

lambda=[0 0.000001 0.001 0.1 1 1e10 ];
for i=1:size(lambda,2)
    [ sol(i)] = coupled_tensors_l( Teeg,Tfmri,R1,R2,'cpd',X,lambda(i),U1,U2);
    U1{1,1}=sol(i).factors.A;
    U1{1,2}=sol(i).factors.B;
    U2{1,1}=sol(i).factors.D;
    U2{1,2}=sol(i).factors.E;    
    [ mcorre2(i),tcorre2(i,:),scorre2(i,:)] = mean_corr( Ueeg,U1 );
    [ mcorrf2(i),tcorrf2(i,:),scorrf2(i,:) ] = mean_corr( Ufmri,U2 );
    U1{1,3}=sol(i).factors.C;
    U2{1,3}=sol(i).factors.F;
end
temp=sprintf('/home/v1cchatz/Dropbox/slall_smooth.mat',i);
save(temp,'sol')
temp=sprintf('/home/v1cchatz/Dropbox/mocorre_smooth.mat',i);
save(temp,'mcorre2')
temp=sprintf('/home/v1cchatz/Dropbox/mcorrf_smmooth.mat',i);
save(temp,'mcorrf2')

