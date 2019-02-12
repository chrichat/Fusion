load('signal_hrf.mat')
hrfn=sig_hrf(:,1)'; %%Change this for different hrf per subject

[tfmri,t_eeg,sig_erp]=Simulate_time(1,0,hrfn);
[sfmri,seeg]=Simulate_space();
t_eeg(957:960,:)=0;

hrf=sig_hrf(:,1)';
Nt=4; %supsampling rate

R1=6;
R2=6;

U1{1}=sfmri';
U1{2}=tfmri';
U2{1}=seeg';
U2{2}=t_eeg;
% U2{3}=sig_erp'; % Add erps

U1{3}=abs(randn(6,6));
% U2{4}=U1{3}; % If erps are used
U2{3}=U1{3}; % If erps are not used


%%%%% Create a matrix for the coupling%%%%%
%Toeplitz matrix created from hrf
sizet=size(U2{1,2},1)/Nt;
toepl_hrf=toeplitz([hrf(1) zeros(1,sizet-1)],[hrf zeros(1,sizet-1)]);

%%Downsampling matrx
d=eye(size(U2{1,2},1));
downs=d(1:Nt:end,:);

%%downsampl+hrf matrix
X=toepl_hrf'*downs;

clearvars -except U1 U2 R1 R2 X
lambda=1;
for j=1:10
    j
    T1=cpdgen(U1);
    T1=T1+mean(mean(std(T1)))*sqrt(10)*randn(size(T1));

    T2=cpdgen(U2);
    T2=T2+mean(mean(mean(std(T2))))*sqrt(1)*randn(size(T2));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    T1=T1/max(max(max(T1))); %Normalize both tensors 
    T2=T2/max(max(max(max(T2))));

    % [U,S]=mlsvd(T2,[size(T2,2)*size(T2,3),size(T2,2),size(T2,3)]);
    % [U,S]=mlsvd(T2,[R,R,R]);

    [iU1,output] = cpd_gevd(T1,R1);
    [iU2,output] = cpd_gevd(T2,R2);

    %% Only T1
    model1 = struct;
    model1.variables.a = iU1{1};
    model1.variables.b = iU1{2};
    model1.variables.c = iU1{3};
    model1.factors.A= 'a';
    model1.factors.B ='b';
    model1.factors.C ='c';
    model1.factorizations.T.data = T1;
    % model3.factorizations.T.data = {U,S};
    model1.factorizations.T.cpd  = {'A','B','C'};
    sdf_check(model1,'print');
    [sol1,output] = sdf_nls(model1,'Display', 100,'CGMaxIter', 200,'MaxIter', 200);
    U11{1,1}=sol1.factors.A;
    U11{1,2}=sol1.factors.B;
    [fc1(j),et1,es1]=mean_corr(U1,U11);

    %% Only T2
    model3 = struct;
    model3.variables.d = iU2{1};
    model3.variables.e = iU2{2};
    model3.variables.f = iU2{3};
    model3.factors.D = 'd';
    model3.factors.E ='e';
    model3.factors.F ='f';
%     model3.factors.G='g';
    model3.factorizations.T.data = T2;
    % model3.factorizations.T.data = {U,S};
    model3.factorizations.T.cpd  = {'D','E','F'};
    sdf_check(model3,'print');
    [sol3,output] = sdf_nls(model3,'Display', 100,'CGMaxIter', 200,'MaxIter', 200);
    U21{1,1}=sol3.factors.D;
    U21{1,2}=sol3.factors.E;
    [ec1(j),t1,s1]=mean_corr(U2,U21);


    %%FIX permutation and
    %%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    tempcorr=abs(corr(X*iU2{1,2},iU1{1,2}));
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
    iU1{1,1}(:,ind)=iU1{1,1};
    iU1{1,2}(:,ind)=iU1{1,2};
    iU1{1,3}(:,ind)=iU1{1,3};

    %Fix signs
    tempcorr=abs(corr(X*iU2{1,2},iU1{1,2}));
    for i=1:size(iU1{1,2},2)
        if tempcorr(i,i)<0

            iU1{1,2}(:,i)=-iU1{1,2}(:,i);
        end
    end
    %%%


    %% Soft coupling
    model2 = struct;
    % {A, N}
    model2.variables.a = iU1{1};
    model2.variables.e = {iU2{2},randn(size(T2,2),R1)*1e-2};
    model2.variables.c = iU1{3};
    model2.variables.d = iU2{1};
    model2.variables.f= iU2{3};
    model2.factors.A = 'a';
    model2.factors.B ={'e',@struct_plus, @(z,task) struct_matvec(z,task,X)};
    model2.factors.C ='c';
    model2.factors.D = 'd';
    model2.factors.E={'e',@(z,task) struct_select(z,task,1)};
    model2.factors.F='f';
    model2.factors.N={'e',@(z,task) struct_select(z,task,2)};
    model2.factorizations.T.data = T1;
    model2.factorizations.T.cpd  = {'A','B','C'};
    model2.factorizations.S.data = T2;
    model2.factorizations.S.cpd  = {'D','E','F'};
    model2.factorizations.N.regL2={'N'};
    model2.factorizations.N.relweight=lambda;
    sdf_check(model2,'print');
    [sol2,output] = sdf_nls(model2,'Display', 100,'CGMaxIter', 200,'MaxIter', 200);
    U23{1,1}=sol2.factors.D;
    U23{1,2}=sol2.factors.E;
    [ec3(j,:),t3,s3]=mean_corr(U2,U23);
    U13{1,1}=sol2.factors.A;
    U13{1,2}=sol2.factors.B;
    [fc3(j),et3,es3]=mean_corr(U1,U13);
    
    
    %% Group ICA
    icafmri=reshape(T1,size(T1,1),[]);
    [S,A]=icaML(icafmri,R1);
    imag=icafmri*pinv(S);
    Uicaf{1,1}=imag;
    Uicaf{1,2}=mean(reshape(S,6,[],6),3)';
    [fcica(j),tica2,sica2]=mean_corr(U1,Uicaf);
    icaeeg=reshape(permute(T2,[2,1,3]),size(T2,2),[]);
    [S,A]=icaML(icaeeg,R1);  
    timeeg=icaeeg*pinv(S);
    Uicae{1,1}=mean(reshape(S,6,[],6),3)';
    Uicae{1,2}=timeeg;
    [ecica(j),tica2,sica2]=mean_corr(U2,Uicae);
end

