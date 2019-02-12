for j=1:30
    load('signal_hrf.mat')
    hrf=sig_hrf(:,1)';

    [tfmri,t_eeg,sig_erp]=Simulate_time(1,0,hrf);
    [sfmri,seeg]=Simulate_space();
    t_eeg(957:960,:)=0;
    
    Nt=4;
    
    R1=6;
    R2=6;
    
    U1{1}=sfmri';
    U1{2}=tfmri';
    U2{1}=seeg';
    U2{2}=t_eeg;
    U2{3}=sig_erp';
    
    U1{3}=abs(randn(6,6));
    U2{4}=U1{3};
    %%%%% Create a matrix for the coupling%%%%%
    %Toeplitz matrix created from hrf
    sizet=size(U2{1,2},1)/Nt;
    toepl_hrf=toeplitz([hrf(1) zeros(1,sizet-1)],[hrf zeros(1,sizet-1)]);
    
    %%Downsampling matrx
    d=eye(size(U2{1,2},1));
    downs=d(1:Nt:end,:);
    
    d=zeros(1,340000);

    %%downsampl+hrf matrix
    X=toepl_hrf'*downs;


    T1=cpdgen(U1);
    T1=T1+mean(mean(std(T1)))*sqrt(30)*randn(size(T1));

    T2=cpdgen(U2);
    T2=T2+mean(mean(mean(std(T2))))*sqrt(10)*randn(size(T2));


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
    [ec1(j),et1,es1]=mean_corr(U1,U11);

    %% Only T2
    model3 = struct;
    model3.variables.d = iU2{1};
    model3.variables.e = iU2{2};
    model3.variables.f = iU2{3};
    model3.variables.g=iU2{4};
    model3.factors.D = 'd';
    model3.factors.E ='e';
    model3.factors.F ='f';
    model3.factors.G='g';
    model3.factorizations.T.data = T2;
    % model3.factorizations.T.data = {U,S};
    model3.factorizations.T.cpd  = {'D','E','F','G'};
    sdf_check(model3,'print');
    [sol3,output] = sdf_nls(model3,'Display', 100,'CGMaxIter', 200,'MaxIter', 200);
    U21{1,1}=sol3.factors.D;
    U21{1,2}=sol3.factors.E;
    [c1(j),t1,s1]=mean_corr(U2,U21);


    %%FIX permutatin and                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
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


    %% Hard coupling
    model2 = struct;
    % {A, N}
    model2.variables.a = iU1{1};
    model2.variables.c = iU1{3};
    model2.variables.d = iU2{1};
    model2.variables.e = iU2{2};
    model2.variables.f= iU2{3};
    model2.variables.g = iU2{4};
    model2.factors.A = 'a';
    model2.factors.B ={'e',@(z,task) struct_matvec(z,task,X)};
    model2.factors.C ='c';
    model2.factors.D = 'd';
    model2.factors.E='e';
    model2.factors.F='f';
    model2.factors.G='g';
    model2.factorizations.T.data = T1;
    model2.factorizations.T.cpd  = {'A','B','C'};
    model2.factorizations.S.data = T2;
    % model3.factorizations.S.data = {U,S};
    model2.factorizations.S.cpd  = {'D','E','F','G'};
    sdf_check(model2,'print');
    [sol2,output] = sdf_nls(model2,'Display', 100,'CGMaxIter', 200,'MaxIter', 200);
    U22{1,1}=sol2.factors.D;
    U22{1,2}=sol2.factors.E;
    [c2(j),t2,s2]=mean_corr(U2,U22);
    U12{1,1}=sol2.factors.A;
    U12{1,2}=sol2.factors.B;
    [ec2(j),et2,es2]=mean_corr(U1,U12);

    %% Soft coupling
    lambda= [ 0.0001 0.01 0.1 1 10 100 10000 1000000];
    model2 = struct;
    model2.variables.a = iU1{1};
    model2.variables.e = {iU2{2},randn(size(T2,2),R1)*1e-2};
    model2.variables.c = iU1{3};
    model2.variables.d = iU2{1};
    model2.variables.f= iU2{3};
    model2.variables.g= iU2{4};
    model2.factors.A = 'a';
    model2.factors.B ={'e',@struct_plus, @(z,task) struct_matvec(z,task,X)};
    model2.factors.C ='c';
    model2.factors.D = 'd';
    model2.factors.E={'e',@(z,task) struct_select(z,task,1)};
    model2.factors.F='f';
    model2.factors.G='g';
    model2.factors.N={'e',@(z,task) struct_select(z,task,2)};
    model2.factorizations.T.data = T1;
    model2.factorizations.T.cpd  = {'A','B','C'};
    model2.factorizations.S.data = T2;
    model2.factorizations.S.cpd  = {'D','E','F','G'};
    model2.factorizations.N.regL2={'N'};
    for i=1:size(lambda,2)
        fprintf('Testing lamba %i\n', i)
        model2.factorizations.N.relweight=lambda(i);
        sdf_check(model2,'print');
        [sol2,output] = sdf_nls(model2,'Display', 100,'CGMaxIter', 200,'MaxIter', 200);
        U23{1,1}=sol2.factors.D;
        U23{1,2}=sol2.factors.E;
        [c3(j,i,:),t3,s3]=mean_corr(U2,U23);
        U13{1,1}=sol2.factors.A;
        U13{1,2}=sol2.factors.B;
        [ec3(j,i,:),et3,es3]=mean_corr(U1,U13);
    end
   
temp=sprintf('diffhrf%d.mat',j);
save(temp,'c1','ec1','c2','ec2','c3','ec3')
end

