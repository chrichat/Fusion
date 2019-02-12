function [ sol, sol2, sol3] = coupled_tensors( T1,T2,R1,R2,Nt,method,hrf,lambda1,lambda2 )
%Coupled tensor decomposition methods.
% -T1 is the tensor of fMRI and R1 its estimated rank.
% -T2 is the tensor of EEG and R2 its estimated rank.
% -Nt is the downsampling factor from the eeg to fmri, has to do with the
% frequency of the sampling in EEG and the TR in fMRI. 
% -method can be 'cpd' for coupled cpds and 'cpdparafac2' for cpd for fMRI and pfac2 for EEG .
% -hrf is the hrf function that will be used - in a latter version it shall
% be also optimized

if method=='cpd'
    m=1;
elseif method=='cpdparafac2'
    m=2;
else
    error ('The specified method for coupled tensors is not applicable. Options are cpd and parafac2');
end


T1=T1/max(max(max(T1))); %Normalize both tensors 
T2=T2/max(max(max(T2)));

%Toeplitz matrix created from hrf 
sizet=size(T1,2)/Nt;
toepl_hrf=toeplitz([hrf(1) zeros(1,sizet-1)],[hrf zeros(1,sizet-1)]);

%%Downsampling matrx
d=eye(size(T1,2));
downs=d(1:Nt:end,:);

%%downsampl+hrf matrix
X=toepl_hrf'*downs;

[U1,output] = cpd_gevd(T1,R1);
[U2,output] = cpd_gevd(T2,R2);

%% Only T1
model2 = struct;
model2.variables.a = U1{1};
model2.variables.b = U1{2};
model2.variables.c = U1{3};
model2.factors.A = 'a';
model2.factors.B ='b';
model2.factors.C ='c';
model2.factorizations.T.data = T1;
model2.factorizations.T.cpd  = {'A','B','C'};
[sol2,output]=sdf_nls(model2,'Display',100,'CGMaxIter',500,'TolFun', eps^2, 'TolX', eps, 'MaxIter', 300);

%% Only T2
model3 = struct;
model3.variables.d = U2{1};
model3.variables.e = U2{2};
model3.variables.f = U2{3};
model3.factors.D = 'd';
model3.factors.E ='e';
model3.factors.F ='f';
model3.factorizations.T.data = T2;
model3.factorizations.T.cpd  = {'D','E','F'};
sdf_check(model3,'print');
[sol3,output]=sdf_nls(model3,'Display',100,'CGMaxIter',500,'TolFun', eps^2, 'TolX', eps, 'MaxIter', 300);


if R1==R2
    R=R1;
    model=struct;
    model.variables.a=U1{1};
%     model.variables.b=randn(size(T1,2),R); %Uncomment for hard coupling
    model.variables.c=randn(size(T1,3),R);
    model.variables.be={U1{2},randn(size(T1,2),R)*1e-2}; %comment for hard
    model.variables.d=U2{1};
    model.variables.e=randn(size(T2,2),R);
    model.factors.A='a';
    model.factors.B='b'; %Those two lines for hard coupling
    model.factors.C='c';
%     model.factors.B={'be',@(z,task) struct_select(z,task,1)};
%     model.factors.C={'cf',@(z,task) struct_select(z,task,1)};
    model.factors.D='d';
    model.factors.F='c';  %% Those two lines for hard coupling
    model.factors.E={'b',@(z,task) struct_matvec(z,task,X)};
%     model.factors.E={'be',@struct_plus, @(z,task) struct_matvec(z,task,X)};
%     model.factors.E='e';
%     model.factors.F={'cf',@struct_plus};
%     model.factors.N1={'be',@(z,task) struct_select(z,task,2)};
%     model.factors.N2={'cf',@(z,task) struct_select(z,task,2)};
    model.factorizations.T.data=T1;
    model.factorizations.T.cpd={'A','B','C'};
    model.factorizations.M.data=T2;
    model.factorizations.M.cpd={'D','E','F'};
%     model.factorization.N1.regL2={'N1'};
%     model.factorization.N2.regL2={'N2'};
    [sol,output]=sdf_nls(model,'Display',100,'CGMaxIter',500,'TolFun', eps^2, 'TolX', eps, 'MaxIter', 300);
elseif R1>R2
    model=struct;
    model.variables.a=randn(size(T1,1),R1);
    model.variables.be={randn(size(T1,2),R2),randn(size(T1,2),R2)*1e-2};
    model.variables.b2={randn(size(T1,2),R1-R2)};
    model.variables.cf={randn(size(T1,3),R2),randn(size(T1,3),R2)*1e-2};
    model.variables.c2={randn(size(T1,3),R1-R2)};
    model.variables.d=randn(size(T2,R2),5);
    model.factors.A='a';
    model.factors.B={{'be',@(z,task) struct_select(z,task,1)},{'b2'}};
    model.factors.C={{'cf',@(z,task) struct_select(z,task,1)},{'c2'}};
    model.factors.D='d';
    model.factors.E={'be',@struct_plus, @(z,task) struct_matvec(z,task,X)};
    model.factors.F={'cf',@struct_plus};
    model.factors.N1={'be',@(z,task) struct_select(z,task,2)};
    model.factors.N2={'cf',@(z,task) struct_select(z,task,2)};
    model.factorizations.T.data=T1;
    model.factorizations.T.cpd={'A','B','C'};
    model.factorizations.M.data={U,S};
    model.factorizations.M.cpd={'D','E','F'};
    model.factorization.N1.regL2={'N1'};
    model.factorization.N2.regL2={'N2'};
    [sol,output]=sdf_nls(model,'Display',10,'CGMaxIter',500,'TolFun', eps^2, 'TolX', eps, 'MaxIter', 300);
else
    model=struct;
    model.variables.a=randn(size(T1,1),R1);
    model.variables.be={randn(size(T1,2),R1),randn(size(T1,2),R1)*1e-2};
    model.variables.b2={randn(size(T1,2),R2-R1)};
    model.variables.cf={randn(size(T1,3),R2),randn(size(T1,3),R2)*1e-2};
    model.variables.c2={randn(size(T1,3),R2-R1)};
    model.variables.d=randn(size(T2,R2),5);
    model.factors.A='a';
    model.factors.B={'be',@(z,task) struct_select(z,task,1)};
    model.factors.C={'cf',@(z,task) struct_select(z,task,1)};
    model.factors.D='d';
    model.factors.E={{'be',@struct_plus, @(z,task) struct_matvec(z,task,X)},'b2'};
    model.factors.F={{'cf',@struct_plus},'c2'};
    model.factors.N1={'be',@(z,task) struct_select(z,task,2)};
    model.factors.N2={'cf',@(z,task) struct_select(z,task,2)};
    model.factorizations.T.data=T1;
    model.factorizations.T.cpd={'A','B','C'};
    model.factorizations.M.data=T2;
    model.factorizations.M.cpd={'D','E','F'};
    model.factorization.N1.regL2={'N1'};
    model.factorization.N2.regL2={'N2'};
    [sol,output]=sdf_nls(model,'Display',10,'CGMaxIter',500,'TolFun', eps^2, 'TolX', eps, 'MaxIter', 200);
end

end

