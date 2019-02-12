function [ sol] = coupled_tensors( T1,T2,R1,R2,method,X,lambda1,U1,U2 )
%Coupled tensor decomposition methods.
% -T1 is the tensor of fMRI and R1 its estimated rank.
% -T2 is the tensor of EEG and R2 its estimated rank.
% -Nt is the downsampling factor from the eeg to fmri, has to do with the
% frequency of the sampling in EEG and the TR in fMRI. 
% -method can be 'cpd' for coupled cpds and 'parafac2' for coupled pfac2s and 'cpdparafac2' for cpd for EEG and pfac2 for fMRI .
% -hrf is the hrf function that will be used - in a latter version it shall
% be also optimized

if method=='cpd'
    m=1;
elseif method=='parafac2'
    m=2;
elseif method=='cpdparafac2'
    m=3;
else
    error ('The specified method for coupled tensors is not applicable. Options are cpd and parafac2');
end

T1=T1/max(max(max(abs(T1)))); %Normalize both tensors 
T2=T2/max(max(max(abs(T2))));

R=R1;
model=struct;
model.variables.a=U1{1};
%     model.variables.b=randn(size(T1,2),R); %Those two lines for hard
model.variables.c=U1{3};
model.variables.be={U1{2},randn(size(T1,2),R)*1e-2};
%     model.variables.cf={U1{3},randn(size(T1,3),R)*1e-2};
model.variables.d=U2{1};
%     model.variables.e=randn(size(T2,2),R);
model.variables.f=U2{3};
model.factors.A='a';
%     model.factors.B='b'; %Those two lines for hard coupling
model.factors.C='c';
model.factors.B={'be',@(z,task) struct_select(z,task,1)};
%     model.factors.C={'cf',@(z,task) struct_select(z,task,1)};
model.factors.D='d';
model.factors.F='f';
%     model.factors.F='c';  %% Those two lines for hard coupling
%     model.factors.E={'b',@(z,task) struct_matvec(z,task,X)};
model.factors.E={'be',@struct_plus, @(z,task) struct_matvec(z,task,X)};
%     model.factors.E='e';
%     model.factors.F={'cf',@struct_plus};
model.factors.N1={'be',@(z,task) struct_select(z,task,2)};
%     model.factors.N2={'cf',@(z,task) struct_select(z,task,2)};
model.factorizations.T.data=T1;
model.factorizations.T.cpd={'A','B','C'};
model.factorizations.T.relweight=1;
model.factorizations.M.data=T2;
model.factorizations.M.cpd={'D','E','F'};
model.factorizations.M.relweight=1;
model.factorizations.N1.regL2={'N1'};
model.factorizations.N1.relweight=lambda1;
%     model.factorization.N2.regL2={'N2'};
[sol,output]=sdf_nls(model,'Display',100,'CGMaxIter',500,'TolFun', eps^2, 'TolX', eps, 'MaxIter', 300);


end
