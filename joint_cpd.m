function [ sol ] = joint_cpd( data1,R )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
sizet=size(data1);
options.Display = 0; % View convergence progress every 5 iterations.
options.TolFun=1e-6;
options.Tolx=eps;

model.variables.a=randn(size(data1,1),R);
model.variables.b=randn(size(data1,2),R);
model.variables.c=randn(size(data1,3),R);
model.factors.A='a';
model.factors.B='b';
model.factors.C='c';
model.factorizations.mycpd.data=data1;
model.factorizations.mycpd.cpd={'A','B','C'};
model.factorizations.myreg1.regL1 = {'B'};
model.factorizations.myreg1.relweight = 0.5;
sol = sdf_nls(model,options);

end



