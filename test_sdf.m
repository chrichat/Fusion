I=15;
R=4;
U=cpd_rnd([I I I],R);
U{1}(1:R,1:R)=toeplitz(rand(4,1));
T=cpdgen(U);
model.variables.atop=randn(2*R-1,1);
orthq=@(z,task)struct_orth(z,task,[I,R]);
model.variables.abtm=randn(I-R,R);
model.variables.b=randn(I,R);
model.factors.A={{'atop',orthq};{'abtm'}};
model.factors.B='b';
model.factors.C=U{3};
model.factorizations.myfac.data=T;
model.factorizations.myfac.cpd={'A','B','C'};
options.Display=0;
sol=sdf_nls(model,options)
