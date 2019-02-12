function [ sol ] = joint_btd( data1,data2,R,L,varargin )
%   Detailed explanation goes here
    sizet=size(data1);

    if sizet~=size(data2)
       error('Datasets of not the same size.')    
    end

    options.Display = 0; % View convergence progress every 5 iterations.
    options.TolFun=1e-9;
    options.Tolx=eps;
    nargin
    if nargin<3
        error('The minimum number of inputs is 3 (data,L and R).')
    end
    
    for i=1:R
        model.variables.(sprintf('A%d',i)) = randn(sizet(1),L);
        model.variables.(sprintf('B%d',i)) = randn(sizet(2),L);
        model.variables.(sprintf('c%d',i))= randn(sizet(3),1); 
        model.variables.(sprintf('kA%d',i)) = randn(sizet(1),L);
        model.variables.(sprintf('kB%d',i)) = randn(sizet(2),L);
    end
    model.factors.A1 = strsplit(sprintf('A%d ',1:R));
    model.factors.A1=model.factors.A1(1:end-1);
    model.factors.B1 = strsplit(sprintf('B%d ',1:R));
    model.factors.B1=model.factors.B1(1:end-1);
    model.factors.kA1 = strsplit(sprintf('kA%d ',1:R));
    model.factors.kA1=model.factors.kA1(1:end-1);
    model.factors.kB1 = strsplit(sprintf('kB%d ',1:R));
    model.factors.kB1=model.factors.kB1(1:end-1);
    model.factors.C1 = strsplit(sprintf('c%d ',ceil((1:L*R)/L)));
    model.factors.C1=model.factors.C1(1:end-1);    
    model.factorizations.imag.data = data1;
    model.factorizations.imag.cpd = {'A1','B1','C1'};
    model.factorizations.kspace.data = data2;
    model.factorizations.kspace.cpd = {'kA1','kB1','C1'};
    % Solve the SDF model.
    sol = sdf_nls(model,options);
end


