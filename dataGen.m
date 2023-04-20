function [A,B,C] = dataGen(opt)
%rng(opt.seed,'twister')
rng('default');
%Generating A,B,C
A = zeros(opt.dim,opt.dim,opt.numOfConstraints);
B = zeros(opt.numOfConstraints,1);
C = zeros(opt.dim);
y = rand(opt.numOfConstraints,1);
dummy = rand(opt.dim);
X_0 = ((dummy * dummy') .* opt.sparsityPat);
X_0 = X_0 + eye(opt.dim) * abs(min(eig(X_0)));
dummy = rand(opt.dim); 
Z_0 = ((dummy * dummy') .* opt.sparsityPat);
Z_0 = Z_0 + eye(opt.dim) * abs(min(eig(Z_0)));
for i = 1:size(A,3)
    dummyA = rand(opt.dim);
    A(:,:,i) = (dummyA + dummyA') .* opt.sparsityPat;
    B(i) = trace(X_0 * A(:,:,i)');
end
for i=1:opt.numOfConstraints
    C = C + y(i) * A(:,:,i);
end
C = C + Z_0;
end