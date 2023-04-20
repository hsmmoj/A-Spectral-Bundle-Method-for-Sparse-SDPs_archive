function [P_star, rho] = OriginalProblem(A,B,C,Ec)
addpath('C:\Users\Hesam\Desktop\UCSD\Research\codes\FengYiCode\SBMP-Decomposed-main\SeDuMi_1_3\SeDuMi_1_3');
yalmip('clear')
dualZ = {};
dim = size(A,1);
x_var = sdpvar(dim,dim);
objective = trace(C' * x_var);
constraints = [];
for i=1:length(B)
    constraints = [constraints, trace(A(:,:,i)' * x_var) == B(i)];
end

for i=1:length(Ec)
    constraints = [constraints, Ec{i} * x_var * Ec{i}' >= 0];
end
options = sdpsettings('verbose',0,'solver','mosek');
sol = optimize(constraints,objective,options);
if sol.problem == 0
 P_star = value(objective) ;
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end
for i=1:length(Ec)
    dualZ{i} = trace(dual(constraints(length(B)+i)));
    Z{i} = dual(constraints(length(B)+i));
end
dualZ = cell2mat(dualZ);
rho = max(dualZ);
end