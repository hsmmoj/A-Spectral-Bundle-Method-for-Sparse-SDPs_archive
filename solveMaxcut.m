function [P_star] = solveMaxcut(A,B,C)
yalmip('clear')
dim = size(A,1);
x_var = sdpvar(dim,dim);
objective = trace(C' * x_var);
constraints = [x_var>=0];
for i=1:length(B)
    constraints = [constraints, trace(A(:,:,i)' * x_var) == B(i)];
end
options = sdpsettings('verbose',1,'solver','mosek');
sol = optimize(constraints,objective,options);
if sol.problem == 0
 P_star = value(objective) ;
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end

end