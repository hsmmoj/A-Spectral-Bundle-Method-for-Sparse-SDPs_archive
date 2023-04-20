function [P_star, rho, Ec] = computeRho(A,B,C,opts)
    addpath('C:\Users\SOC-LAB\FengYi\Solver\CDCS-master\packages\+cdcs_utils');
    import cdcs_utils.*
    K.l = 0;
    K.q = 0;
    K.f = 0;
    K.s = sqrt(length(opts.At));
    [~,~,Ats,Cs] = svecData(opts.At,opts.c,K);
    temp          = chordalDecomposition(Ats,Cs,K);
    cliques = temp{1}.Set;
    for i=1:length(cliques)
        cliques{i} = cliques{i} - 1;
    end


    %Extract cliques and make "E" matrix  
    %cliques          = pyrunfile("clique.py","res",E=opts.sparsityPat,dim=opts.dim);
    %cliques          = cellfun(@double,cell(cliques),'UniformOutput',false);
    %cliques = {[0,1,2], [3,4,5,6]};
    %cliques = {linspace(0,99,100)};
    %cliques = opts.cliques;
    
    
    numOfCliques     = length(cliques);
    Ec = {};
    for k=1:length(cliques)
        dummy        = zeros(length(cliques{k}),opts.dim);
        for i=1:length(cliques{k})
            dummy(i,cliques{k}(i)+1) = 1;
        end
        Ec           = [Ec; dummy];
    end
    %Solve the problem using convetional solvers, i.e. Mosek, Sedumi. 
  [P_star, rho]           = OriginalProblem(A,B,C,Ec);