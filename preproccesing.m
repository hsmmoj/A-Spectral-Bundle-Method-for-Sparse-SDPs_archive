function [cliques,W,omega, P,numOfConstraints, numOfCliques, overlaps, subA, subC, N,invQ33] = preproccesing(A,B,C,opts)
    addpath('C:\Users\SOC-LAB\FengYi\Solver\CDCS-master\packages\+cdcs_utils');
    import cdcs_utils.*
   
    %printHeader(opts);
    %opts.dim              = size(A,1);
    numOfConstraints = size(B,1);            % Will be added later. 
    %Extract cliques and make "E" matrix  
%     cliques          = pyrunfile("clique.py","res",E=opts.sparsityPat,dim=opts.dim);
%     cliques          = cellfun(@double,cell(cliques),'UniformOutput',false);
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
    overlaps = {};
    subA = {};
    subC = {};
    for i=1:length(cliques)
        for j=1:length(cliques)
            [~,ia,ib] = intersect(cliques{i},cliques{j});
            overlaps{end+1} = [ia ib];
        end
    end

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
   %[P_star, ~]           = OriginalProblem(A,B,C,Ec);
    % Generate Omega_0, P_0, and W_0
    omega            = eye(opts.dim);
    P = {};
    W = {};
    V = {};
    for i=1:length(cliques)
        cliques{i} = cliques{i} + 1;
        omegaSmall = subMatrixExt(omega, cliques{i});
        [eigVec, eigVal] = eig(-omegaSmall);
        [d,ind]      = sort(diag(eigVal), 'descend');
        Vs           = eigVec(:,ind);
        V{i}         = Vs(:,1:opts.r_c);
        P{i}         = Vs(:,1:opts.r_p+opts.r_c);
        W{i}         = Vs(:,1) * Vs(:,1)';
        tempA = zeros(length(cliques{i})^2, size(A,3));
        for j=1:size(A,3)
            tempA(:,j) = vec(subMatrixExt(A(:,:,j), cliques{i}));
        end
        subA{i} = tempA;
        subC{i} = vec(subMatrixExt(C, cliques{i}));
        
    end
    N = createN(A);
    
    U = chol(Ats.'*Ats);
    inv_U = inv(U);
    invQ33 = inv_U*inv_U.';



end