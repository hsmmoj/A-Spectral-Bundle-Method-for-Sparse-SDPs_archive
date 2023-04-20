function [At,b,c,K, OutSpa] = interconnectedCones(m,nCones,nBlk,BlkSize,overlap, ArrowHead)

% Setup problem for block-arrow sdp with multiple cones. Inputs:
% m        : number of equality constraints
% nCones   : number of cones with same block-arrow sparsity pattern
% nBlk     : number of diagonal blocks in block-arrow matrix data (vector of length nCones)
% BlkSize  : size of each diagonal block (vector of length nCones)
% ArrowHead: size of head of arrow pattern (vector of length nCones)

% cone
K.f = 0;
K.l = 0;
K.q = 0;
K.s = zeros(1,nCones);

% Sparsity pattern of each cone
Spa = cell(nCones,1);
for k = 1:nCones
    n = nBlk(k)*BlkSize(k) - (nBlk(k) - 1) * overlap(k) + ArrowHead(k);
    Spa{k} = zeros(n);
    for i = 1:nBlk(k)
        Spa{k}((i-1)*BlkSize(k)+1 - (i-1)*overlap(k):BlkSize(k)*i -(i-1)*overlap(k), ...
               (i-1)*BlkSize(k)+1 - (i-1)*overlap(k):BlkSize(k)*i -(i-1)*overlap(k)) = ones(BlkSize(k));
    end
    Spa{k}(nBlk(k)*BlkSize(k)+1 -(nBlk(k) - 1) * overlap(k):n,:) = 1;
    Spa{k}(:,nBlk(k)*BlkSize(k)+1 -(nBlk(k) - 1) * overlap(k):n) = 1;
    
    % Set cone size
    K.s(k) = n;
end
OutSpa = Spa{1};

% Data
At = [];
for i = 1:m
    Ai = [];
    for k = 1:nCones
        M = 100*sprandsym(Spa{k});  % random symmetric data with given sparsity pattern
        Ai = [Ai; M(:)];        % concatenate
    end
    At = [At, Ai(:)];
end


% stictly feasible primal point
X = cell(nCones,1);
for k = 1:nCones
    Temp = 10*sprandsym(Spa{k});
    Temp = Temp + (-min(eig(full(Temp)))+1)*speye(size(Spa{k}));
    X{k} = Temp(:);
end
b = At'*vertcat(X{:});

% stictly feasible dual point
y = rand(m,1);
S = cell(nCones,1);
for k = 1:nCones
    Temp = 10*sprandsym(Spa{k});
    Temp = Temp + (-min(eig(full(Temp)))+1)*speye(size(Spa{k}));
    S{k} = Temp(:);
end
c = vertcat(S{:}) + At*y;

end