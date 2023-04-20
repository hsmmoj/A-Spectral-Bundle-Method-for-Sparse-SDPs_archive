% Parameters
m = 100;                        % # constraints
nCones = 1;                     % # cones with block-arrow sparsity pattern
nBlk = [15, 20];                % # diagonal blocks for each cone
BlkSize = [10, 5];              % block size for each cone
ArrowHead = [5, 5];             % arrow head size for each PSD cone

% Setup
fprintf('\nSetting up random block-arrow SDP, m=%i...',m);
tsetup = tic;
[At,b,c,K] = blockArrowMultCones(m,nCones,nBlk,BlkSize,ArrowHead);
tsetup = toc(tsetup);
fprintf('done in %.2f seconds. \n',tsetup);





% -------------------
% blockArrowMultCones
% -------------------
function [At,b,c,K,OutSpa] = blockArrowMultCones(m,nCones,nBlk,BlkSize,ArrowHead)

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
    n = nBlk(k)*BlkSize(k)+ArrowHead(k);
    Spa{k} = zeros(n);
    for i = 1:nBlk(k)
        Spa{k}((i-1)*BlkSize(k)+1:BlkSize(k)*i, ...
               (i-1)*BlkSize(k)+1:BlkSize(k)*i) = ones(BlkSize(k));
    end
    Spa{k}(nBlk(k)*BlkSize(k)+1:n,:) = 1;
    Spa{k}(:,nBlk(k)*BlkSize(k)+1:n) = 1;
    
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

Spa{1}((end-ArrowHead(1)+1):end,:) = 0;
Spa{1}(:,(end-ArrowHead(1)+1):end) = 0;

% stictly feasible primal point
X = cell(nCones,1);
for k = 1:nCones
    Temp = 10*sprandsym(Spa{k});
    Temp = Temp + (-min(eig(full(Temp)))+1)*speye(size(Spa{k}));
%     Temp((end-ArrowHead(1)+1):end,:) = 0;
%     Temp(:,(end-ArrowHead(1)+1):end) = 0;
    X{k} = Temp(:);
end

b = At'*vertcat(X{:});

[h,w] = size(mat(X{1}));
Spa{1} = zeros(h,w);
Spa{1}(end-ArrowHead(1)+1:end,end-ArrowHead(1)+1:end) = 1;


% stictly feasible dual point
y = rand(m,1);
S = cell(nCones,1);
for k = 1:nCones
    Temp = 10*sprandsym(Spa{k});
    Temp = Temp + (-min(eig(full(Temp)))+1)*speye(size(Spa{k}));
    Temp(1:(end-ArrowHead(1)),1:(end-ArrowHead(1))) = 0;
    S{k} = Temp(:);
end
c = vertcat(S{:}) + At*y;

end

