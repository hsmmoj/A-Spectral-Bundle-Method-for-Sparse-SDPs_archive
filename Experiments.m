clear;
clc;

% Set problem parameters 
opts.numOfConstraints = 1000;
opts.seed = 123;
opts.dynamicAlpha = true;
opts.alpha = 0.5;
opts.maxIter =50 ;
opts.eps =  10^-5 ;
opts.beta = 0.1;
opts.mu = 0.6;
opts.ml = 0.001;
opts.r_p = 0;
opts.r_c = 4;
  
              
nCones = 1;               % # cones with block-arrow sparsity pattern
nBlk = [40];              % # diagonal blocks for each cone
BlkSize = [50];           % block size for each cone
ArrowHead = [4];          % arrow head size for each PSD cone
iternconnection = [2];    % # overlapping elements for each two neighbour cones   

fprintf('\nSetting up random block-arrow SDP, m=%i...');
tsetup = tic;
%[At,b,c,K,sparsityPat] = interconnectedCones(opts.numOfConstraints, ...
%    nCones,nBlk,BlkSize,iternconnection,ArrowHead);
[At,b,c,K,sparsityPat,X] = GenStrictData(opts.numOfConstraints,nCones, ...
    nBlk,BlkSize,ArrowHead);
opts.sparsityPat = sparsityPat;
opts.dim = length(opts.sparsityPat);
A = reshape((full(At)), [opts.dim,opts.dim,opts.numOfConstraints]);
B = full(b);
C = reshape(full(c), [opts.dim,opts.dim]);

X = mat(full(X{1}));
X(end-ArrowHead+1:end,end-ArrowHead+1:end) = zeros(ArrowHead);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CDCS DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load("qap9.mat");
% opts.dim = sqrt(length(At));
% A = reshape((full(At)), [opts.dim,opts.dim,opts.numOfConstraints]);
% B = full(b);
% C = reshape(full(c), [opts.dim,opts.dim]);
% opts.sparsityPat = zeros(opts.dim);
% for i=1:size(A,3)
%     opts.sparsityPat = opts.sparsityPat | A(:,:,i) ~= 0;
% end
% opts.sparsityPat = opts.sparsityPat | C ~= 0;


%%%%%%%%%%%%%%%%%%%%%%%%% MAX CUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maxCUTdata = load('G11.mat');
% %C = full(maxCUTdata.Problem.A);
% C = -maxCUTdata.Problem.A;
% C = 0.25 * (diag(C *ones(size(C,1),1))- C);
% c_sdp = vec(C) ;
% [~,W] = size(maxCUTdata.Problem.A);
% At_sdp = sparse(W,W^2);
% for i = 1:W
%     At_sdp(i,(i-1)*W+i) = 1;
% end
% b_sdp = ones(W,1);
% K.s = W;
% [res1,mosektime]= SolveMosek(At_sdp,b_sdp,c_sdp,K);
% P_star = -res1.sol.itr.pobjval;
% disp('wait');
% opts.At = At_sdp';
% opts.b = b_sdp;
% opts.c = c_sdp;
% A = zeros(length(C),length(C),length(C));
% B = ones(length(C),1);
% for i=1:length(C)
%     A(i,i,i) = 1;
% end
%  opts.sparsityPat = eye(length(C)) | C~=0 ;
% opts.dim = length(opts.sparsityPat);
%[P_star] = solveMaxcut(A,B,C);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%POP DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('Broyden_sphere15_R1_Tran');
% At = A_new';
% b = b_new;
% c = c_new;
% opts.dim = sqrt(length(At));
% A = reshape((full(At)), [opts.dim,opts.dim,opts.numOfConstraints]);
% B = full(b);
% C = reshape(full(c), [opts.dim,opts.dim]);
% opts.sparsityPat = zeros(opts.dim);
% for i=1:size(A,3)
%     opts.sparsityPat = opts.sparsityPat | A(:,:,i) ~= 0;
% end
% opts.sparsityPat = opts.sparsityPat | C ~= 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


opts.At = At;
opts.b = b;
opts.c = c;
%[P_star, rho, Ec] = computeRho(A,B,C,opts);
P_star = trace(C'*X);
opts.rho =  trace(C);
[X_star, count, a] = BundleDecomposedSDP(A,B,C,opts);
semilogy(1:opts.maxIter-1,abs((a.Obj(2:end)-P_star))/abs(P_star))
Out.gap = a.gap;
Out.costOpt = (a.Obj(end) - P_star)/abs(P_star);
Out.affineOpt = a.dualAffine;
Out.semiOpt = 0;
Out.P_star = P_star;
Out.obj = a.Obj;
for i=1:length(Ec)
 temp = Ec{i} * X_star * Ec{i}'; 
 minEig = min(eig(temp));
 if minEig < Out.semiOpt
    Out.semiOpt = minEig;
 end
end
