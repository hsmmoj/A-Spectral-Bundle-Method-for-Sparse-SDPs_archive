function [A,B,C,sparsityPat] = makeStrictData(nblk, blkSize, m, rankX)
    rng('default');
    dim = nblk * blkSize;
    sparsityPat = zeros(dim,dim);
    for i=1:nblk
        sparsityPat((i-1)*blkSize+1:(i)*blkSize, (i-1)*blkSize+1:(i)*blkSize) = ones(blkSize);
    end
    A = rand(dim,dim, m);
    B = zeros(m,1);
    C= zeros(dim,dim);
    y = rand(m,1);
    for i=1:m
        tmp = rand(dim);
        A(:,:,i) =  A(:,:,i) - trace(A(:,:,i)) * tmp/trace(tmp);
        A(:,:,i) = (A(:,:,i) + A(:,:,i)') .*sparsityPat;
        X = zeros(nblk*blkSize);
    end
    Z = {};
    for i=1:nblk
        tmp = rand(blkSize);
        tmp = tmp + tmp';
        tmp = tmp + 2* eye(blkSize) * abs(min(eig(tmp)));

     [v d] = eig(tmp);
     x_tmp = v * [d(:,1:rankX) zeros(blkSize, blkSize- rankX)] * v';
     x_tmp = x_tmp + x_tmp';
     z_tmp = v * [zeros(blkSize,rankX) d(:,rankX+1:end)] * v';
     z_tmp = z_tmp + z_tmp;
     X((i-1)*blkSize+1:(i)*blkSize, (i-1)*blkSize+1:(i)*blkSize) = x_tmp;
     Z{i} = z_tmp;
    end
    for i=1:m
        B(i) = trace(A(:,:,i)' * X);
        C  = C + y(i)*A(:,:,i);
    end
    for i=1:nblk
        C = C + matrixExp(Z{i},linspace((i-1)*blkSize+1,i*blkSize,blkSize),dim);

    end

end