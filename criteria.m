function [decent,stop,realValAtOmega,opts] = criteria(omega,X_star,C,numOfCliques,cliques,W_val,opts,iter)

    realValAtOmega = trace(C'*omega);
    realValAtXstar = trace(C'*X_star);
    appValAtXstar = realValAtXstar;
    for i=1:numOfCliques
        smallX = subMatrixExt(X_star, cliques{i});
        [eigVecX, eigValX] = eig(-smallX,'vector');
        smallOmega = subMatrixExt(omega, cliques{i});
        [eigVecOmega, eigValOmega] = eig(-smallOmega,'vector');
        realValAtOmega = realValAtOmega +  opts.rho * max(0,max(eigValOmega));
        realValAtXstar =  realValAtXstar + opts.rho * max(0,max(eigValX));
          appValAtXstar = appValAtXstar - trace(W_val{i}' * smallX);
    end
    stop = 0;
    decent = 0;
    null = 0;
    if opts.beta * (realValAtOmega - appValAtXstar) <= realValAtOmega - realValAtXstar
        decent = 1;
        if  (opts.mu * (realValAtOmega - appValAtXstar) <= realValAtOmega - realValAtXstar)  && opts.dynamicAlpha && null>5
            opts.alpha = max(10^-5,opts.alpha/1.01);
            null = 0;
        end
%         if iter > 100 
%             opts.alpha = 0.0001;
%         elseif null == 10 
%             opts.alpha = 0.1;
%             null = 0;
%         end
    else
        if opts.ml * (realValAtOmega - appValAtXstar) >= (realValAtOmega - realValAtXstar)  && opts.dynamicAlpha
           opts.alpha = min(opts.alpha*1.01,1);
        end
        null = null + 1;
    end
    if (realValAtOmega - appValAtXstar) <= opts.eps
        stop = 1;
    end
    fprintf('%7.6f \n',realValAtOmega- appValAtXstar);
end
