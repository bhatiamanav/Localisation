function [X,D_edm] = Stress_Alg(S,M,dim)
    %S = 2 x m or 3 x m Array
    %M = 2 x n or 3 x n Array
    %dim = Dimension for Alternating Descent
    sz1 = size(S);
    sz2 = size(M);
    m = sz1(2);
    n = sz2(2);
    D_crt = zeros(n+m,n+m);
    iter_num = 60;
    iter = 0;
    min_err = 1e-6;
    
    while iter < iter_num
        D_prev = D_crt;
        
        for ii = (n+1) : (n+m)
            for jj = 1 : n
                D_crt(ii,jj) = (norm(S(:,ii - n)-M(:,jj)));
                D_crt(jj,ii) = (norm(S(:,ii - n)-M(:,jj)));
            end
        end
        
        [X, D] = alternating_descent_sumanth(D_crt, dim);
        err = norm(D_prev - D, 'fro');
        
        if err < min_err
            break;
        end
        
        D_crt = D;
        iter = iter + 1;
    end
    
   D_edm = D_crt;
   X = classic_mds(D_edm, dim);
end