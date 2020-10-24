function [R, S, D] = crocco_method(sqT, c)
    [M, K] = size(sqT);

    T   = (sqT * c).^2;
    T   = bsxfun(@minus, T, T(:, 1)); 
    T   = bsxfun(@minus, T, T(1, :));
    T   = T(2:end, 2:end);
    D   = (sqT * c).^2;
    [U,Sigma,V] = svd(T);
    Sigma = Sigma(1:3, 1:3);
    U     = U(:, 1:3);
    V     = V(:, 1:3);
    a1 = sqT(1,1) * c;
    opt = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e6);
    MAX_ITER = 0;
    [Cbest, costbest] = fminsearch(@(C) costC2(C, U, Sigma, V, D, a1), randn(3), opt);
    for i = 1:MAX_ITER
        [C, costval] = fminsearch(@(C) costC2(C, U, Sigma, V, D, a1), randn(3), opt);
        if costval < costbest
            costbest = costval;
            Cbest = C;
        end
    end
    C = Cbest;

    tilde_R = (U*C)';
    tilde_S = -1/2 * C\(Sigma*V');
    R = [ [ 0 0 0]' tilde_R ];
    tilde_S(1, :) = tilde_S(1, :) + a1;
    S = [ [a1 0 0]' tilde_S ];
    D = edm([R S], [R S]);

end
