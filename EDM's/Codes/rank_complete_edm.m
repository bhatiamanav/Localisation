function D = rank_complete_edm(t_D, W, dim, gram)
    %If gram=1,we rank threshold the gram matrix, else the EDM.
    %Alterntes between rank enforcing and known points enforcing
    N = 1000;%Max. number of iterations for a continuously decreasing sequence of errors
    Min_err  = 1e-7;%Criteria to stop if error is less than this(case of convergence)

    W = logical(W);%masking matrix
    n = size(t_D, 1);
    mean_d = mean(mean(sqrt(t_D)));
    D = t_D;%Incomplete EDM
    D(W) = mean_d;%Enforcing known points
    I = logical(eye(n));%Identity matrix
    J = I - (1/n)*ones(n);

    while true
        for i = 1:N
            D_prev = D;
            %Doing rank thresholding
            if gram == 1
                G = -1/2 * J * D * J;
                [U, S, V] = svd(G);
                S(dim+1:end, dim+1:end) = 0;%Utilising rank property of Gram matrices 
                G = U * S * V';
                D = diag(G)*ones(1, n) + ones(n, 1)*diag(G)' - 2*G;
            else
                [U, S, V] = svd(D);
                S(dim+3:end, dim+3:end) = 0;%Utilising rank property of EDM's
                D = U * S * V';
            end

            error1 = norm(D_prev - D, 'fro');%Calculate the error in present matrix
            
            D(W) = t_D(W);%Enforce known entries
            D(I) = 0;%Zero out the diagonal
            D(D<0) = 0;%Enforce positivity of matrix.

            error2 = norm(D_prev - D, 'fro');%Error in enforced matrix

            if (error1 < Min_err) && (error2 < Min_err)%If error is less than min error, matrix is achieved
                break;
            end
        end

        if gram < 2
            break;
        else
            gram = 1;%Gram always becomes 1 in the end
        end
    end
end