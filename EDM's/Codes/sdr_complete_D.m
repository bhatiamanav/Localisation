function D = sdr_complete_D(t_D, W, dim)
    n = size(t_D, 1);%Size of the EDM
    
    x = -1/(n + sqrt(n));
    y = -1/sqrt(n);%Initialising the given variables
    V = [y*ones(1, n-1); x*ones(n-1) + eye(n-1)];
    e = ones(n, 1);

    cvx_quiet('true');
    cvx_begin sdp
        variable G(n-1, n-1) symmetric;
        maximize trace(G);
        B = V*G*V';
        E = diag(B)*e' + e*diag(B)' - 2*B;
        subject to
            E .* W == t_D .* W;
            G >= 0;
    cvx_end
    cvx_quiet('false');

    B = V*G*V';
    D = diag(B)*e' + e*diag(B)' - 2*B;
end
