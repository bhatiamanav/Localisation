function D = sdr_complete_edm_noise(t_D, W, dim, l)
    n = size(t_D, 1);
    x = -1/(n + sqrt(n));
    y = -1/sqrt(n);
    V = [y*ones(1, n-1); x*ones(n-1) + eye(n-1)];
    e = ones(n, 1);

    cvx_quiet('true');
    cvx_begin sdp
        variable G(n-1, n-1) symmetric;
        B = V*G*V';
        E = diag(B)*e' + e*diag(B)' - 2*B;
        maximize (trace(G) - l * norm(W.*(E -  t_D), 'fro'));
        subject to
            G >= 0;
    cvx_end
    cvx_quiet('false');

    B = V*G*V';
    D = diag(B)*e' + e*diag(B)' - 2*B;
end
%A similar CVX code to previous with the exception that the equality of Hadamard product is 
%removed by maximising the new quantity of frobenius norm subtracted from
%the trace
%Uses Algorithm 5 of the Dokmanic paper