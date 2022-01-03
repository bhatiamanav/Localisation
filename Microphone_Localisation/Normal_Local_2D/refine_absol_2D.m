function [Si_quote_cap, Si_cap, Mj_cap] = refine_absol_2D(tij,ti,v)
    ti_dia = diag(ti);

    sz = size(tij);
    a_mat = sz(1);

    ti_sq = ti.*ti;
    v_sq = (v * v);
    b_mat = v_sq*ti_sq;

    inp_2 = (tij .^ 2) - 2.*(ti_dia*tij);
    inp =zeros(sz(1),sz(2));
    disp(inp - inp_2);
    for ii = 1:sz(1)
        for jj = 1:sz(2)
            inp(ii,jj) = (tij(ii,jj) * tij(ii,jj)) - (2 * tij(ii,jj) * ti_dia(ii,ii)); 
        end
    end
    inp_v = v_sq .* inp;

    [U,S,V] = svd(inp_v);

    k = 4;
    c = zeros(k-1,a_mat);
    A = [b_mat;c];

    Si_quote_cap = (U(:,1:k) * sqrt(S(1:k,1:k))).';
    Si_cap = Si_quote_cap + A;
    Mj_cap = (V(:,1:k) * sqrt(S(1:k,1:k))).';
end