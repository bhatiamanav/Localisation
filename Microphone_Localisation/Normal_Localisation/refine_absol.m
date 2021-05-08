function [Si_cap,Mj_cap] = refine_absol(tij,ti,v)
    ti_dia = diag(ti);
    inp = ((tij .* tij) - 2.*(ti_dia*tij)).*(v.^2);
    [U,S,V] = svd(inp);
    disp(size(V));
    k = 5;
    Si_cap = (U(:,1:k) * sqrt(S(1:k,1:k))).';
    Mj_cap = sqrt(S(1:k,1:k)) * (V(1:k,:));
end

