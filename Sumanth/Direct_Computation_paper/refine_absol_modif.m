function [Si_cap, Mj_cap] = refine_absol_modif(tij,ti,v)  
    n = size(tij,2);
    v_sq = (v * v);
    ones_1_n = ones(1,n);
    prop_mat = tij - ti'*ones_1_n;
    prop_mat_sq = (prop_mat).*(prop_mat);
    inp_v = v_sq*prop_mat_sq;
    
    [U,S,V] = svd(inp_v);
    
    k=4;
    Si_cap = (U(:,1:k) * sqrt(S(1:k,1:k))).';
    Mj_cap = (V(:,1:k) * sqrt(S(1:k,1:k))).';
end