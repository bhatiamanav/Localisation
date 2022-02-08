function [Si_cap, Mj_cap] = refine_absol_2D_sumanth_tdoa(tij,ti)
    v = 340;
    n = size(tij,2);
    v_sq = (v * v);
    ones_1_n = ones(1,n);
    prop_mat = tij - ti'*ones_1_n;
    prop_mat_sq = (prop_mat).*(prop_mat);
    inp_v = v_sq*prop_mat_sq;
    
    [U,S,V] = svd(inp_v);
      
    k=4;
    S_Trans = S';

    Si_cap = (U(:,:) * sqrt(S(:,1:k))).';
    Mj_cap = (V(:,:) * sqrt(S_Trans(:,1:k))).';
end