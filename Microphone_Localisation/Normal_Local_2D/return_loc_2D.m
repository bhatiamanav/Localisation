function [SX_hat,SY_hat,SZ_hat,MX_hat,MY_hat, ti, eps] = return_loc_2D(tij)
    v = 340;
    
    V = v*tij;
    [ti, ~] = ctod_2D(tij);
    disp(ti);
    [Si_cap,Mj_cap] = refine_absol_2D(tij,ti,v);
    [Si,Mj] = compute_locations_2D(Si_cap,Mj_cap);

    disp(size(Si));
    disp(size(Mj));

    eps = 0;
    SX_hat = zeros(1,size(Si,2));
    SY_hat = zeros(1,size(Si,2));
    SZ_hat = zeros(1,size(Si,2));

    MX_hat = zeros(1,size(Mj,2));
    MY_hat = zeros(1,size(Mj,2));

    for i = 1:size(Si,2)
        SX_hat(i) = (-0.5) * Si(2,i);
        SY_hat(i) = (-0.5) * Si(3,i);
        SZ_hat(i) = sqrt(si(1,i)- 0.25*si(2,i)*si(2,i)-0.25*si(3,i)*si(3,i));
    end

    for i = 1:size(Mj,2)
        MX_hat(i) = Mj(2,i);
        MY_hat(i) = Mj(3,i);
    end
    
    for i = 1:size(Si,2)
        for j = 1:size(Mj,2)
            eps = eps + (V(i,j) - (ti(i)*v + distance_2D(SX_hat(i),SY_hat(i),SZ_hat(i),MX_hat(j),MY_hat(j)))).^2;
        end
    end
end