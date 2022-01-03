function [Si,Mj] = compute_locations_2D(Si_cap,Mj_cap)
    %The Transformation function
    sz1 = size(Si_cap);
    sz2 = size(Mj_cap);
    O1 = ones(1,sz2(2));
    h_m = linsolve(Mj_cap',O1');
    I = eye(3);
    O2 = zeros(3,1);
    O22 = [O2, I];
    H_m = [h_m.'; O22];

    O3 = ones(1,sz1(2));
    h_s = linsolve(Si_cap.',O3.');
    O33 = [I ; O2.'];
    H_s_inv = [O33 , h_s];
    %H_s = inv(H_s_inv);
    
    Mj_cap_dash = (H_s_inv\H_m) * Mj_cap;
    Mj_cap_dash_T = Mj_cap_dash';
    w = sz2(2);
    C = zeros(w,6);
    D = zeros(w,1);

    for i = 1:w
        C(i,1) = Mj_cap_dash_T(i,1).^2;
        C(i,2) = Mj_cap_dash_T(i,1)*Mj_cap_dash_T(i,2)*2;
        C(i,3) = Mj_cap_dash_T(i,1)*Mj_cap_dash_T(i,3)*2;
        C(i,4) = Mj_cap_dash_T(i,2).^2;
        C(i,5) = Mj_cap_dash_T(i,2)*Mj_cap_dash_T(i,3)*2;
        C(i,6) = Mj_cap_dash_T(i,3).^2;
        D(i,1) = Mj_cap_dash_T(i,1)*Mj_cap_dash_T(i,4);
    end

    E = pinv(C) * D;
    %E = lsqr(C,D);
    %disp(E);
    
    Q = [
        E(1) E(2) E(3) -1/2; 
        E(2) E(4) E(5) 0; 
        E(3) E(5) E(6) 0;
        -1/2 0 0 0
        ];

    Q_mid = Q(2:3,2:3);
    disp(Q_mid);
    %[V,D] = eig(Q_mid);
    %disp(D);
    %K_abs = V * abs(D) * V.';
    %[~,D] = eig(K_abs);
    %disp(D);
    %K = chol(K_abs, 'lower');
    [L,D,P] = ldl(Q_mid);
    K = ((P.')\L) * sqrt(D);
    disp(abs(K));
    
    R = eye(2);
    t = [0 0];
    
    R1 = [1 0 0 0];
    R23 = [t',R*abs(K),zeros(2,1)];
    R4 = [(t*t')-Q(1,1),2*(t*R*abs(K) - [Q(1,2) Q(1,3)]),ones(1,1)];
    H_Q = [R1;R23;R4];

    H = H_Q * (H_s_inv\H_m);
    H_inv = inv(H);
    Si = (H_inv.')*Si_cap;
    Mj = H * Mj_cap;
    Si = real(Si);
    Mj = real(Mj);
end