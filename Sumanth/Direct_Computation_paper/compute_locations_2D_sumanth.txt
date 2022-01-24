function [Si,Mj] = compute_locations_2D_sumanth(Si_cap,Mj_cap,R,t)

    %The Transformation function
    sz1 = size(Si_cap);
    sz2 = size(Mj_cap);
    O1 = ones(1,sz2(2));
    h_m = linsolve(Mj_cap',O1');
    I = eye(3);
    O2 = zeros(3,1);
    O22 = [O2, I];
    H_m = [h_m'; O22];
    
    O3 = ones(1,sz1(2));
    h_s_sub = linsolve(Si_cap(2:sz1(1),:)',O3');
    h_s = [0; h_s_sub];
    O33 = [I ; O2'];
    H_s_inv = [O33 , h_s];
    H_s = inv(H_s_inv);

    Mj_cap_dash = (H_s_inv\H_m) * Mj_cap; 
    Mj_cap_dash_T = Mj_cap_dash';
    w = sz2(2);
    C = zeros(w,6);
    D = zeros(w,1);
    
    for i = 1:w
        C(i,1) = (Mj_cap_dash_T(i,1))^2;
        C(i,2) = 2*(Mj_cap_dash_T(i,1)*Mj_cap_dash_T(i,2)) ;
        C(i,3) = 2*(Mj_cap_dash_T(i,1)*Mj_cap_dash_T(i,3));
        C(i,4) = (Mj_cap_dash_T(i,2))^2;
        C(i,5) = 2*(Mj_cap_dash_T(i,2)*Mj_cap_dash_T(i,3));
        C(i,6) = (Mj_cap_dash_T(i,3))^2;
        D(i,1) = Mj_cap_dash_T(i,1)*Mj_cap_dash_T(i,4);
    end
    
    E = linsolve(C,D);
    %E = pinv(C) * D;
    %E = lsqr(C,D);
    %disp(E);
    
    %E = linsolve(C,D);
    Q = [E(1) E(2) E(3) -1/2; E(2) E(4) E(5) 0; E(3) E(5) E(6) 0;  -1/2 0 0 0];
    disp('Q:');
    disp(Q);

    Q_mid = Q(2:3,2:3);
    disp('Q_mid:');
    disp(Q_mid);
     
    [L,DD] = ldl(Q_mid);
     K = (L*sqrt(DD))';
     disp('K:');
     disp(K);
     Q_mid_1 = conj(K')*K;
     disp('Q_mid_1:');
    disp(Q_mid_1);
     
     %K1 = chol(Q_mid);
      %disp('K1:');
     %disp(K1);
    %Q_mid_2 = (K1')*K1;
    %disp('Q_mid_2:');
    %disp(Q_mid_2);
     
    %R = eye(2);
    %t = [0 0];
    R1 = [1 0 0 0];
    R23 = [t',R*K,zeros(2,1)];
    R4 = [(t*t')-Q(1,1),2*(t*R*K - [Q(1,2) Q(1,3)]),ones(1,1)];
    H_Q = [R1;R23;R4];
    
     disp('H_Q:');
    disp(H_Q);
    disp('H_s:');
    disp(H_s);
    disp('H_m:');
    disp(H_m);
    
    H = H_Q * (H_s_inv\H_m);
    disp('H:');
    disp(H);
    H_inv = inv(H);
     disp('H_inv:');
    disp(H_inv);
    H_inv_x_H = H_inv*H;
    disp('H_inv_x_H');
    disp(H_inv_x_H);
    Si = H_inv'*Si_cap;
    Mj = H * Mj_cap;
    Si = real(Si);
    Mj = real(Mj);
    
    disp('Si_cap:');
    disp(Si_cap);
    disp('Mj_cap:');
    disp(Mj_cap);
    
end