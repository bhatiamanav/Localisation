function [Si,Mj] = compute_locations(Si_cap,Mj_cap)
    Si = zeros(size(Si_cap));
    Mj = zeros(size(Mj_cap));
    %The Transformation function
    sz1 = size(Si_cap);
    sz2 = size(Mj_cap);
    O1 = ones(1,sz2(2));
    h_m = linsolve(Mj_cap.',O1.');
    I = eye(4);
    O2 = zeros(4,1);
    H_m = vertcat(h_m.',horzcat(O2,I));
    
    O3 = ones(1,sz1(2));
    h_s = linsolve(Si_cap.',O3.');
    H_s_inv = horzcat(vertcat(I,O2.'),h_s);
    H_s = inv(H_s_inv);
    Mj_cap_dash = H_s * H_m * Mj_cap; 
    Mj_cap_dash_T = Mj_cap_dash.';
    A = zeros(10,10);
    B = zeros(10,1);
    for i = 1:10
        A(i,1) = Mj_cap_dash_T(i,1).^2;
        A(i,5) = Mj_cap_dash_T(i,2).^2; 
        A(i,8) = Mj_cap_dash_T(i,3).^2;
        A(i,10) = Mj_cap_dash_T(i,4).^2;
    end
    for i = 1:10
        A(i,2) = Mj_cap_dash_T(i,1).*Mj_cap_dash_T(i,2).*2; 
        A(i,3) = Mj_cap_dash_T(i,1).*Mj_cap_dash_T(i,3).*2; 
        A(i,4) = Mj_cap_dash_T(i,1).*Mj_cap_dash_T(i,4).*2; 
        A(i,6) = Mj_cap_dash_T(i,2).*Mj_cap_dash_T(i,3).*2;
        A(i,7) = Mj_cap_dash_T(i,2).*Mj_cap_dash_T(i,4).*2;
        A(i,9) = Mj_cap_dash_T(i,3).*Mj_cap_dash_T(i,4).*2;
    end
    for i = 1:10
        B(i,1) = Mj_cap_dash_T(i,1).*Mj_cap_dash_T(i,5);
    end
    ans = linsolve(A,B);
    Q = [ans(1) ans(2) ans(3) ans(4) -1/2;
        ans(2) ans(5) ans(6) ans(7) 0;
        ans(3) ans(6) ans(8) ans(9) 0;
        ans(4) ans(7) ans(9) ans(10) 0;
        -1/2 0 0 0 0];
    Q_mid = Q(2:4,2:4);
    K = chol(abs(Q_mid));
    R = eye(3);
    t = [0 0 0];
    R1 = [1 0 0 0 0];
    R24 = horzcat(t.',R*K,zeros(3,1));
    R5 = horzcat((t*t.')-Q(1,1),2.*(t*K - [Q(1,2) Q(1,3) Q(1,4)]),ones(1,1));
    H_Q = vertcat(R1,R24,R5);
    
    H = H_Q * H_s * H_m;
    Si = (Si_cap.')*inv(H);
    Mj = H * Mj_cap;
end

