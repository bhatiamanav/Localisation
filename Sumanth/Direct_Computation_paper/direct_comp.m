V = 340*[0 0; 8 0; 8 8; 0 8];
s = 340*[3 2];
M = 340*[1 6; 2 5; 4 5; 4 4; 5 1; 6 6; 7 1; 7 3];
Ti = [2;2;2;2;2];

[tij,S] = compute_tij(V,s,M,Ti);
disp('True tij:');
disp(tij);

[ti] = calctod_2D_sumanth(tij);
disp('Obtained ti:');
disp(ti);

[Si_cap, Mj_cap] = refine_absol_2D_sumanth_tdoa(tij,ti);
disp('Si_cap:');
disp(Si_cap);
disp('Mj_cap:');
disp(Mj_cap);

R = [1 0; 0 1];
t = [0 0];

[Si,Mj] = compute_locations_2D_sumanth(Si_cap,Mj_cap,R,t);
disp('Si:');
disp(Si);
disp('Mj:');
disp(Mj);

[si,mj] = obtain_source_mic_locations(Si,Mj);
disp('si:');
disp(si);
disp('mj:');
disp(mj);

S = [S'; zeros(1,size(S,1))];
M = [M'; zeros(1,size(M,1))];

disp('S:');
disp(S);
disp('M:');
disp(M);

point_set_X_a = [si , mj];
point_set_Y = [S , M];

edm_X_a = zeros(size(point_set_X_a));
edm_Y = zeros(size(point_set_Y));

n1 = size(point_set_X_a,2);
n2 = size(point_set_Y,2);

[k, o] = find(triu(ones(n1),1));

disp(size(point_set_X_a));
disp(size(point_set_Y));

disp(k);
%disp(o);

edm_X_a(k + n1 * (o - 1)) = sum((point_set_X_a(:,k) - point_set_X_a(:,o)).^2, 1);
edm_X_a(o + n1 * (k - 1)) = edm_X_a(k + n1 * (o - 1));

[r, u] = find(triu(ones(n2),1));
edm_Y(r + n2 * (u - 1)) = sum((point_set_Y(:,r) - point_set_Y(:,u)).^2, 1);
edm_Y(u + n2 * (r - 1)) = edm_Y(r + n * (u - 1));

disp('edm_X_a:');
disp(edm_X_a);
disp('edm_Y:');
disp(edm_Y);