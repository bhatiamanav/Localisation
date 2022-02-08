clear;
clc;

m = 4;
n = 8;

[V,s,M,Ti] = random_room_generator_2d(m,n);

tij = compute_tij(V,s,M,Ti);
v = 340;

[Si_cap, Mj_cap] = refine_absol_modif(tij,Ti,v);
[Si,Mj] = compute_locations_2D_sumanth(Si_cap,Mj_cap);
%disp(M);
%disp(Mj(2:3,:));
   