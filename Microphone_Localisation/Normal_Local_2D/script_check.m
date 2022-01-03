clc;
clear;

V = [0 0; 8 0; 8 8; 0 8] * 340;
s= [3 2] * 340;
M = [1 6; 2 5; 4 5; 4 4; 5 1; 6 6] * 340;

Ti = [2;2;2;2;2];

tij = compute_tij(V,s,M,Ti);
ti_est = ctod_2D(tij);

disp(tij);
disp(ti_est);