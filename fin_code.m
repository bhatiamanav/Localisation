pks = zeros(10,7);

pks(1,:) = [0 46 61 85 125 240 265];
pks(2,:) = [0 86 121 127 147 280 450];
pks(3,:) = [0 85 160 181 197 262 566];
pks(4,:) = [0 33 123 190 197 241 691];
pks(5,:) = [0 42 103 184 206 398 727];
pks(6,:) = [0 20  60 142 205 402 629];
pks(7,:) = [0 66  94 132 151 245 520];
pks(8,:) = [0 60  90  92 112 144 157];
pks(9,:) = [0 11  83 108 125 319 327];
pks(10,:) = [0 23  51  73  93 145 193];

pks = pks.';
Si = zeros(1,1);
Mj = zeros(10,10);
eps_best = Inf;

unsorted_echoes = zeros(10,10);
unsorted_echoes(1,:) = [0.00000e+00 1.11875e-02 8.93750e-03 3.24375e-02 1.19375e-02 2.36250e-02 1.86875e-02 9.25000e-03 7.50000e-04 6.25000e-04];
unsorted_echoes(2,:) = [1.06250e-03 0.00000e+00 3.93750e-03 5.31250e-03 1.16250e-02 1.46250e-02 8.18750e-03 1.37500e-02 6.50000e-03 1.01250e-02];
unsorted_echoes(3,:) = [1.03750e-02 3.12500e-03 0.00000e+00 1.53750e-02 1.48125e-02 7.56250e-03 6.25000e-05 3.93750e-03 3.06250e-03 3.75000e-04];
unsorted_echoes(4,:) = [2.25000e-03 1.71250e-02 0.00000e+00 0.00000e+00 0.00000e+00 3.62500e-03 7.81250e-03 4.04375e-02 6.81250e-03 4.32500e-02];
unsorted_echoes(5,:) = [2.21875e-02 1.02500e-02 0.00000e+00 5.62500e-04 0.00000e+00 2.37500e-03 1.01250e-02 1.58125e-02 8.12500e-04 1.21875e-02];
unsorted_echoes(6,:) = [3.62500e-03 3.75000e-04 3.75000e-04 3.81250e-03 4.50000e-03 0.00000e+00 3.12500e-04 2.18125e-02 1.25000e-03 2.68750e-02];
unsorted_echoes(7,:) = [1.81250e-03 6.25000e-05 1.12500e-03 6.37500e-03 3.50000e-03 6.43750e-03 0.00000e+00 2.97500e-02 2.56250e-03 9.93750e-03];
unsorted_echoes(8,:) = [0.00000e+00 7.75000e-03 2.46250e-02 3.50000e-03 2.75625e-02 1.46875e-02 0.00000e+00 0.00000e+00 1.32500e-02 3.87500e-03];
unsorted_echoes(9,:) = [4.93750e-03 6.25000e-05 1.05625e-02 2.21875e-02 2.76250e-02 2.03125e-02 1.22500e-02 1.68750e-03 0.00000e+00 1.87500e-04];
unsorted_echoes(10,:) = [1.53125e-02 1.80625e-02 3.48750e-02 7.37500e-03 3.78750e-02 1.63125e-02 2.65000e-02 2.81250e-03 2.14375e-02 0.00000e+00];

unsorted_echoes = unsorted_echoes .* 340;

[~, D] = alternating_descent(unsorted_echoes, 10);

for ii = 1:7

sortedEchoes = sort_echoes_local(unsorted_echoes, pks(ii,:), 3, 10);

[~, ~, microphones] = dist_opt_mex_3d(unsorted_echoes, 3);

estimatedImages  = zeros(3, 10);

for i = 1:10
    [estimatedImages(:, i), ~] = trilaterate_beck(microphones, sqrt(sortedEchoes(i, :))');
end

ti_ini = zeros(1,10);

tij = zeros(10,7);

tij(1,:) = [20  34  66  81 105 145 260 285];
tij(2,:) = [250 336 371 377 397 530 700];
tij(3,:) = [190 275 350 371 387 452 756];
tij(4,:) = [189 223 313 380 387 431 881];
tij(5,:) = [148 190 251 332 354 546 875];
tij(6,:) = [193 213 253 335 398 595 822];
tij(7,:) = [251 317 345 383 402 496 771];
tij(8,:) = [420 480 510 512 532 564 577];
tij(9,:) = [325 336 408 433 450 644 652];
tij(10,:) = [449 472 500 522 542 594 642];

tij = tij ./ 1000;

[Six, Mjx, ti, eps] = return_loc(tij);
if eps < eps_best
    Si = Six;
    Mj = Mjx;
    eps_best = eps;
end
end

disp(Mj(:,2));
disp(Mj(:,3));
disp(Mj(:,4));

d1 = zeros(len(Mj,1),len(Mj,1),3);
d2 = zeros(len(Si,1),len(Mj,1),3);
d3 = zeros(len(Mj,1),len(Si,1),3);
d4 = zeros(len(Si,1),len(Si,1),3);

for i =1 : len(Mj,1)
    for j = 1:len(Mj,1)
        d1(i,j,1) = Mj(i,2) - Mj(j,2);
        d1(i,j,2) = Mj(i,3) - Mj(j,3);
        d1(i,j,3) = Mj(i,4) - Mj(j,4);
    end
end

d1 = d1.^2;

for i =1 : len(Si,1)
    for j = 1:len(Mj,1)
        d2(i,j,1) = (-0.5*Si(i,2)) - Mj(j,2);
        d3(j,i,1) = (-0.5*Si(i,2)) - Mj(j,2);
        
        d2(i,j,2) = (-0.5*Si(i,3)) - Mj(j,3);
        d3(j,i,2) = (-0.5*Si(i,3)) - Mj(j,3);
        
        d2(i,j,3) = (-0.5*Si(i,4)) - Mj(j,4);
        d3(j,i,3) = (-0.5*Si(i,4)) - Mj(j,4);
    end
end

d2 = d2.^2;
d3 = d3.^2;

for i =1 : len(Si,1)
    for j = 1:len(Si,1)
        d4(i,j,1) = (-0.5*Si(i,2)) - (-0.5*Si(j,2));
        d4(i,j,2) = (-0.5*Si(i,3)) - (-0.5*Si(j,3));
        d4(i,j,3) = (-0.5*Si(i,4)) - (-0.5*Si(j,4));
    end
end

d4 = d4.^2;

Dx = vertcat(horzcat(d1(:,:,1),d2(:,:,1)),horzcat(d3(:,:,1),d4(:,:,1)));
Dy = vertcat(horzcat(d1(:,:,2),d2(:,:,2)),horzcat(d3(:,:,2),d4(:,:,2)));
Dz = vertcat(horzcat(d1(:,:,3),d2(:,:,3)),horzcat(d3(:,:,3),d4(:,:,3)));

[X,~] = alternating_descent(Dx,10);
[Y,~] = alternating_descent(Dy,10);
[Z,~] = alternating_descent(Dz,10);

disp(X);
disp(Y);
disp(Z);
