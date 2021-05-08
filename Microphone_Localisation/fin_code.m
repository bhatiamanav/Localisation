[y, fs] = audioread('notebooks_arctic_a0010.wav');
rir_small = load('rir.mat');

mult_sig = rir_small .* y;
acf = autocorr(mult_sig);
pks = findpeaks(acf);

unsorted_echoes_x = (bsxfun(@minus, x,x)).^2;
unsorted_echoes_y = (bsxfun(@minus,y,y)).^2;
unsorted_echoes_z = (bsxfun(@minus,z,z)).^2;

unsorted_echoes = sqrt(unsorted_echoes_x + unsorted_echoes_y + unsorted_echoes_z);

[echoCombinations, scores] = sort_echoes(unsorted_echoes, pks, 3, 6);

[~, ~, microphones] = dist_opt_mex_3d(D, 3);

estimatedImages  = zeros(3, 6);

for i = 1:6
    [estimatedImages(:, i), ~] = trilaterate_beck(microphones, sqrt(sortedEchoes(i, :))');
end

ti_ini = zeros(1,10);

eps_best = Inf;

[Si, Mj, ti, eps] = return_loc(ti_ini,estimatedImages(1,:), estimatedImages(2,:), estimatedImages(3,:), x, y, z);
print(Mj(:,2));
print(Mj(:,3));
print(Mj(:,4));

d1 = zeros(len(Mj,1),len(Mj,1),3);
d2 = zeros(len(Si,1),len(Mj,1),3);
d3 = zeros(len(Mj,1),len(Si,1),3);
d4 = zeros(len(Si,1),len(Si,1),3);

for i =1 : len(Mj,1)
    for j = 1:len(Mj,1)
        d1(i,j,1) = sqrt(Mj(i,2),Mj(j,2));
        d1(i,j,2) = sqrt(Mj(i,3),Mj(j,3));
        d1(i,j,3) = sqrt(Mj(i,4),Mj(j,4));
    end
end

for i =1 : len(Si,1)
    for j = 1:len(Mj,1)
        d2(i,j,1) = sqrt(-0.5*Si(i,2),Mj(j,2));
        d3(i,j,1) = sqrt(-0.5*Si(i,2),Mj(j,2));
        
        d2(i,j,2) = sqrt(-0.5*Si(i,3),Mj(j,3));
        d3(i,j,2) = sqrt(-0.5*Si(i,3),Mj(j,3));
        
        d2(i,j,3) = sqrt(-0.5*Si(i,4),Mj(j,4));
        d3(i,j,3) = sqrt(-0.5*Si(i,4),Mj(j,4));
    end
end


for i =1 : len(Si,1)
    for j = 1:len(Si,1)
        d4(i,j,1) = sqrt(-0.5*Si(i,2),-0.5*Si(j,2));
        d4(i,j,2) = sqrt(-0.5*Si(i,3),-0.5*Si(j,3));
        d4(i,j,3) = sqrt(-0.5*Si(i,4),-0.5*Si(j,4));
    end
end

Dx = vertcat(horzcat(d1(:,:,1),d2(:,:,1)),horzcat(d3(:,:,1),d4(:,:,1)));
Dy = vertcat(horzcat(d1(:,:,2),d2(:,:,2)),horzcat(d3(:,:,2),d4(:,:,2)));
Dz = vertcat(horzcat(d1(:,:,3),d2(:,:,3)),horzcat(d3(:,:,3),d4(:,:,3)));

[X,~] = alternating_descent(Dx,10);
[Y,~] = alternating_descent(Dy,10);
[Z,~] = alternating_descent(Dz,10);

print(X);
print(Y);
print(Z);
