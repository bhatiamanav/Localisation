cities = {'Lausanne','Geneva','Zurich','Neuchatel','Bern'};%Tags of cities
      
D = [ 0   33  128   40   66     
      0    0  158   64  101    
      0    0    0   88   56     
      0    0    0    0   34     
      0    0    0    0    0   ];%The Distance Matrix
  
 n = size(D,1);
 D = D+D';
 D = D.^2;%The Euclidean distance matrix
 
 X = classic_mds(D,2);%classical Multidimensional scaling with dim=2
 
 X(1, :) = -X(1, :);%Translation
 
phi = 130 * pi/180;
R = [ cos(phi) sin(phi)
     -sin(phi) cos(phi) ];%Rotation

X = R*X;

%Plotting

figure(1);
clf;
hold on;
for j = 1:n
    scatter(X(1, j), X(2, j), 'filled');
    text(X(1, j), X(2, j), cities{j});
end
grid on;
grid minor;
axis equal;
axis tight;