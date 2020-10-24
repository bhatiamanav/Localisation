m=20;
d=3;
k=5;
X = rand(d, m);        % Microphones
Y = rand(d, k);        % Acoustic events
D = edm(X,Y);
disp(size(D));