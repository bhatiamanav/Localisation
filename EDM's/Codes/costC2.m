function C = costC2(C, U, Sigma, V, D, a1)
X_tilde = (U*C)';
Y_tilde = -1/2*inv(C)*Sigma*V';
X = [ [0 0 0]' X_tilde ];
Y = [ [0 0 0]' Y_tilde ];
Y(1, :) = Y(1, :) + a1;
C = norm(edm(X, Y) - D, 'fro')^2;
end