function D = edm(X, Y)
    dist_X2 = sum(X.^2); %Norm of X
    dist_Y2 = sum(Y.^2); %Norm of Y

    D = bsxfun(@plus, dist_X2', dist_Y2) - 2*X'*Y;%Using the formula for an EDM
end