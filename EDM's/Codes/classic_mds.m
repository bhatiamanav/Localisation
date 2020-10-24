function X = classic_mds(D, dim)

    n = size(D, 1);
    I = eye(n);
    J = I - 1/n*ones(n);

    [U, S, V] = svd(-J*D*J/2);
    S = S(1:dim, :);
    X = sqrt(S)*V';
end
%Classical Multi-dimensional scaling.
%Given a EDM and an embedding dimension, it returns a list of the centered
%points for the given EDM.