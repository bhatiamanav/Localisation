function W = random_deletion_mask(n, n_del)
    k = randperm(n*(n - 1)/2, n_del);
    W = zeros(n);
    i = find(tril(ones(n), -1));
    W(i(k)) = 1;
    W = logical(ones(n) - (W + W'));
end
%Given some number of points and some distances to remove, it generates
%a mask to randomly generate those n_del points and delete them, that is,
%set them to zero and known distances to 1 in the masking matrix.It chooses
%the randomly to be removen points uniformly in the matrix.