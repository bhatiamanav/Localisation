n = 20; % Number of points in the set
d = 2;  % Embedding dimension

% Minimal and maximal number of deletions
n_del_min  = 0;
n_del_step = 5;
n_del_max  = 150;
n_del      = n_del_min:n_del_step:n_del_max;

% Number of random point sets per deletion level
n_random_p_sets = 500;

% Number of random masks for every point configuration
n_mask = 1;

methods = {'No Noise', 'Noise in range [-0.05,0.05]','Noise in range [-0.10,0.10]','Noise in range [-0.15,0.15]'};err = zeros(4, numel(n_del));
error = zeros(4, numel(n_del));

parfor i = 1:numel(n_del)
    err_in = zeros(4, 1);
    for l = 1:n_random_p_sets
        for j = 1:4
        
        fprintf('#(Deletions) %d in the range %d-%d, configuration %d/%d\n', n_del(i), n_del_min, n_del_max, l, n_random_p_sets);

        X = rand(d, n);      % Point set
        G = X'*X;            % Gramian
        D = edm(X, X); % EDM
        D_nois = D + ((-j+1)*0.05 + (j-1)*0.10*rand(n,n));
        
        for a = 1:n_mask
            W = random_deletion_mask(n, n_del(i));
            E = rank_complete_edm(D_nois, W, d, 0);
            err_in(j) = err_in(j) + norm(E - D, 'fro');
            err_in(j) = err_in(j) / norm(D,'fro');
        end
        end
    end
    err(:, i) = err_in;
end

err = err / n_random_p_sets / n_mask;

figure(1);
clf;

plot(n_del, err', 'LineWidth', 2);
ylabel('Relative Error');
xlabel('Number of deletions');
legend(methods, 'Location','SouthWest');
title('Relative error vs Number of deletions for Rank Complete Method');

axis tight;
grid on;