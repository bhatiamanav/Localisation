n = 20; % Number of points in the set
d = 2;  % Embedding dimension

% Minimal and maximal number of deletions
n_del_min  = 0;
n_del_step = 5;
n_del_max  = 150;
n_del      = n_del_min:n_del_step:n_del_max;

% Number of random point sets per deletion level
n_rand_p_sets = 500;

% Number of random masks for every point configuration
n_mask = 1;

methods = {'Alternating Descent', 'Rank Alternation','Semidefinite Relaxation'};
success = zeros(3, numel(n_del));%numel calculates the number of elements in a array
MIN_ERR = 1e-2;%Criteria to measure success

parfor i = 1:numel(n_del)%Parfor runs the algo on parallel cores to increase its speed for otherwise bulky algorithms
    success_in = zeros(3, 1);
    for j = 1:n_rand_p_sets
        
        fprintf('#(Deletions) %d in the range %d-%d, configuration %d/%d\n', n_del(i), n_del_min, n_del_max, j, n_rand_p_sets);

        X = rand(d, n);%Random Point set      
        G = X'*X; %Gramanian Matrix           
        D = edm(X, X);%Euclidean Distance Matrix
        
        %Success is 1 if ||E-D||f/||D||f is less than minimum error
        for k = 1:n_mask
            W = random_deletion_mask(n, n_del(i));
            
            [~, E] = alternating_descent(D .* W, d);
            success_in(1) = success_in(1) + (norm(E - D, 'fro') < MIN_ERR*norm(D, 'fro'));
         
            E = rank_complete_edm(D, W, d, 0);
            success_in(2) = success_in(2) + (norm(E - D, 'fro') < MIN_ERR*norm(D, 'fro'));
            
            E = sdr_complete_edm_noise(D, W, d);
            success_in(3) = success_in(3) + (norm(E - D, 'fro') < MIN_ERR*norm(D, 'fro'));
        end
    end
    success(:, i) = success_in;
end

success = success / n_rand_p_sets / n_mask;
figure(1);
clf;

plot(n_del, success', 'LineWidth', 2);
ylabel('Success percentage');
xlabel('Number of deletions');
legend(methods, 'Location','SouthWest');

axis tight;
grid on;