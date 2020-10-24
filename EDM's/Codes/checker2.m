m = 20; % Number of microphones
d = 3;  % Embedding dimension

% Minimal and maximal number of acoustic events
n_events_min =  5;
n_events_step = 1;
n_events_max = 30;
n_events = n_events_min:n_events_step:n_events_max;

% Number of random point sets per event number
n_rand_p_sets = 100;

methods = {'A Method by Crocco', 'Alternating Descent', 'Rank alternation', 'Semidefinite Relaxation'};
success = zeros(4, numel(n_events));
MIN_ERR = 1e-2;

parfor i = 1:numel(n_events)
    
    k = n_events(i);
    
    err_in = zeros(4, 1);
    success_in = zeros(4, 1);
    for j = 1:n_rand_p_sets
        
        fprintf('#(Events) %d in the range %d-%d, configuration %d/%d\n', k, n_events_min, n_events_max, j, n_rand_p_sets);

        X = rand(d, m);        % Microphones
        Y = rand(d, k);        % Acoustic events
        D = edm([X Y], [X Y]); % EDM
        W = mdu_mask(m, k);%Mask for Acoustic Methods

        [~, ~, E] = crocco_method(sqrt(edm(X, Y)), 1);
        success_in(1) = success_in(1) + (norm(E(1:m, 1:m) - D(1:m, 1:m), 'fro') < MIN_ERR*norm(D(1:m, 1:m), 'fro'));

        [~, E] = alternating_descent(D .* W, d);
        success_in(2) = success_in(2) + (norm(E(1:m, 1:m) - D(1:m, 1:m), 'fro') < MIN_ERR*norm(D(1:m, 1:m), 'fro'));

        E = rank_complete_edm(D, W, d, 0);
        success_in(3) = success_in(3) + (norm(E(1:m, 1:m) - D(1:m, 1:m), 'fro') < MIN_ERR*norm(D(1:m, 1:m), 'fro'));

        E = sdr_complete_edm_noise(D, W, d);
        success_in(4) = success_in(4) + (norm(E(1:m, 1:m) - D(1:m, 1:m), 'fro') < MIN_ERR*norm(D(1:m, 1:m), 'fro'));
    end
    success(:, i) = success_in;
end

success = success / n_rand_p_sets;

figure(1);

plot(n_events, success', 'LineWidth', 2);
ylabel('Success percentage');
xlabel('Number of acoustic events');
legend(methods, 'Location','East');

axis tight;
grid on;