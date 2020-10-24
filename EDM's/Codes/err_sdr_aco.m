m = 20; % Number of microphones
d = 3;  % Embedding dimension

% Minimal and maximal number of acoustic events
n_events_min =  5;
n_events_step = 1;
n_events_max = 30;
n_events = n_events_min:n_events_step:n_events_max;

% Number of random point sets per event number
n_random_p_sets = 100;

methods = {'No Noise', 'Noise in range [-0.05,0.05]','Noise in range [-0.10,0.10]','Noise in range [-0.15,0.15]'};
err = zeros(4, numel(n_events));
parfor i = 1:numel(n_events)
    k = n_events(i);
    err_in = zeros(4, 1);
    for l = 1:n_random_p_sets
        for j=1:4
        
        fprintf('#(Events) %d in the range %d-%d, configuration %d/%d\n', k, n_events_min, n_events_max, l, n_random_p_sets);

        X = rand(d, m);        % Microphones
        Y = rand(d, k);        % Acoustic events
        D = edm([X Y], [X Y]); % EDM
        D_nois = D + ((-j+1)*0.05 + (j-1)*0.10*rand(m+k,m+k));
        W = mdu_mask(m, k);

       
        E = sdr_complete_edm_noise(D_nois, W, d,1);
        err_in(j) = err_in(j) + norm(E(1:m, 1:m) - D(1:m, 1:m), 'fro');
        err_in(j) = err_in(j)/(norm(D(1:m,1:m),'fro'));
        end
    end
    err(:, i) = err_in;
end

err = err / n_random_p_sets;

figure(1);

plot(n_events, err', 'LineWidth', 2);
ylabel('Relative Error');
xlabel('Number of acoustic events');
legend(methods, 'Location','East');
title('Relative Error vs number of Acoustic Events for SDR completion');

axis tight;
grid on;