function [X, D] = alternating_descent(t_D, dim)
    N = 50;%Max. number of iterations
    [n, m] = size(t_D);%Extracting the matrix size
    if (n ~= m)%Does not work on rectangular matrices
        error('The inpur matrix t_D needs to be a square matrix!')
    end

    L = eye(n) - 1/n*ones(n,1)*ones(n,1).';%L = I - 1/n*e*e'

    X0 = zeros(n, dim);
    X = X0;

    %Extract the information related to connectivity
    conn_nodes = cell(1, n);%Connected Nodes
    conn_vector = cell(1, n);%Connectivity vector
    sens_list = 1:n;%List of sensors from 1 to n

    for i = 1 : n
        sensor_ind = sens_list(i);
        conn_nodes{i} = find(t_D(sensor_ind, :) ~= 0);%Finding the sensors with known distances
        conn_vector{i} = t_D(sensor_ind, conn_nodes{i});%Feeding known information to connectivity vector
        conn_vector{i} = conn_vector{i}(:);
    end

    for iter_ind = 1 : N%Repeating for all N points
        for i = 1 : n%Repeating for all k coordinates

            X_conn = X(conn_nodes{i}, :);
            for cor_ind = 1 : dim%Reapeating for all dimensions
                %Finding derivative of polynomial number (26) of the Dokmanic paper
                a = 4 * size(X_conn, 1);
                b = 12 * sum((X(i,cor_ind) -  X_conn(:,cor_ind)));%Creating polynomial coefficients for the equation
                c = 4 * sum(sum((repmat(X(i,:), size(X_conn, 1), 1) - X_conn).^2, 2) + 2 * (X(i,cor_ind) -  X_conn(:,cor_ind)).^2 - conn_vector{i});
                d = 4 * sum( (X(i,cor_ind) -  X_conn(:,cor_ind)) .* (sum((repmat(X(i,:), size(X_conn, 1), 1)- X_conn).^2, 2) - conn_vector{i}));

                roots = cubicfcnroots(a,b,c,d);%Finding the roots for the equation
                roots = real(roots(abs(imag(roots)) < 1e-15));%Condition for convergence

                if isempty(roots)%If no roots, skip the chance
                    continue;
                end

                cc = zeros(1, length(roots));
                %Finding the values of derivatives at all roots
                for j = 1 : length(roots)
                    cc(j) = sum(((X(i,cor_ind)+roots(j)-X_conn(:,cor_ind)).^2 + sum((repmat(X(i,:), size(X_conn, 1), 1) - X_conn).^2, 2) - (X(i, cor_ind) - X_conn(:, cor_ind)).^2 - conn_vector{i}).^2);
                end
                [~, cc_min_ind] = min(cc);%Finding the global minimizer
                min_value = roots(cc_min_ind);
                X(i,cor_ind) = X(i,cor_ind) + min_value;%Replacing the value for the global minimizer

                X = L * X;
            end
        end
    end
    
    %Constructing the EDM
    X = L * X;
    [k, o] = find(triu(ones(n),1));%Diagonal remains zero

    D = zeros(n);
    D(k + n * (o - 1)) = sum((X(k,:) - X(o,:)).^2, 2);
    D(o + n * (k - 1)) = D(k + n * (o - 1));%Enforcing symmetry about the diagonal
end
