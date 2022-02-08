function [ti] = calctod_2D_sumanth(tij)

    sz = size(tij);
    %ti = zeros(1,sz(1));

    A = tij .^ 2;
    B = -2 * tij;

    O = ones(1,sz(2));

    n1 = 4;
    s1 = size(tij, 1);

    comb = combntns(1:s1,n1);
    %disp(sz);
    %disp(comb);
    sz_comb = size(comb,1);
    Ti_vec = zeros(sz_comb,s1);
    tin = zeros(1,sz(1));

    for i = 1:sz_comb
    j = comb(i,:);
    A_bar = zeros(n1,sz(2));
    B_bar = zeros(n1,sz(2));
    tii = zeros(1,s1);
    for k = 1:n1
    A_bar(k,:) = A(j(k),:);
    B_bar(k,:) = B(j(k),:);
    end

    Y =[A_bar; B_bar];

    X = O * pinv(Y);

    tin(1) = X(5)/X(1);
    tin(2) = X(6)/X(2);
    tin(3) = X(7)/X(3);
    tin(4) = X(8)/X(4);
    %end
    for l = 1:n1
    %disp(j(l));
    tii(:,j(l)) = tin(l);
    end
    Ti_vec(i,:) = tii;
    end

    Ti_int = cast(Ti_vec,'uint8');
    ti_mode = mode(Ti_int);
    ti = double(ti_mode);

end