function [ti,s] = ctod_2D(tij)
    sz = size(tij);
    ti = zeros(1,sz(1));

    A = tij .^ 2;
    B = -2 * tij;
    O = ones(1,sz(2));
    n = 4;
    s1 = size(tij, 1);
    for i = 1:s1/n
        A_bar = [A((i-1)*n+1,:); A((i-1)*n+2,:); A((i-1)*n+3,:); A((i-1)*n+4,:)];
        B_bar = [B((i-1)*n+1,:); B((i-1)*n+2,:); B((i-1)*n+3,:); B((i-1)*n+4,:)];
        Y =[A_bar; B_bar];

        s = O * pinv(Y);
        ti((i-1)*n + 1) = s(5)/s(1);
        ti((i-1)*n + 2) = s(6)/s(2);
        ti((i-1)*n + 3) = s(7)/s(3);
        ti((i-1)*n + 4) = s(8)/s(4);
    end
end