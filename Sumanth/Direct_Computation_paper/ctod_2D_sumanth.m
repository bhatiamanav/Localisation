function [ti,X] = ctod_2D_sumanth(tij)
    sz = size(tij);
    ti = zeros(1,sz(1));
    
    A = tij .^ 2;
    B = -1 * 2 * tij;
    O = ones(1,sz(2));
    n = 4;             
    s1 = size(tij, 1);      
    for i = 1:s1/n
            A_bar = [A((i-1)*n+1,:); A((i-1)*n+2,:); A((i-1)*n+3,:); A((i-1)*n+4,:)];
            B_bar = [B((i-1)*n+1,:); B((i-1)*n+2,:); B((i-1)*n+3,:); B((i-1)*n+4,:)];
            Y =[A_bar; B_bar];
            %Y = vertcat(A_bar, B_bar);
            disp(rank(Y));
            %X = O * pinv(Y);
            X = linsolve(Y',O');
            
            ti((i-1)*n + 1) = X(5)/X(1);
            ti((i-1)*n + 2) = X(6)/X(2);
            ti((i-1)*n + 3) = X(7)/X(3);
            ti((i-1)*n + 4) = X(8)/X(4);
    end
end