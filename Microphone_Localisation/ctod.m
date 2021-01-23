function [ti,ans] = ctod(tij)
    sz = size(tij);
    ti = zeros(1,sz(1));
    mnj = zeros(1,sz(1));
    for i = 1:sz(1)
        for j = 1:sz(2)
            mnj(i) = mnj(i) + tij(i,j);
        end
    end
    mnj = mnj / sz(1);
    for i = 1:sz(1)
        for j = 1:sz(2)
            tij(i,j) = tij(i,j) - mnj(i);
        end
    end
    A = tij .^ 2;
    B = -1 * 2 * tij;
    O = ones(1,sz(2));
    n = 5;             
    s1 = size(tij, 1);      
    for i = 1:s1/n
            %A_bar = horzcat(A(:,(i-1)*n+1),A(:,(i-1)*n+2),A(:,(i-1)*n+3),A(:,(i-1)*n+4),A(:,(i-1)*n+5));
            %B_bar = horzcat(B(:,(i-1)*n+1),B(:,(i-1)*n+2),B(:,(i-1)*n+3),B(:,(i-1)*n+4),B(:,(i-1)*n+5));
            A_bar = horzcat(A((i-1)*n+1,:).',A((i-1)*n+2,:).',A((i-1)*n+3,:).',A((i-1)*n+4,:).',A((i-1)*n+5,:).').';
            B_bar = horzcat(B((i-1)*n+1,:).',B((i-1)*n+2,:).',B((i-1)*n+3,:).',B((i-1)*n+4,:).',B((i-1)*n+5,:).').';
            Y = vertcat(A_bar,B_bar);
            disp(size(Y));
            sq = Y * Y.';
            ans = O * Y.' * inv(sq);
            disp(ans);
            ti((i-1)*n + 1) = ans(6)/ans(1);
            ti((i-1)*n + 2) = ans(7)/ans(2);
            ti((i-1)*n + 3) = ans(8)/ans(3);
            ti((i-1)*n + 4) = ans(9)/ans(4);
            ti((i-1)*n + 5) = ans(10)/ans(5);
    end

end



