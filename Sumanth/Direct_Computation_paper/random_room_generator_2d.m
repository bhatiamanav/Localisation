function [V,s,M,Ti] = random_room_generator_2d(m,n)
    rng(1, 'philox');
    sr = rng;
    disp(sr);
    r1 = randi([5, 20], 1 , 2);
    r2 = randi([-10, 10], 1 , 2);
    trans = randi([0, 10], 2 , 1);
    ang = randi([0, 360], 1 , 1);
    Ti = randi([1, 10], 1 , m+4);

    side1 = r1(1);
    side2 = r1(2);

    %disp('side1 length:');
    %disp(side1);
    %disp('side2 length:');
    %disp(side2);
    %disp('First vertex of the room:');
    %disp(r2);

    s_ini = zeros(2,m);
    %M_ini = zeros(2,n);
    V_ini = zeros(2,4);
    all_matrix = zeros(2,m+n);

    V_ini(1,1) = r2(1);
    V_ini(2,1) = r2(2);

    V_ini(1,2) = r2(1) + side1;
    V_ini(2,2) = r2(2);

    V_ini(1,3) = r2(1) + side1;
    V_ini(2,3) = r2(2) + side2;

    V_ini(1,4) = r2(1);
    V_ini(2,4) = r2(2) + side2;

    Rot = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];
    R3 = [trans(1) 0; 0 trans(2)];
    ones_2_4 = ones(2,4);
    ones_2_m = ones(2,m);
    ones_2_n = ones(2,n);

    V = 340*(Rot*V_ini + R3*ones_2_4)';

    x1 = randi([r2(1) + 1, r2(1) + side1 - 1]);
    y1 = randi([r2(2) + 1, r2(2) + side2 - 1]);
    s_ini(1,1) = x1;
    s_ini(2,1) = y1;

    if m>=2
        for i = 2:m
            for k = 1:i-1
                x1 = randi([r2(1) + 1, r2(1) + side1 - 1]);
                y1 = randi([r2(2) + 1, r2(2) + side2 - 1]);
                s_ini(1,i) = x1;
                s_ini(2,i) = y1;
                while 1
                    if s_ini(:,i) == s_ini(:,k)
                        x1 = randi([r2(1) + 1, r2(1) + side1 - 1]);
                        y1 = randi([r2(2) + 1, r2(2) + side2 - 1]);
                        s_ini(1,i) = x1;
                        s_ini(2,i) = y1;
                    else
                        break;  
                    end
                end
            end
        end
    end

    all_matrix(:,1:m) = s_ini;

    for j = m+1:m+n
        for l = 1:j-1
            x1 = randi([r2(1) + 1, r2(1) + side1 - 1]);
            y1 = randi([r2(2) + 1, r2(2) + side2 - 1]);
            all_matrix(1,j) = x1;
            all_matrix(2,j) = y1;
            while 1
                if all_matrix(:,j) == all_matrix(:,l)
                    x1 = randi([r2(1) + 1, r2(1) + side1 - 1]);
                    y1 = randi([r2(2) + 1, r2(2) + side2 - 1]);
                    all_matrix(1,j) = x1;
                    all_matrix(2,j) = y1;
                else
                    break;  
                end
            end
        end
    end

    M_ini = all_matrix(:,m+1:m+n);
    s = 340*(Rot*s_ini + R3*ones_2_m)';
    M = 340*(Rot*M_ini + R3*ones_2_n)';

end