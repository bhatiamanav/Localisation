function tij = compute_tij(V,s,M,Ti)
    v = 340;
    a = size(V,1);
    V_s = [0 0]; %Virtual source
    V = [V ; V(1,:)]; %Vertices of room
    %disp(V);
    for i = 1:a
        s1 = V(i,:);
        s2 = V(i+1,:);
        slop1 = (s2(2) - s1(2))/(s2(1)-s1(1));
        slop2 = -1/slop1;
        if slop2==0
            c1 = -s1(1);
            v_s_x = -s(1)-2*c1;
            v_s_y = s(2);
            v_s = [v_s_x, v_s_y];
        elseif slop1==0
            Md = zeros(1,2);
            Md(1) = s(1);
            Md(2) = s1(2);
            v_s = 2 .* Md  - s;
        else
            c1 = s1(2) - slop1*s1(1);
            c2 = s(2) - slop2*s(1);
            x_mid = (c2 - c1)*slop1/((slop1^2)+1);
            y_mid = (c2*(slop1^2)+c1)/((slop1^2)+1);
            mid = [x_mid y_mid];
            v_s = 2*mid - s;
        end
        V_s = [V_s ; v_s];
    end
    
    V_s = V_s(2:a+1,:);
    %disp(V_s/340);
    V = V(1:a,:);
    S = [s ; V_s];
    %disp(S/340);
    sz_S = size(S,1);
    sz_M = size(M,1);
    prop = zeros(sz_S,sz_M);
    ones_vec = ones(1,sz_M);
    for j=1:sz_S
        for k=1:sz_M
            prop(j,k) = (norm(S(j,:)-M(k,:)))/v;
        end
    end
    %disp(prop);
    tij = Ti * ones_vec + prop;
end

