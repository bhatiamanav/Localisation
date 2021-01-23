ti_ini = [0.125,0.200,0.345,0.515,0.690];
tij = zeros(5,10);
SX = [0,0,0,1,0];
SY = [0,0,1,0,1];
SZ = [0,1,0,0,1];
MX = [1,1,1,1/3,0,0,1/2,0,4/9,1/5];
MY = [0,1,1,0,3/4,0,2/3,3/7,0,4/7];
MZ = [1,0,1,0,0,2/5,0,3/5,6/11,5/11];
v = 340;
for i = 1:5
    for j = 1:10
        tij(i,j) = (ti_ini(i).*v) + distance_3D(SX(i),SY(i),SZ(i),MX(j),MY(j),MZ(j));
    end
end
tij = tij .* (1/v);
[ti,Y] = ctod(tij);
disp(ti);
[Si_cap,Mj_cap] = refine_absol(tij,ti,v);
[Si,Mj] = compute_locations(Si_cap,Mj_cap);
