function [Si, Mj, ti, eps] = return_loc(tij)
	v = 340;

	V = tij;
	[ti,Y] = ctod(tij);

	[Si_cap,Mj_cap] = refine_absol(tij,ti,v);
	[Si,Mj] = compute_locations(Si_cap,Mj_cap);
	eps = 0;
	SX_hat = zeros(1,len(Si,1));
	SY_hat = zeros(1,len(Si,1));
	SZ_hat = zeros(1,len(Si,1));

	MX_hat = zeros(1,len(Mj,1));
	MY_hat = zeros(1,len(Mj,1));
	MZ_hat = zeros(1,len(Mj,1));
	for i = 1:len(Si,1)
	    SX_hat(i) = (-0.5) * Si(i,2);
	    SY_hat(i) = (-0.5) * Si(i,3);
	    SZ_hat(i) = (-0.5) * Si(i,4);
	end

	for i = 1:len(Mj,1)
	    MX_hat(i) = Mj(i,2);
	    MY_hat(i) = Mj(i,3);
	    MZ_hat(i) = Mj(i,4);
	end
	for i = 1:len(Si,1)
	    for j = 1:len(Mj,1)
		eps = eps + (V(i,j) - (ti(i) .* v + distance_3D(SX_hat(i),SY_hat(i),SZ_hat(i),MX_hat(j),MY_hat(j),MZ_hat(j)))).^2;
	    end
	end
end
