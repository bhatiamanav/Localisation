function [dist] = distance_3D(X1,Y1,Z1,X2,Y2,Z2)
    sq_dist = (((X1 - X2).^2) + ((Y1-Y2).^2) + ((Z1-Z2).^2));
    dist = sqrt(sq_dist);
end

