function [si,mj] = obtain_source_mic_locations(Si,Mj)

    row3_1 = sqrt(Si(1,:)- 0.25*(Si(2,:)).*(Si(2,:))-0.25*(Si(3,:)).*(Si(3,:)));
    row3_2 = zeros(1,size(Mj,2));


    si = [-0.5*Si(2,:); -0.5*Si(3,:); row3_1 ];
    mj = [-0.5*Mj(2,:); -0.5*Mj(3,:); row3_2];

    si = real(si);
    mj = real(mj);

end