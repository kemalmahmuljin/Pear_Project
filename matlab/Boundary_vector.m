function [f1, f2] = Boundary_vector(boundaries, coords, C, R)
     
    R_u = R(1,1);
    R_v = R(2,1);    
    Cuamb = C(1,1);
    Cvamb = C(2,1); 
    
    n = size(coords,1);
    f1 = zeros(n,1);
    f2 = zeros(n,1);

    for i = 1:1:size(boundaries,1)
        if boundaries(i,3) == 0 
            n1 = boundaries(i,1)+1;
            n2 = boundaries(i,2)+1;
            %r1,r2
            r1 =  coords(n1,1);
            r2 =  coords(n2,1);
            z1 =  coords(n1,2);
            z2 =  coords(n2,2);
            %length
            len = sqrt((r2-r1)^2+(z2-z1)^2);             
            k1 = (r1/3 + r2/6)*len;
            k2 = (r2/3 + r1/6)*len;
            % Updating f_vector
            f1(n1) = f1(n1) - R_u*Cuamb*k1;
            f1(n2) = f1(n2) - R_u*Cuamb*k2;
            f2(n1) = f2(n1) - R_v*Cvamb*k1;
            f2(n2) = f2(n2) - R_v*Cvamb*k2;             
        end
    end
end