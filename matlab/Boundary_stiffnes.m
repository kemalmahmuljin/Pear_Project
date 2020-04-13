%% Boundary
function [matrixb1, matrixb2] = Boundary_stiffnes(boundaries, coords, R)
     R_u = R(1,1);
     R_v = R(2,1);
     n = size(coords,1);
     matrixb1 = zeros(n,n);
     matrixb2 = zeros(n,n);

     for i = 1:1:size(boundaries,1)
         if boundaries(i,3) == 0 % so if \elem \Gamma_2
            n1 = boundaries(i,1)+1;
            n2 = boundaries(i,2)+1;
            %r1,r2
            r1 =  coords(n1,1);
            r2 =  coords(n2,1);
            z1 =  coords(n1,2);
            z2 =  coords(n2,2);
            %length
            len = sqrt((r2-r1)^2+(z2-z1)^2);
            %s11 etc
            s11bu = R_u* len * ( r1/4 + r2/12);
            s12bu = R_u* len * ( r1 + r2)/12;
            s22bu = R_u* len * ( r2/4 + r1/12);

            matrixb1(n1,n1) = matrixb1(n1,n1)+s11bu;
            matrixb1(n1,n2) = matrixb1(n1,n2)+s12bu;
            matrixb1(n2,n1) = matrixb1(n2,n1)+s12bu;
            matrixb1(n2,n2) = matrixb1(n2,n2)+s22bu;

            s11bv = R_v* len * ( r1/4 + r2/12);
            s12bv = R_v* len * ( r1 + r2)/12;
            s22bv = R_v* len * ( r2/4 + r1/12);

            matrixb2(n1,n1) = matrixb2(n1,n1)+s11bv;
            matrixb2(n1,n2) = matrixb2(n1,n2)+s12bv;
            matrixb2(n2,n1) = matrixb2(n2,n1)+s12bv;
            matrixb2(n2,n2) = matrixb2(n2,n2)+s22bv;
         end 
     end
end