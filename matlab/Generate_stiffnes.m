function [matrix1,matrix2]=Generate_stiffnes(elements, coords, S) 
    S_ur = S(1,1);
    S_uz = S(1,2);
    S_vr = S(2,1);
    S_vz = S(2,2);


    %% update stifness matrix
    m = size(coords,1);
    matrix1 = zeros(m,m);
    matrix2 = zeros(m,m);
    for i = 1:1:size(elements, 1)
            %nodes
            n1 = elements(i,1)+1;
            n2 = elements(i,2)+1;
            n3 = elements(i,3)+1;       
            % r1,r2,r3
            r1 =  coords(n1,1);
            r2 =  coords(n2,1);
            r3 =  coords(n3,1);
            % z1,z2,z3
            z1 =  coords(n1,2);
            z2 =  coords(n2,2);
            z3 =  coords(n3,2);  
            % Jacobian
            J = abs( (r2 - r1)*(z3 - z1) - (r3 - r1)*(z2 - z1) );
            % k
            k = (r1+r2+r3)/(6*J); 
            z12 = z1 - z2;
            z23 = z2 - z3;
            z31 = z3 - z1;
            
            r21 = r2 - r1;
            r32 = r3 - r2;
            r13 = r1 - r3;            
            
            % s11,s12 etc
            s11u = (S_ur*z23^2   + S_uz*r32^2  ) * k;       
            s12u = (S_ur*z23*z31 + S_uz*r32*r13) * k;        
            s13u = (S_ur*z23*z12 + S_uz*r32*r21) * k;   
            s22u = (S_ur*z31^2   + S_uz*r13^2  ) * k;
            s23u = (S_ur*z31*z12 + S_uz*r13*r21) * k;
            s33u = (S_ur*z12^2   + S_uz*r21^2  ) * k; 

            matrix1(n1,n1) = matrix1(n1,n1)+s11u;
            matrix1(n1,n2) = matrix1(n1,n2)+s12u;
            matrix1(n1,n3) = matrix1(n1,n3)+s13u;
            matrix1(n2,n1) = matrix1(n2,n1)+s12u;
            matrix1(n2,n2) = matrix1(n2,n2)+s22u;
            matrix1(n2,n3) = matrix1(n2,n3)+s23u;
            matrix1(n3,n1) = matrix1(n3,n1)+s13u;
            matrix1(n3,n2) = matrix1(n3,n2)+s23u;
            matrix1(n3,n3) = matrix1(n3,n3)+s33u;

            %start for matrix 2
            s11v = (S_vr*z23^2   + S_vz*r32^2  ) * k;       
            s12v = (S_vr*z23*z31 + S_vz*r32*r13) * k;        
            s13v = (S_vr*z23*z12 + S_vz*r32*r21) * k;   
            s22v = (S_vr*z31^2   + S_vz*r13^2  ) * k;
            s23v = (S_vr*z31*z12 + S_vz*r13*r21) * k;
            s33v = (S_vr*z12^2   + S_vz*r21^2  ) * k; 
            
       
            matrix2(n1,n1) = matrix2(n1,n1)+s11v;
            matrix2(n1,n2) = matrix2(n1,n2)+s12v;
            matrix2(n1,n3) = matrix2(n1,n3)+s13v;
            matrix2(n2,n1) = matrix2(n2,n1)+s12v;
            matrix2(n2,n2) = matrix2(n2,n2)+s22v;
            matrix2(n2,n3) = matrix2(n2,n3)+s23v;
            matrix2(n3,n1) = matrix2(n3,n1)+s13v;
            matrix2(n3,n2) = matrix2(n3,n2)+s23v;
            matrix2(n3,n3) = matrix2(n3,n3)+s33v;
    end
end