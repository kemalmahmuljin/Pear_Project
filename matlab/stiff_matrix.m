%% update stiffness with linearized integral
Vmu = 2.39 * 10 ^ (-4);
Kmu = 0.4103;
Kmv = 27.2438;
RESP = 0.97;
matrix3 = zeros(m,m);
matrix4 = zeros(m,m);
for i = 1 : m
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
        
        klin = -Vmu/Kmv;
        
        Jlin = (r2 - r1)*(z3 - z1) - (r3 - r1)*(z2 - z1);
        
        klin = klin * Jlin;
        
        s11_lin = ((3 * r1 + r2 + r3) * klin) / 60;
        s12_lin = ((2 * r1 + 2 * r2 + r3) * klin) / 120;
        s13_lin = ((3 * r1 + r2 + 2 * r3) * klin) / 120; 
        
        s21_lin = ((2 * r1 + 2 * r2 + r3) * klin) / 120;
        s22_lin = (( r1 + 3 * r2 + r3) * klin) / 60;
        s23_lin = (( r1 + 2 * r2 + 2 * r3) * klin) / 120;
        
        s31_lin = ((2 * r1 + r2 + 2 * r3) * klin) / 120;
        s32_lin = ((r1 + 2 * r2 + 2 * r3) * klin) / 120;
        s33_lin = ((r1 + r2 + 3 * r3) * klin) / 60;
        
        matrix3(n1,n1) = matrix3(n1,n1)+s11_lin;
        matrix3(n1,n2) = matrix3(n1,n2)+s12_lin;
        matrix3(n1,n3) = matrix3(n1,n3)+s13_lin;
        matrix3(n2,n1) = matrix3(n2,n1)+s12_lin;
        matrix3(n2,n2) = matrix3(n2,n2)+s22_lin;
        matrix3(n3,n1) = matrix3(n3,n1)+s13_lin;
        matrix3(n3,n3) = matrix3(n3,n3)+s33_lin;
        
        
        matrix4(n1,n1) = matrix4(n1,n1)+s11_lin*RESP;
        matrix4(n1,n2) = matrix4(n1,n2)+s12_lin*RESP;
        matrix4(n1,n3) = matrix4(n1,n3)+s13_lin*RESP;
        matrix4(n2,n1) = matrix4(n2,n1)+s12_lin*RESP;
        matrix4(n2,n2) = matrix4(n2,n2)+s22_lin*RESP;
        matrix4(n3,n1) = matrix4(n3,n1)+s13_lin*RESP;
        matrix4(n3,n3) = matrix4(n3,n3)+s33_lin*RESP;
end


Boundary_stiffnes(boundaries);
Generate_stiffnes(); 
jacobian();


