S_ur = 2.8 * 10 ^ (-10);
S_uz = 1.1 * 10 ^ (-9);
S_vr = 2.32 * 10 ^ (-9);
S_vz = 6.97 * 10 ^ (-9);

%% update stifness matrix
m = size(elements,1);
matrix1 = zeros(m,m);
matrix2 = zeros(m,m);
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
        % k
        k = (r1+r2+r3)/6;   
        % Jacobian
        J = abs((r2 - r1)*(z3 - z1) - (r3 - r1)*(z2 - z1));
        % s11,s12 etc
        s11u = (S_ur + S_uz) * k ;       
        s12u = -S_ur * k;        
        s13u = -S_uz * k;
        s22u =  S_ur * k;
        s33u =  S_uz * k; 
        
        matrix1(n1,n1) = matrix1(n1,n1)+s11u;
        matrix1(n1,n2) = matrix1(n1,n2)+s12u;
        matrix1(n1,n3) = matrix1(n1,n3)+s13u;
        matrix1(n2,n1) = matrix1(n2,n1)+s12u;
        matrix1(n2,n2) = matrix1(n2,n2)+s22u;
        matrix1(n3,n1) = matrix1(n3,n1)+s13u;
        matrix1(n3,n3) = matrix1(n3,n3)+s33u;
        
        %start for matrix 2
        %s11,s12 etc
        s11v = (S_vr + S_vz) * k ;       
        s12v = -S_vr * k;        
        s13v = -S_vz * k;
        s22v =  S_vr * k;
        s33v =  S_vz * k;
        
        matrix2(n1,n1) = matrix2(n1,n1)+s11v;
        matrix2(n1,n2) = matrix2(n1,n2)+s12v;
        matrix2(n1,n3) = matrix2(n1,n3)+s13v;
        matrix2(n2,n1) = matrix2(n2,n1)+s12v;
        matrix2(n2,n2) = matrix2(n2,n2)+s22v;
        matrix2(n3,n1) = matrix2(n3,n1)+s13v;
        matrix2(n3,n3) = matrix2(n3,n3)+s33v;
end
% matrix_store = matrix1;
%spy(matrix1)
%spy(matrix2)
        

%% Boundary
 R_u = 7 * 10 ^ (-7);
 R_v = 7.5 * 10 ^ (-7);
 n = size(boundaries,1);
 matrixb1 = zeros(n,n);
 matrixb2 = zeros(n,n);

 for i = 1 : n
     if boundaries(i,3) == 0 
         n1 = boundaries(i,1)+1;
         n2 = boundaries(i,2)+1;
      
      %r1,r2
      r1 =  coords(n1,1);
      r2 =  coords(n2,1);
      %set the initial coordinates to (0,0)
      r11=0;
      r22=0;
      %length
      len = sqrt((r2-r22)^2+(r1-r11)^2);
      %set the new starting coordinate to calculate the length next time
      r11=r1;
      r22=r2;
      %s11 etc
      s11bu = R_u* len * ( r1/4 + r2/12);
      s12bu = R_u* len * ( r1 + r2)/12;
      s22bu = R_u* len * ( r2/4 + r1/12);
      
      matrixb1(n1,n1) = matrixb1(n1,n1)+s11bu;
      matrixb1(n1,n1) = matrixb1(n1,n2)+s12bu;
      matrixb1(n2,n1) = matrixb1(n2,n1)+s12bu;
      matrixb1(n2,n2) = matrixb1(n1,n1)+s22bu;
      
      s11bv = R_v* len * ( r1/4 + r2/12);
      s12bv = R_v* len * ( r1 + r2)/12;
      s22bv = R_v* len * ( r2/4 + r1/12);
      
      matrixb2(n1,n1) = matrixb2(n1,n1)+s11bv;
      matrixb2(n1,n1) = matrixb2(n1,n2)+s12bv;
      matrixb2(n2,n1) = matrixb2(n2,n1)+s12bv;
      matrixb2(n2,n2) = matrixb2(n1,n1)+s22bv;
     end 
 end
% figure;
% spy(matrix_store);
% figure;
% spy(matrix1 - matrix_store)

%% update stiffness with linear integral

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
        
        Jlin = abs((r2 - r1)*(z3 - z1) - (r3 - r1)*(z2 - z1));
        
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