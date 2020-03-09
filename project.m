sigma_ur = 2.8 * 10 ^(-10);
sigma_uz = 1.1 * 10 ^(-9);
sigma_vr = 2.32 * 10 ^(-9);
sigma_vz = 6.97 * 10 ^(-9);
Tref = 293.15;
Vmu_ref = 2.39 * 10 ^ (-4);
Eumu = 80200;
Vmfv_ref = 1.61 * 10 ^ (-4);
Evmfu = 56700;
Kmu = 0.4103;
Kmv = 27.2438;
Kmfu = 0.1149;
rq = 0.97;
ro_u = 7 * 10 ^(-7);
ro_v = 7.5 * 10 ^(-7);
Patm = 101300;
%Tcel
Tcel = 25;
%nu,nv
nu = 20.8;
nv = 0.04;

T = Tcel + 273.15;
Cu_amb = (Patm*nu)/R*T;
Cv_amb = (Patm*nv)/R*T;







Ru = (Vmu*Cu)/((Kmu+Cu)*(1+Cv/Kmv));
RV = rq* Ru + (Vmfu/(1+Cu/Kmfu));


elements= [ -(sigma_ur+sigma_uz) -sigma_vr  - sigma_vz;
           -sigma_vr             sigma_ur       0;
           -sigma_vz              0            sigma_vz];
       
       
           