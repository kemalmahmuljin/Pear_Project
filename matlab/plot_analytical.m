%% Importing everything
coords = dlmread('../output/coordsLarge', ' ', 1, 0);
coords = coords(:, 1:2);

%% Declaring some constants
% Tuneable parameters
% case (1) of Orchard
TEMP = 25; % Temperature in °C
nu   = 20.8/100.0; % Fraction of O2 in atmosphere
nv   = 0.04/100.0; % Fraction of CO2 in atmosphere
% General constants
Rg           = 8.314; % universal gas constant
TREF         = 293.15; % reference temperature
T0           = 273.15; % 0 °C
V_MU         = 2.39e-4 * exp( ( 80200 / Rg ) * ( 1 / TREF - 1 / ( T0 + TEMP ) ) ); % Maximum oxygen consumption rate
K_MV         = 27.2438; % Michaelis Menten non competitive carbon inhibition
K_MU         = 0.4103; % Michaelis Menten oxygen consumption
K_MFU        = 0.1149; % Michaelis Menten oxygen inhibition on fermentative carbon production
MAX_FERM_CO2 = 1.61e-4 * exp( ( 56700 / Rg ) * ( 1 / TREF - 1 / ( T0 + TEMP ) ) ); % Maximum fermentative carbon production rate
RESP_Q       = 0.97; % Respiration Quotient
patm         = 101300;  % Atmospheric pressure
% Radial and Axial diffusivity of 02/CO2 in pear tissue
S      = zeros(2,2);
S(1,1) = 2.8 * 10 ^ (-10); % sur
S(1,2) = 1.1 * 10 ^ (-9); % suz 
S(2,1) = 2.32 * 10 ^ (-9); % svr
S(2,2) = 6.97 * 10 ^ (-9); % svz
% Convective Mass transfer coefficients
R      = zeros(2,1);
R(1,1) = 7 * 10 ^ (-7); % rho_u
R(2,1) = 7.5 * 10 ^ (-7); % rho_v
% Ambient O2/CO2 concentrations
C      = zeros(2,1);
C(1,1) = patm * nu / ( Rg * ( T0 + TEMP ) ); %C_U_AMB
C(2,1) = patm * nv / ( Rg * ( T0 + TEMP ) ); %C_V_AMB


%% calculating constants
Xcoords = coords(1:size(coords,1), 1);
Ycoords = coords(1:size(coords,1), 2);

radius = 0.06;
a = V_MU / ( K_MU * S(1,1) );
C1 = R(1,1) * C(1,1) / ( R(1,1) * besselj(0,sqrt(-a)*radius) + S(1,1) * sqrt(-a) * besselj(1,sqrt(-a)*radius) );
C3 = C1 * RESP_Q * S(1,1) * (- sqrt(-a) * besselj(1,sqrt(-a)*radius) + R(2,1) * besselj(0, sqrt(-a)*radius) / S(2,1) ) / R(2,1) + C(2,1);

MC_lin = zeros(2*size(coords, 1), 1);
for i=1:1:size(coords, 1)
    r = sqrt(Xcoords(i)^2 + ( Ycoords(i) - 0.07 )^2);
 
    cuu = C1 * besselj(0 ,sqrt(-a)*r);
    cvv =- C1 * RESP_Q * S(1,1) * besselj(0,sqrt(-a)*r) / S(2,1) + C3;
    MC_lin(i) = cuu;
    MC_lin(size(coords, 1) + i) = cvv;
end

% writing results to file
delete('../output/MC_lin');
dlmwrite('../output/MC_lin',4,'delimiter',' ','precision',12,'-append');
dlmwrite('../output/MC_lin',MC_lin','delimiter',' ','precision',12,'-append');