clear all; close all
%% Importing needed matrices
%% These are configs for GMSH
%{
boundaries = dlmread('../output/boundaries', ' ', 1, 0);
boundaries = boundaries(:, 1:3);
elements = dlmread('../output/elements', ' ', 1, 0);
elements = elements(:, 1:3);
coords = dlmread('../output/coords', ' ', 1, 0);
coords = coords(:, 1:2);
%}
%% These are configs for Matlab Mesh
%{
boundaries = dlmread('../output/M_boundaries');
boundaries = boundaries(:, 1:3);
elements = dlmread('../output/M_elements');
elements = elements(:, 1:3);
coords = dlmread('../output/M_coords');
coords = coords(:, 1:2);
%}
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
% Setting VAR (variables)
VAR    = zeros(6,1);
VAR(1) = V_MU;
VAR(2) = K_MV;
VAR(3) = K_MU;
VAR(4) = K_MFU;
VAR(5) = MAX_FERM_CO2;
VAR(6) = RESP_Q;


%% Start of program

[matrix1, matrix2] = Generate_stiffnes(elements, coords, S);
[matrixb1, matrixb2] = Boundary_stiffnes(boundaries, coords, R);
[f1, f2] = Boundary_vector(boundaries, coords, C, R);
f_vector = - [f1;f2];
matrix_amb = [( matrix1 + matrixb1 ), zeros( size(coords,1) ); zeros( size(coords,1) ),( matrix2 + matrixb2 )];
C_start = matrix_amb\f_vector; % should be Camb
% writing results to file
delete('../output/MC_lin');
dlmwrite('../output/MC_lin',4,'delimiter',' ','precision',12,'-append');
dlmwrite('../output/MC_lin',C_start','delimiter',' ','precision',12,'-append');

% Calculating contribution of nonlinear term around ambient concentrations
% Jacobian = jacobian_integrandE(elements, coords, V_MU, K_MU, RESP_Q); 
Jacobian = jacobian_integrand(elements, coords, 0.5, 0  , C_start, VAR) ...
         + jacobian_integrand(elements, coords, 0  , 0.5, C_start, VAR) ...
         + jacobian_integrand(elements, coords, 0.5, 0.5, C_start, VAR);
Jacobian = Jacobian/6;
F =        Linearized_f_integrand(elements, coords, 0.5, 0  , C_start, VAR) ...
         + Linearized_f_integrand(elements, coords, 0  , 0.5, C_start, VAR) ...
         + Linearized_f_integrand(elements, coords, 0.5, 0.5, C_start, VAR);
F = F/6;

% Creating new linear system that contains linearization
mat_lin =  matrix_amb + Jacobian;
% F(1:size(coords,1)) = - 0.5*F(1:size(coords,1));
C_lin = mat_lin\(f_vector - F);
% writing results to file
delete('../output/MC_lin');
dlmwrite('../output/MC_lin',4,'delimiter',' ','precision',12,'-append');
dlmwrite('../output/MC_lin',C_lin','delimiter',' ','precision',12,'-append');
% dlmwrite('../output/MC_lin',(F-f_vector)','delimiter',' ','precision',12,'-append');
% dlmwrite('../output/MC_lin',(Jacobian*C_lin)','delimiter',' ','precision',12,'-append');

% Executing the nonlinear solver with the calculated starting coefficient
options = optimoptions(@fsolve,'Display','iter',...
    'Algorithm','trust-region',...
    'SpecifyObjectiveGradient',true,'PrecondBandWidth',0 , ...
    'FunctionTolerance', 1e-30, 'OptimalityTolerance', 4e-22)%, ...
    %'PlotFcn', 'optimplotresnorm');
model_func = @(coefficients) model(elements, coords, coefficients, VAR, ...
	matrix_amb, -f_vector)

[x,fval,exitflag,output] = fsolve(model_func,C_lin,options);

% writing results to file
delete('../output/MC_lin');
dlmwrite('../output/MC_lin',4,'delimiter',' ','precision',12,'-append');
dlmwrite('../output/MC_lin',x','delimiter',' ','precision',12,'-append');




% All main tests are in C++ 
%% Test cases to check boundary
%{
% TEST coeff
CU = ones(529,1);
CV = ones(529,1);
CU = CU*10;
CV = CV/10;
TESTC = [CU;CV];
% Ambient Coeff
U = ones(529,1)*C(1,1);
V = ones(529,1)*C(2,1);
CC = [U;V];

%stiffness 1
[matrix1, matrix2] = Generate_stiffnes(elements, coords, S);
[matrixb1, matrixb2] = Boundary_stiffnes(boundaries, coords, R);
mat = [(matrix1 + matrixb1), zeros(size(coords,1));zeros(size(coords,1)),(matrix2 + matrixb2) ];
%F 1
FLIN = mat*TESTC;
[f1, f2] = Boundary_vector(boundaries, coords, TESTC, R);
FLIN = [f1; f2];
dlmwrite('../output/f_vector',4,'delimiter',' ','precision',12,'-append');
dlmwrite('../output/f_vector',FLIN','delimiter',' ','precision',12,'-append');

%stiffness 2
Jacobian = jacobian_integrand(elements, coords, 0.5, 0, CC, VAR) ...
         + jacobian_integrand(elements, coords, 0, 0.5, CC, VAR) ...
         + jacobian_integrand(elements, coords, 0.5, 0.5, CC, VAR);
Jacobian = Jacobian/6; 
mat =Jacobian;
TESTF = mat*TESTC;
dlmwrite('../output/f_vector_lin',4,'delimiter',' ','precision',12,'-append');
dlmwrite('../output/f_vector_lin',TESTF','delimiter',' ','precision',12,'-append');
%}

%% Another test by adding tiny constant to linearized f_vector
%{
[matrix1, matrix2] = Generate_stiffnes(elements, coords, S);
[matrixb1, matrixb2] = Boundary_stiffnes(boundaries, coords, R);
mat = [(matrix1 + matrixb1), zeros(size(coords,1));zeros(size(coords,1)),(matrix2 + matrixb2) ];
[f1, f2] = Boundary_vector(boundaries, coords, C, R);
CONSTANT = ones(2*529,1);
CONSTANT = 3*CONSTANT*10^-19;
F = [f1;f2] + CONSTANT;
C_lin = - mat\F;
dlmwrite('../output/MC_lin',4,'delimiter',' ','precision',12,'-append');
dlmwrite('../output/MC_lin',C_lin','delimiter',' ','precision',12,'-append');
%}
%% Generating linearly increasing Coefficients in R and testing magnitude
%{
x0 = dlmread('../output/calculated_coeff');
% x0 = (3*x0);
[matrix1, matrix2] = Generate_stiffnes(elements, coords, S);
mat = [(matrix1), zeros(size(coords,1));zeros(size(coords,1)),(matrix2) ];
F_diffusion = mat*x0';
dlmwrite('../output/MC_lin',4,'delimiter',' ','precision',12,'-append');
dlmwrite('../output/MC_lin',F_diffusion','delimiter',' ','precision',12,'-append');
[f1, f2] = Boundary_vector(boundaries, coords, C, R);
f_vector = - [f1;f2];
dlmwrite('../output/MC_lin',4,'delimiter',' ','precision',12,'-append');
dlmwrite('../output/MC_lin',f_vector','delimiter',' ','precision',12,'-append');
[matrixb1, matrixb2] = Boundary_stiffnes(boundaries, coords, R);
matb = [(matrixb1), zeros(size(coords,1));zeros(size(coords,1)),(matrixb2) ];
F_diffusion = matb*x0';
% dlmwrite('../output/MC_lin',4,'delimiter',' ','precision',12,'-append');
% dlmwrite('../output/MC_lin',F_diffusion','delimiter',' ','precision',12,'-append');
%}