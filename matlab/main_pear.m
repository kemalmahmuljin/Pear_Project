%% Importing needed matrices
boundaries = dlmread('../output/boundaries', ' ', 1, 0);
boundaries = boundaries(:, 1:3);
elements = dlmread('../output/elements', ' ', 1, 0);
elements = elements(:, 1:3);
coords = dlmread('../output/coords', ' ', 1, 0);
coords = coords(:, 1:2);



%% Declaring some constants
% Tuneable parameters
TEMP = 25;
nu = 20.8/100.0;
nv = 0.04/100.0;
% General constants
Rg = 8.314;
TREF = 293.15; % found another bug
T0 = 273.15;
V_MU = 2.39e-4*exp((80200/Rg)*(1/TREF - 1/(T0 + TEMP)));
K_MV = 27.2438;
K_MU = 0.4103;
K_MFU = 0.1149;
MAX_FERM_CO2 = 1.61e-4*exp((56700/Rg)*(1/TREF - 1/(T0 + TEMP)));
RESP_Q = 0.97;
patm = 101300; 
% sigma's
S  = zeros(2,2);
S(1,1) = 2.8 * 10 ^ (-10); % sur
S(1,2) = 1.1 * 10 ^ (-9); % suz 
S(2,1) = 2.32 * 10 ^ (-9); % svr
S(2,2) = 6.97 * 10 ^ (-9); % svz
% respirations
R = zeros(2,1);
R(1,1) = 7 * 10 ^ (-7);
R(2,1) = 7.5 * 10 ^ (-7);
% Ambient concentrations
C = zeros(2,1);
C(1,1) = patm*nu/(Rg*(T0 + TEMP));
C(2,1) = patm*nv/(Rg*(T0 + TEMP));
% Setting VAR
VAR = zeros(6,1);
VAR(1) = V_MU;
VAR(2) = K_MV;
VAR(3) = K_MU;
VAR(4) = K_MFU;
VAR(5) = MAX_FERM_CO2;
VAR(6) = RESP_Q;

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
%% Start of program

[matrix1, matrix2] = Generate_stiffnes(elements, coords, S);
[matrixb1, matrixb2] = Boundary_stiffnes(boundaries, coords, R);
[f1, f2] = Boundary_vector(boundaries, coords, C, R);
Initial_C1 = -(matrix1 + matrixb1)\f1;
Initial_C2 = -(matrix2 + matrixb2)\f2;
C_start = [Initial_C1; Initial_C2];
% dlmwrite('../output/MC_lin',4,'delimiter',' ','precision',12,'-append');
% dlmwrite('../output/MC_lin',C','delimiter',' ','precision',12,'-append');
% Calculating contribution of nonlinear term around ambient concentrations
Jacobian = jacobian_integrand(elements, coords, 0.5, 0  , C_start, VAR) ...
         + jacobian_integrand(elements, coords, 0  , 0.5, C_start, VAR) ...
         + jacobian_integrand(elements, coords, 0.5, 0.5, C_start, VAR);
Jacobian = Jacobian/6;    
F =        Linearized_f_integrand(elements, coords, 0.5, 0  , C_start, VAR) ...
         + Linearized_f_integrand(elements, coords, 0  , 0.5, C_start, VAR) ...
         + Linearized_f_integrand(elements, coords, 0.5, 0.5, C_start, VAR);
F = F/6;
mat1 = [(matrix1), zeros(size(coords,1));zeros(size(coords,1)),(matrix2 ) ];

% Creating new linear system that contains linearization
mat = [(matrix1 + matrixb1), zeros(size(coords,1));zeros(size(coords,1)),(matrix2 + matrixb2) ];
mat_lin = mat + Jacobian;
f_vector = - [f1;f2];
F_lin =  f_vector - F;
C_lin = mat_lin\F_lin;
%print this
% C_lin =F;
dlmwrite('../output/MC_lin',4,'delimiter',' ','precision',12,'-append');
dlmwrite('../output/MC_lin',C_lin,'delimiter',' ','precision',12,'-append');

% Executing the nonlinear solver with the calculated starting coefficient
options = optimoptions(@fsolve,'Display','iter',...
    'Algorithm','trust-region',...
    'SpecifyObjectiveGradient',true,'PrecondBandWidth',0 , ...
    'FunctionTolerance', 1e-20, 'OptimalityTolerance', 5e-19);%, ...
    %'PlotFcn', 'optimplotstepsize');

model_func = @(coefficients) model(elements, coords, coefficients, VAR, ...
	mat, -f_vector)
[x,fval,exitflag,output] = fsolve(model_func,C_lin,options);
dlmwrite('../output/MC_lin',4,'delimiter',' ','precision',12,'-append');
dlmwrite('../output/MC_lin',x','delimiter',' ','precision',12,'-append');

