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
TREF = 273.15;
V_MU = 2.39e-4*exp((80200/Rg)*(1/TREF - 1/(TREF + TEMP)));
K_MV = 27.2438;
K_MU = 0.4103;
K_MFU = 0.1149;
MAX_FERM_CO2 = 1.61e-4*exp((56700/Rg)*(1/TREF - 1/(TREF + TEMP)));
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
C(1,1) = patm*nu/(Rg*(TREF + TEMP));
C(2,1) = patm*nv/(Rg*(TREF + TEMP));
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

%% Another test
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
%% Start of program

[matrix1, matrix2] = Generate_stiffnes(elements, coords, S);
[matrixb1, matrixb2] = Boundary_stiffnes(boundaries, coords, R);
[f1, f2] = Boundary_vector(boundaries, coords, C, R);
Initial_C1 = -(matrix1 + matrixb1)\f1;
Initial_C2 = -(matrix2 + matrixb2)\f2;
C = [Initial_C1; Initial_C2];

% Calculating contribution of nonlinear term around ambient concentrations
Jacobian = jacobian_integrand(elements, coords, 0.5, 0, C, VAR) ...
         + jacobian_integrand(elements, coords, 0, 0.5, C, VAR) ...
         + jacobian_integrand(elements, coords, 0.5, 0.5, C, VAR);
Jacobian = Jacobian/6;    
F = Linearized_f_integrand(elements, coords, 0.5, 0, C, VAR) ...
         + Linearized_f_integrand(elements, coords, 0, 0.5, C, VAR) ...
         + Linearized_f_integrand(elements, coords, 0.5, 0.5, C, VAR);
F = F/6;
CONSTANT = ones(2*529,1);
CONSTANT = 3*CONSTANT;
% Creating new linear system that contains linearization
mat = [(matrix1 + matrixb1), zeros(size(coords,1));zeros(size(coords,1)),(matrix2 + matrixb2) ];
mat_lin = mat + Jacobian;
f_vector = - [f1;f2];
F_lin =  f_vector - CONSTANT;
C_lin = mat\F_lin;
dlmwrite('../output/MC_lin',4,'delimiter',' ','precision',12,'-append');
dlmwrite('../output/MC_lin',C_lin','delimiter',' ','precision',12,'-append');

% Executing the nonlinear solver with the calculated starting coefficient
options = optimoptions(@fsolve,'Display','iter',...
    'Algorithm','trust-region',...
    'SpecifyObjectiveGradient',true,'PrecondBandWidth',0, ...
    'FunctionTolerance', 10e-21, 'OptimalityTolerance', 10e-35, ...
    'PlotFcn', 'optimplotstepsize');

model_func = @(coefficients) model(elements, coords, coefficients, VAR, ...
	mat, -f_vector)
[x,fval,exitflag,output] = fsolve(model_func,C_lin,options);


