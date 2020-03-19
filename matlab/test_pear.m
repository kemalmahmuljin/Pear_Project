clear all
clc

global elements boundaries coords fvector fvector_lin stiffness stiffness_lin


boundaries = dlmread('../src/boundaries', ' ', 1, 0);
boundaries = boundaries(:, 1:3);
elements = dlmread('../src/elements', ' ', 1, 0);
elements = elements(:, 1:3);
coords = dlmread('../src/coords', ' ', 1, 0);
coords = coords(:, 1:2);

stiffness = dlmread('../src/stiff', ' ', 1, 0);
stiffness = stiffness(:, 1:end - 1);
stiffness_lin = dlmread('../src/stiff_lin', ' ', 1, 0);
stiffness_lin = stiffness_lin(:, 1:end - 1);
fvector = dlmread('../src/f_vector', ' ', 1, 0);
fvector = fvector(:, 1:end - 1)';
fvector_lin = dlmread('../src/f_vector_lin', ' ', 1, 0);
fvector_lin = fvector_lin(:, 1:end - 1)';

x0 = dlmread('../src/initial_coeff', ' ', 1, 0);
x0= x0(:, 1:end - 1)';

model_function(x0)
options = optimoptions('OptimalityTolerance',1e-9)
[x,fval,exitflag,out] = fsolve(@model_function, x0);
function int_val = integrand(coefficients, epsilon, eta)
global elements boundaries coords fvector fvector_lin stiffness stiffness_lin
	TEMP = 25;
	V_MU = 2.39e-4*exp((80200/8.32)*(1/(273+TEMP) - 1/293));
	K_MV = 27.2438;
	K_MU = 0.4103;
	K_MFU = 0.1149;
	MAX_FERM_CO2 = 1.61e-4*exp((56700/8.32)*(1/(273+TEMP) - 1/293));
	RESP_Q = 0.97; 

	nodes = length(coefficients)/2;
	[elem_num, npe] = size(elements);
	int_val = zeros(2*nodes, 1);
	for i=1:elem_num
		elem = elements(i, :);
		coord1 = coords(elem(1) + 1, :);
		r1 = coord1(1);
		z1 = coord1(2);
		coord2 = coords(elem(2) + 1, :);
		r2 = coord2(1);
		z2 = coord2(2);
		coord3 = coords(elem(3) + 1, :);
		r3 = coord3(1);
		z3 = coord3(2);
		r = r1 + (r2 - r1)*epsilon + (r3 - r1)*eta;
		jac = abs((r2 - r1)*(z3 - z1) + (r3 - r1)*(z2 - z1));

		phi_1 = 1 - epsilon - eta;
		phi_2 = epsilon;
		phi_3 = eta;
		c_u1 = coefficients(elem(1) + 1);
		c_u2 = coefficients(elem(2) + 1);
		c_u3 = coefficients(elem(3) + 1);
		c_u = c_u1*phi_1 + c_u2*phi_2 + c_u3*phi_3;

		c_v1 = coefficients(elem(1) + nodes + 1);
		c_v2 = coefficients(elem(2) + nodes + 1);
		c_v3 = coefficients(elem(3) + nodes + 1);
		c_v = c_v1*phi_1 + c_v2*phi_2 + c_v3*phi_3;

		ru = V_MU*c_u/((K_MU + c_u)*(1 + c_v/K_MV));
		rv = RESP_Q*ru + MAX_FERM_CO2/(1 + c_u/K_MFU);

		int_val(elem(1) + 1) = r*ru*phi_1;
		int_val(elem(2) + 1) = r*ru*phi_2;
		int_val(elem(3) + 1) = r*ru*phi_3;
		int_val(elem(1) + nodes + 1) = -r*rv*phi_1;
		int_val(elem(2) + nodes + 1) = -r*rv*phi_2;
		int_val(elem(3) + nodes + 1) = -r*rv*phi_3;
	end
end

function val = model_function(x)
global elements boundaries coords fvector fvector_lin stiffness stiffness_lin
	val = (integrand(x, 0, 0.5) + ...
	integrand(x, 0.5, 0) + ...
	integrand(x, 0.5, 0.5))/6.0 + ...
	stiffness*x + fvector;
end


