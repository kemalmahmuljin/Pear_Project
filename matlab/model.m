function [f, J] = model(elements, coords, coefficients, VAR, stiffness, f_vector)

% Fuction evaluation
f = nonlinear_integrand(elements, coords, 1/3, 1/3, coefficients, VAR);
f = f + stiffness*coefficients + f_vector; 

% Jacobian
J = jacobian_integrand(elements, coords, 1/3, 1/3, coefficients, VAR);
J = J + stiffness;

end

