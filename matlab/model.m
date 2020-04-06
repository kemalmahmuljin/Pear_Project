function [f, J] = model(elements, coords, coefficients, VAR, stiffness, f_vector)

% Fuction evaluation

f = nonlinear_integrand(elements, coords, 0, 0.5, coefficients, VAR) + ...
    nonlinear_integrand(elements, coords, 0.5, 0, coefficients, VAR) + ...
    nonlinear_integrand(elements, coords, 0.5, 0.5, coefficients, VAR);

f = f/6.0 + stiffness*coefficients + f_vector; 

% Jacobian

J = jacobian_integrand(elements, coords, 0, 0.5, coefficients, VAR) + ...
    jacobian_integrand(elements, coords, 0.5, 0, coefficients, VAR) + ...
    jacobian_integrand(elements, coords, 0.5, 0.5, coefficients, VAR);
J = J/6.0 + stiffness;
end

