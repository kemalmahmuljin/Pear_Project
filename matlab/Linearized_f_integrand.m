function [F] = Linearized_f_integrand(elements, coords, epsilon, eta, coefficients, VAR )
    % Read necessary variables
	V_MU = VAR(1);
	K_MV = VAR(2);
	K_MU = VAR(3);
	K_MFU = VAR(4);
	MAX_FERM_CO2 = VAR(5);
	RESP_Q = VAR(6);
    
    % Determine number of nodes
    nodes = size(coords,1);
    % Initialize F
    F = zeros(2*nodes,1);
    
    for i=1:1:size(elements,1)
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
        
        % Necessary variables
		r = r1 + (r2 - r1)*epsilon + (r3 - r1)*eta;
		jac = (r2 - r1)*(z3 - z1) - (r3 - r1)*(z2 - z1);

        % Basis functions in triangle
		phi_1 = 1 - epsilon - eta;
		phi_2 = epsilon;
		phi_3 = eta;
        
        % Getting previous coefficients values
		c_u1 = coefficients(n1);
		c_u2 = coefficients(n2);
		c_u3 = coefficients(n3);
        
		c_v1 = coefficients(n1 + nodes);
		c_v2 = coefficients(n2 + nodes);
		c_v3 = coefficients(n3 + nodes);
        
        c_u = c_u1*phi_1 + c_u2*phi_2 + c_u3*phi_3;
		c_v = c_v1*phi_1 + c_v2*phi_2 + c_v3*phi_3;

        % Evaluating the function in the coefficients
        Fcu = K_MV*V_MU*c_u/ ((K_MV + c_v)*(K_MU + c_u) ); 
        Fcv = RESP_Q*Fcu + K_MFU*MAX_FERM_CO2/(K_MV + c_v);
        
        % Updating F vector
        F(n1) = F(n1) + r*Fcu*phi_1*jac;
        F(n2) = F(n2) + r*Fcu*phi_2*jac;
        F(n3) = F(n3) + r*Fcu*phi_3*jac;
        F(n1 + nodes) = F(n1 + nodes) - r*Fcv*phi_1*jac;
        F(n2 + nodes) = F(n2 + nodes) - r*Fcv*phi_2*jac;
        F(n3 + nodes) = F(n3 + nodes) - r*Fcv*phi_3*jac;

        
        % Partial derivatives of respirations
		ruu = K_MU*K_MV*V_MU/( (K_MV + c_v)*(K_MU + c_u)^2 );                
		ruv = - c_u*K_MV*V_MU/( (K_MU + c_u)*(K_MV + c_v)^2 );
		rvu =  RESP_Q*ruu - K_MFU*MAX_FERM_CO2/( (K_MFU + c_u)^2 );
		rvv =  RESP_Q*ruv;

        % updating F vector with partial derivatives
        F(n1) = F(n1) + r*(ruu*c_u + ruv*c_v)*phi_1*jac;
        F(n2) = F(n2) + r*(ruu*c_u + ruv*c_v)*phi_2*jac;
        F(n3) = F(n3) + r*(ruu*c_u + ruv*c_v)*phi_3*jac;
        F(n1 + nodes) = F(n1 + nodes) - r*(rvu*c_u + rvv*c_v)*phi_1*jac;
        F(n2 + nodes) = F(n2 + nodes) - r*(rvu*c_u + rvv*c_v)*phi_2*jac;
        F(n3 + nodes) = F(n3 + nodes) - r*(rvu*c_u + rvv*c_v)*phi_3*jac;       
    end
end