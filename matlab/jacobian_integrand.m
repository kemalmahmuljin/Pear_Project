function jacob = jacobian_integrand(elements, coords, epsilon, eta, coefficients, VAR)
    % Read necessary variables
	V_MU = VAR(1);
	K_MV = VAR(2);
	K_MU = VAR(3);
	K_MFU = VAR(4);
	MAX_FERM_CO2 = VAR(5);
	RESP_Q = VAR(6);
    % Determine number of nodes
	nodes = size(coefficients,1)/2;
    % Initialize Jacobian
    jacob = zeros(2*nodes, 2*nodes);
    
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

        % Partial derivatives of respirations
		ruu = K_MU*K_MV*V_MU/( (K_MV + c_v)*(K_MU + c_u)^2 );                
		ruv = - c_u*K_MV*V_MU/( (K_MU + c_u)*(K_MV + c_v)^2 );
		rvu =  RESP_Q*ruu - K_MFU*MAX_FERM_CO2/( (K_MFU + c_u)^2 );
		rvv =  RESP_Q*ruv;

        % Appending the values to the Jacobian
		% u_u Quadrant
		jacob(n1, n1) = jacob(n1, n1) + r*ruu*phi_1*phi_1*jac;
		jacob(n1, n2) = jacob(n1, n2) + r*ruu*phi_1*phi_2*jac;
		jacob(n1, n3) = jacob(n1, n3) + r*ruu*phi_1*phi_3*jac;

		jacob(n2, n1) = jacob(n2, n1) + r*ruu*phi_2*phi_1*jac;
		jacob(n2, n2) = jacob(n2, n2) + r*ruu*phi_2*phi_2*jac;
		jacob(n2, n3) = jacob(n2, n3) + r*ruu*phi_2*phi_3*jac;

		jacob(n3, n1) = jacob(n3, n1) + r*ruu*phi_3*phi_1*jac;
		jacob(n3, n2) = jacob(n3, n2) + r*ruu*phi_3*phi_2*jac;
		jacob(n3, n3) = jacob(n3, n3) + r*ruu*phi_3*phi_3*jac;

		% u_v Quadrant
		jacob(n1, n1 + nodes) = jacob(n1, n1 + nodes) + r*ruv*phi_1*phi_1*jac;
		jacob(n1, n2 + nodes) = jacob(n1, n2 + nodes) + r*ruv*phi_1*phi_2*jac;
		jacob(n1, n3 + nodes) = jacob(n1, n3 + nodes) + r*ruv*phi_1*phi_3*jac;
		
		jacob(n2, n1 + nodes) = jacob(n2, n1 + nodes) + r*ruv*phi_2*phi_1*jac;
		jacob(n2, n2 + nodes) = jacob(n2, n2 + nodes) + r*ruv*phi_2*phi_2*jac;
		jacob(n2, n3 + nodes) = jacob(n2, n3 + nodes) + r*ruv*phi_2*phi_3*jac;
	
		jacob(n3, n1 + nodes) = jacob(n3, n1 + nodes) + r*ruv*phi_3*phi_1*jac;
		jacob(n3, n2 + nodes) = jacob(n3, n2 + nodes) + r*ruv*phi_3*phi_2*jac;
		jacob(n3, n3 + nodes) = jacob(n3, n3 + nodes) + r*ruv*phi_3*phi_3*jac;

		% v_u Quadrant
		jacob(n1 + nodes, n1) = jacob(n1 + nodes, n1) - r*rvu*phi_1*phi_1*jac;
		jacob(n1 + nodes, n2) = jacob(n1 + nodes, n2) - r*rvu*phi_1*phi_2*jac;
		jacob(n1 + nodes, n3) = jacob(n1 + nodes, n3) - r*rvu*phi_1*phi_3*jac;

		jacob(n2 + nodes, n1) = jacob(n2 + nodes, n1) - r*rvu*phi_2*phi_1*jac;
		jacob(n2 + nodes, n2) = jacob(n2 + nodes, n2) - r*rvu*phi_2*phi_2*jac;
		jacob(n2 + nodes, n3) = jacob(n2 + nodes, n3) - r*rvu*phi_2*phi_3*jac;

		jacob(n3 + nodes, n1) = jacob(n3 + nodes, n1) - r*rvu*phi_3*phi_1*jac;
		jacob(n3 + nodes, n2) = jacob(n3 + nodes, n2) - r*rvu*phi_3*phi_2*jac;
		jacob(n3 + nodes, n3) = jacob(n3 + nodes, n3) - r*rvu*phi_3*phi_3*jac;

		% v_v Quadrant
		jacob(n1 + nodes, n1 + nodes) = jacob(n1 + nodes, n1 + nodes) - r*rvv*phi_1*phi_1*jac;
		jacob(n1 + nodes, n2 + nodes) = jacob(n1 + nodes, n2 + nodes) - r*rvv*phi_1*phi_2*jac;
		jacob(n1 + nodes, n3 + nodes) = jacob(n1 + nodes, n3 + nodes) - r*rvv*phi_1*phi_3*jac;
		
		jacob(n2 + nodes, n1 + nodes) = jacob(n2 + nodes, n1 + nodes) - r*rvv*phi_2*phi_1*jac;
		jacob(n2 + nodes, n2 + nodes) = jacob(n2 + nodes, n2 + nodes) - r*rvv*phi_2*phi_2*jac;
		jacob(n2 + nodes, n3 + nodes) = jacob(n2 + nodes, n3 + nodes) - r*rvv*phi_2*phi_3*jac;
	
		jacob(n3 + nodes, n1 + nodes) = jacob(n3 + nodes, n1 + nodes) - r*rvv*phi_3*phi_1*jac;
		jacob(n3 + nodes, n2 + nodes) = jacob(n3 + nodes, n2 + nodes) - r*rvv*phi_3*phi_2*jac;
		jacob(n3 + nodes, n3 + nodes) = jacob(n3 + nodes, n3 + nodes) - r*rvv*phi_3*phi_3*jac;	
	end
end