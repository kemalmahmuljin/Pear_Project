function jacob = jacobian_integrandE(elements, coords, V_MU, K_MU, RESP_Q)
    % Determine number of nodes
	nodes = size(coords,1);
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
		jac = abs( (r2 - r1)*(z3 - z1) - (r3 - r1)*(z2 - z1) );
        
        % Partial derivatives of respirations
		ru = V_MU/K_MU;                
        k = jac*ru/120;
        % s11,s12 etc
        u11 = (6*r1 + 2*r2 + 2*r3) * k;       
        u12 = (2*r1 + 2*r2 + r3) * k;        
        u13 = (2*r1 + r2 + 2*r3) * k;   
        %u21 = (2*r1 + 2*r2 + r3) * k;
        u22 = (2*r1 + 6*r2 + 2*r3) *k;
        u23 = (r1 + 2*r2 + 2*r3) * k;
        %u31 = (2*r1 + r2 + 2*r3) * k;
        %u32 = (r1 + 2*r2 + 2*r3) * k;
        u33 = (2*r1 + 2*r2 + 6*r3) * k; 

        % Appending the values to the Jacobian
		% u_u Quadrant
		jacob(n1, n1) = jacob(n1, n1) + u11;
		jacob(n1, n2) = jacob(n1, n2) + u12;
		jacob(n1, n3) = jacob(n1, n3) + u13;

		jacob(n2, n1) = jacob(n2, n1) + u12;
		jacob(n2, n2) = jacob(n2, n2) + u22;
		jacob(n2, n3) = jacob(n2, n3) + u23;

		jacob(n3, n1) = jacob(n3, n1) + u13;
		jacob(n3, n2) = jacob(n3, n2) + u23;
		jacob(n3, n3) = jacob(n3, n3) + u33;

		% v_u Quadrant
		jacob(n1 + nodes, n1) = jacob(n1 + nodes, n1) - RESP_Q*u11;
		jacob(n1 + nodes, n2) = jacob(n1 + nodes, n2) - RESP_Q*u12;
		jacob(n1 + nodes, n3) = jacob(n1 + nodes, n3) - RESP_Q*u13;

		jacob(n2 + nodes, n1) = jacob(n2 + nodes, n1) - RESP_Q*u12;
		jacob(n2 + nodes, n2) = jacob(n2 + nodes, n2) - RESP_Q*u22;
		jacob(n2 + nodes, n3) = jacob(n2 + nodes, n3) - RESP_Q*u23;

		jacob(n3 + nodes, n1) = jacob(n3 + nodes, n1) - RESP_Q*u13;
		jacob(n3 + nodes, n2) = jacob(n3 + nodes, n2) - RESP_Q*u23;
		jacob(n3 + nodes, n3) = jacob(n3 + nodes, n3) - RESP_Q*u33;
	end
end