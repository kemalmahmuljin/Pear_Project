import sympy as sp

def grad(funct, var_list):
    res = []
    for var in var_list:
        res.append(sp.diff(funct, var))
    return sp.Matrix(res)

e = sp.Symbol('e')
n = sp.Symbol('n')
r1 = sp.Symbol('r1')
r2 = sp.Symbol('r2')
r3 = sp.Symbol('r3')
c_u1 = sp.Symbol('c_u1')
c_u2 = sp.Symbol('c_u2')
c_u3 = sp.Symbol('c_u3')
c_v1 = sp.Symbol('c_v1')
c_v2 = sp.Symbol('c_v2')
c_v3 = sp.Symbol('c_v3')
C_a_u = sp.Symbol('C_a_u')
C_a_v = sp.Symbol('C_a_v')

k1 = r1
k2 = r2 - r1
k3 = r3 - r1
s_u_r = sp.Symbol('s_u_r')
s_u_z = sp.Symbol('s_u_z')
s_v_r = sp.Symbol('s_v_r')
s_v_r = sp.Symbol('s_v_z')
J = sp.Symbol('|J|')
Vmu = sp.Symbol('Vmu')
Kmu = sp.Symbol('Kmu')
Kmv = sp.Symbol('Kmv')
Kmfu = sp.Symbol('Kmfu')
Vmfv = sp.Symbol('Vmfv')
rq = sp.Symbol('rq')

coeff = [c_u1, c_u2, c_u3, c_v1, c_v2, c_v3]
phi_1 = sp.Symbol('phi_1')
phi_2 = sp.Symbol('phi_2')
phi_3 = sp.Symbol('phi_3')
r = sp.Symbol('r')

Cu = c_u1*phi_1+c_u2*phi_2+c_u3*phi_3
Cv = c_v1*phi_1+c_v2*phi_2+c_v3*phi_3

Ru = Vmu/(1+Cv/Kmv) 
Rv = rq*Ru + Vmfv/(1+Cu/Kmfu)

Ru_comp = Ru*Cu/(Kmu+Cu)
Rv_comp = rq*Ru_comp + Vmfv/(1+Cu/Kmfu)

function_vector = []
function_vector.append(r*Ru_comp*phi_1)
function_vector.append(r*Ru_comp*phi_2)
function_vector.append(r*Ru_comp*phi_3)
function_vector.append(r*Rv_comp*phi_1)
function_vector.append(r*Rv_comp*phi_2)
function_vector.append(r*Rv_comp*phi_3)

jacobian = []
for func in function_vector:
    jacobian.append(sp.simplify(grad(func, coeff)))

