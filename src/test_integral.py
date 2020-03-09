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

phi_1 = 1-e-n
phi_2 = e
phi_3 = n
r = r1 + (r2-r1)*e + (r3-r1)*n

g_p1 = grad(phi_1, [e, n])
g_p2 = grad(phi_2, [e, n])
g_p3 = grad(phi_3, [e, n])

sig_m = sp.Matrix([[s_u_r, 0], [0, s_u_z]])
integrand = k1 + k2*e + k3*n
integral = sp.integrate(sp.integrate(integrand, (n, 0, 1-e)), (e, 0, 1))

elem_matrix1 = sp.Matrix(integral*sig_m*g_p1).T*g_p1
elem_matrix1 = elem_matrix1.col_insert(1, sp.Matrix(integral*sig_m*g_p1).T*g_p2)
elem_matrix1 = elem_matrix1.col_insert(2, sp.Matrix(integral*sig_m*g_p1).T*g_p3)

elem_matrix2 = sp.Matrix(integral*sig_m*g_p2).T*g_p1
elem_matrix2 = elem_matrix2.col_insert(1, sp.Matrix(integral*sig_m*g_p2).T*g_p2)
elem_matrix2 = elem_matrix2.col_insert(2, sp.Matrix(integral*sig_m*g_p2).T*g_p3)

elem_matrix3 = sp.Matrix(integral*sig_m*g_p3).T*g_p1
elem_matrix3 = elem_matrix3.col_insert(1, sp.Matrix(integral*sig_m*g_p3).T*g_p2)
elem_matrix3 = elem_matrix3.col_insert(2, sp.Matrix(integral*sig_m*g_p3).T*g_p3)

elem_matrix = elem_matrix1.row_insert(1, elem_matrix2)
elem_matrix = elem_matrix.row_insert(2, elem_matrix3)

Cu = c_u1*phi_1+c_u2*phi_2+c_u3*phi_3
Cv = c_v1*phi_1+c_v2*phi_2+c_v3*phi_3

Ru = Vmu/(1+Cv/Kmv) 
Rv = rq*Ru + Vmfv/(1+Cu/Kmfu)

Ru_comp = Ru*Cu/(Kmu+Cu)
Rv_comp = rq*Ru_comp + Vmfv/(1+Cu/Kmfu)

int_sim1_u = r*Ru*phi_1
int_sim2_u = r*Ru*phi_2
int_sim3_u = r*Ru*phi_3

g_i1_u = grad(int_sim1_u, [c_u1, c_u2, c_u3, c_v1, c_v2, c_v3])
g_i2_u = grad(int_sim2_u, [c_u1, c_u2, c_u3, c_v1, c_v2, c_v3])
g_i3_u = grad(int_sim3_u, [c_u1, c_u2, c_u3, c_v1, c_v2, c_v3])

g_i1_u = g_i1_u.subs([(c_v1, 0),(c_v2, 0),(c_v3, 0)])
g_i2_u = g_i2_u.subs([(c_v1, 0),(c_v2, 0),(c_v3, 0)])
g_i3_u = g_i3_u.subs([(c_v1, 0),(c_v2, 0),(c_v3, 0)])

int_1_u = sp.simplify(
        sp.integrate(sp.integrate(g_i1_u, (n, 0, 1-e)), (e, 0, 1)))
int_2_u = sp.simplify(
        sp.integrate(sp.integrate(g_i2_u, (n, 0, 1-e)), (e, 0, 1)))
int_3_u = sp.simplify(
        sp.integrate(sp.integrate(g_i3_u, (n, 0, 1-e)), (e, 0, 1)))

equil_1_u = int_sim1_u.subs([(c_u1, C_a_u),(c_u2, C_a_u),(c_u3, C_a_u),
    (c_v1, 0),(c_v2, 0),(c_v3, 0)])
equil_2_u = int_sim2_u.subs([(c_u1, C_a_u),(c_u2, C_a_u),(c_u3, C_a_u),
    (c_v1, 0),(c_v2, 0),(c_v3, 0)])
equil_3_u = int_sim3_u.subs([(c_u1, C_a_u),(c_u2, C_a_u),(c_u3, C_a_u),
    (c_v1, 0),(c_v2, 0),(c_v3, 0)])

equil_1_u = sp.simplify(sp.integrate(sp.integrate(equil_1_u, 
    (n, 0, 1-e)), (e, 0, 1)))
equil_2_u = sp.simplify(sp.integrate(sp.integrate(equil_2_u, 
    (n, 0, 1-e)), (e, 0, 1)))
equil_3_u = sp.simplify(sp.integrate(sp.integrate(equil_3_u, 
    (n, 0, 1-e)), (e, 0, 1)))

int_sim1_v = r*Rv*phi_1
int_sim2_v = r*Rv*phi_2
int_sim3_v = r*Rv*phi_3

g_i1_v = grad(int_sim1_v, [c_u1, c_u2, c_u3, c_v1, c_v2, c_v3])
g_i2_v = grad(int_sim2_v, [c_u1, c_u2, c_u3, c_v1, c_v2, c_v3])
g_i3_v = grad(int_sim3_v, [c_u1, c_u2, c_u3, c_v1, c_v2, c_v3])

g_i1_v = g_i1_v.subs([(c_u1, C_a_u),(c_u2, C_a_u),(c_u3, C_a_u),
    (c_v1, 0),(c_v2, 0),(c_v3, 0)])
g_i2_v = g_i2_v.subs([(c_u1, C_a_u),(c_u2, C_a_u),(c_u3, C_a_u),
    (c_v1, 0),(c_v2, 0),(c_v3, 0)])
g_i3_v = g_i3_v.subs([(c_u1, C_a_u),(c_u2, C_a_u),(c_u3, C_a_u),
    (c_v1, 0),(c_v2, 0),(c_v3, 0)])

int_1_v = sp.simplify(sp.integrate(
    sp.integrate(g_i1_v, (n, 0, 1-e)), (e, 0, 1)))
int_2_v = sp.simplify(sp.integrate(
    sp.integrate(g_i2_v, (n, 0, 1-e)), (e, 0, 1)))
int_3_v = sp.simplify(sp.integrate(
    sp.integrate(g_i3_v, (n, 0, 1-e)), (e, 0, 1)))

equil_1_v = int_sim1_v.subs([(c_u1, C_a_u),(c_u2, C_a_u),(c_u3, C_a_u),
    (c_v1, 0),(c_v2, 0),(c_v3, 0)])
equil_2_v = int_sim2_v.subs([(c_u1, C_a_u),(c_u2, C_a_u),(c_u3, C_a_u),
    (c_v1, 0),(c_v2, 0),(c_v3, 0)])
equil_3_v = int_sim3_v.subs([(c_u1, C_a_u),(c_u2, C_a_u),(c_u3, C_a_u),
    (c_v1, 0),(c_v2, 0),(c_v3, 0)])

equil_1_v = sp.simplify(sp.integrate(sp.integrate(equil_1_v, 
    (n, 0, 1-e)), (e, 0, 1)))
equil_2_v = sp.simplify(sp.integrate(sp.integrate(equil_2_v, 
    (n, 0, 1-e)), (e, 0, 1)))
equil_3_v = sp.simplify(sp.integrate(sp.integrate(equil_3_v, 
    (n, 0, 1-e)), (e, 0, 1)))

phi_1 = 1-e
phi_2 = e
r = r1 + (r2-r1)*e
Cu = c_u1*phi_1+c_u2*phi_2
Cv = c_v1*phi_1+c_v2*phi_2

s1 = sp.integrate(r*(Cu-C_a_u)*phi_1, (e, 0, 1))
s2 = sp.integrate(r*(Cu-C_a_u)*phi_2, (e, 0, 1))

#f1 = sp.integrate(r*C_a_u*phi_1, (e, 0, 1))
#f2 = sp.integrate(r*C_a_u*phi_2, (e, 0, 1))
