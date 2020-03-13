import numpy as np

sig_u_r = 2.8e-10
sig_u_z = 1.1e-9
sig_v_r = 2.32e-9
sig_v_z = 6.97e-9
sig_u = np.array([[sig_u_r, 0],[0, sig_u_z]])
sig_v = np.array([[sig_v_r, 0],[0, sig_v_z]])

r1 = 55.4328
r2 = 36.2876
r3 = 29.0488
z1 = 47.039
z2 = 61.2636
z3 = 41.3302

jac = abs((r2-r1)*(z3-z1)-(r3-r1)*(z2-z1))

grad_1 = np.array([-1, -1])
grad_2 = np.array([1, 0])
grad_3 = np.array([0, 1])
grads = [grad_1, grad_2, grad_3]

k = jac*(r1 + r2 + r3)/6

