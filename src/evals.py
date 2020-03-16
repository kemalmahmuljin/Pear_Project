import sympy as sp
import math
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







#########################################
#### playground #########################
#########################################

### constant belonging to F_u
'''
T1 = C_a_u*Kmv*Vmu*(4*C_a_u*C_a_v*r1 + 2*C_a_u*C_a_v*r2 + 2*C_a_u*C_a_v*r3 + 2*C_a_u*Kmv*r1 + C_a_u*Kmv*r2 + C_a_u*Kmv*r3 + 2*C_a_v*Kmu*r1 + C_a_v*Kmu*r2 + C_a_v*Kmu*r3)/(24*(C_a_u**2*C_a_v**2 + 2*C_a_u**2*C_a_v*Kmv + C_a_u**2*Kmv**2 + 2*C_a_u*C_a_v**2*Kmu + 4*C_a_u*C_a_v*Kmu*Kmv + 2*C_a_u*Kmu*Kmv**2 + C_a_v**2*Kmu**2 + 2*C_a_v*Kmu**2*Kmv + Kmu**2*Kmv**2))

T2 = C_a_u*Kmv*Vmu*(2*C_a_u*C_a_v*r1 + 4*C_a_u*C_a_v*r2 + 2*C_a_u*C_a_v*r3 + C_a_u*Kmv*r1 + 2*C_a_u*Kmv*r2 + C_a_u*Kmv*r3 + C_a_v*Kmu*r1 + 2*C_a_v*Kmu*r2 + C_a_v*Kmu*r3)/(24*(C_a_u**2*C_a_v**2 + 2*C_a_u**2*C_a_v*Kmv + C_a_u**2*Kmv**2 + 2*C_a_u*C_a_v**2*Kmu + 4*C_a_u*C_a_v*Kmu*Kmv + 2*C_a_u*Kmu*Kmv**2 + C_a_v**2*Kmu**2 + 2*C_a_v*Kmu**2*Kmv + Kmu**2*Kmv**2))

T3 = C_a_u*Kmv*Vmu*(2*C_a_u*C_a_v*r1 + 2*C_a_u*C_a_v*r2 + 4*C_a_u*C_a_v*r3 + C_a_u*Kmv*r1 + C_a_u*Kmv*r2 + 2*C_a_u*Kmv*r3 + C_a_v*Kmu*r1 + C_a_v*Kmu*r2 + 2*C_a_v*Kmu*r3)/(24*(C_a_u**2*C_a_v**2 + 2*C_a_u**2*C_a_v*Kmv + C_a_u**2*Kmv**2 + 2*C_a_u*C_a_v**2*Kmu + 4*C_a_u*C_a_v*Kmu*Kmv + 2*C_a_u*Kmu*Kmv**2 + C_a_v**2*Kmu**2 + 2*C_a_v*Kmu**2*Kmv + Kmu**2*Kmv**2))


c = sp.simplify(T2/(r1 + 2*r2 + r3))
print('result: ', sp.simplify(T3/c))

trial = (C_a_u*C_a_v + C_a_u*Kmv + C_a_v*Kmu + Kmv*Kmu)**2
print('trial: ', sp.simplify(trial*c))
'''
### constant belonging to Delta c1..3 U
'''
T11 = Kmu*Kmv*Vmu*(3*r1 + r2 + r3)/(60*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))
T12 = Kmu*Kmv*Vmu*(2*r1 + 2*r2 + r3)/(120*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))
T13 = Kmu*Kmv*Vmu*(2*r1 + r2 + 2*r3)/(120*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))

T21 = Kmu*Kmv*Vmu*(2*r1 + 2*r2 + r3)/(120*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))
T31 = Kmu*Kmv*Vmu*(2*r1 + r2 + 2*r3)/(120*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))

c = sp.simplify(T11/(3*r1 + r2 + r3))/2
print('result: \n', sp.simplify(c))
print('lets test: ', sp.simplify(T31/c) ) 
print('lets test: ', sp.simplify(T21/c) ) 
trial = (C_a_v +Kmv)*(C_a_u + Kmu)**2
print('trial: ', sp.simplify(trial*c))
'''
### constant belonging to Delta cM1..3 U
'''
T11 = -C_a_u*Kmv*Vmu*(3*r1 + r2 + r3)/(60*C_a_u*C_a_v**2 + 120*C_a_u*C_a_v*Kmv + 60*C_a_u*Kmv**2 + 60*C_a_v**2*Kmu + 120*C_a_v*Kmu*Kmv + 60*Kmu*Kmv**2)
T12 = -C_a_u*Kmv*Vmu*(2*r1 + 2*r2 + r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)
T13 = -C_a_u*Kmv*Vmu*(2*r1 + r2 + 2*r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)


T31 = -C_a_u*Kmv*Vmu*(2*r1 + r2 + 2*r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)
c = sp.simplify(T11/(3*r1 + r2 + r3))/2
print('result: \n', sp.simplify(c))
print('lets test: ', sp.simplify(T31/c) ) 
trial = (C_a_u +Kmu)*(C_a_v + Kmv)**2
print('trial: ', sp.simplify(trial*c))
'''
### constant belonging to F_v
########################################
########################################




#determining constant
#c = sp.simplify(T1/(2*r1 + r2 + r3))
#print('result: \n', sp.simplify(c))


# splitting constant in numerator and denominator
#num = 3*C_a_u**4*C_a_v*Kmv*Vmu*rq + C_a_u**4*Kmv**2*Vmu*rq + 2*C_a_u**3*C_a_v**2*Kmfu*Vmfv + 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmfv + 6*C_a_u**3*C_a_v*Kmfu*Kmv*Vmu*rq + C_a_u**3*C_a_v*Kmu*Kmv*Vmu*rq + 2*C_a_u**3*Kmfu*Kmv**2*Vmfv + 2*C_a_u**3*Kmfu*Kmv**2*Vmu*rq - C_a_u**3*Kmu*Kmv**2*Vmu*rq + C_a_u**2*C_a_v**2*Kmfu**2*Vmfv + 4*C_a_u**2*C_a_v**2*Kmfu*Kmu*Vmfv + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmfv + 3*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmu*rq + 8*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmfv + 2*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmu*rq + C_a_u**2*Kmfu**2*Kmv**2*Vmfv + C_a_u**2*Kmfu**2*Kmv**2*Vmu*rq + 4*C_a_u**2*Kmfu*Kmu*Kmv**2*Vmfv - 2*C_a_u**2*Kmfu*Kmu*Kmv**2*Vmu*rq + 2*C_a_u*C_a_v**2*Kmfu**2*Kmu*Vmfv + 2*C_a_u*C_a_v**2*Kmfu*Kmu**2*Vmfv + 4*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmfv + C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmu*rq + 4*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv*Vmfv + 2*C_a_u*Kmfu**2*Kmu*Kmv**2*Vmfv - C_a_u*Kmfu**2*Kmu*Kmv**2*Vmu*rq + 2*C_a_u*Kmfu*Kmu**2*Kmv**2*Vmfv + C_a_v**2*Kmfu**2*Kmu**2*Vmfv + 2*C_a_v*Kmfu**2*Kmu**2*Kmv*Vmfv + Kmfu**2*Kmu**2*Kmv**2*Vmfv

#denom = ( C_a_u**4*C_a_v**2 + 2*C_a_u**4*C_a_v*Kmv + C_a_u**4*Kmv**2 + 2*C_a_u**3*C_a_v**2*Kmfu + 2*C_a_u**3*C_a_v**2*Kmu + 4*C_a_u**3*C_a_v*Kmfu*Kmv + 4*C_a_u**3*C_a_v*Kmu*Kmv + 2*C_a_u**3*Kmfu*Kmv**2 + 2*C_a_u**3*Kmu*Kmv**2 + C_a_u**2*C_a_v**2*Kmfu**2 + 4*C_a_u**2*C_a_v**2*Kmfu*Kmu + C_a_u**2*C_a_v**2*Kmu**2 + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv + 8*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv + 2*C_a_u**2*C_a_v*Kmu**2*Kmv + C_a_u**2*Kmfu**2*Kmv**2 + 4*C_a_u**2*Kmfu*Kmu*Kmv**2 + C_a_u**2*Kmu**2*Kmv**2 + 2*C_a_u*C_a_v**2*Kmfu**2*Kmu + 2*C_a_u*C_a_v**2*Kmfu*Kmu**2 + 4*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv + 4*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv + 2*C_a_u*Kmfu**2*Kmu*Kmv**2 + 2*C_a_u*Kmfu*Kmu**2*Kmv**2 + C_a_v**2*Kmfu**2*Kmu**2 + 2*C_a_v*Kmfu**2*Kmu**2*Kmv + Kmfu**2*Kmu**2*Kmv**2) *24



#listing them
#num_l = [3*C_a_u**4*C_a_v*Kmv*Vmu*rq , C_a_u**4*Kmv**2*Vmu*rq , 2*C_a_u**3*C_a_v**2*Kmfu*Vmfv , 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmfv , 6*C_a_u**3*C_a_v*Kmfu*Kmv*Vmu*rq , C_a_u**3*C_a_v*Kmu*Kmv*Vmu*rq , 2*C_a_u**3*Kmfu*Kmv**2*Vmfv , 2*C_a_u**3*Kmfu*Kmv**2*Vmu*rq ,- C_a_u**3*Kmu*Kmv**2*Vmu*rq , C_a_u**2*C_a_v**2*Kmfu**2*Vmfv , 4*C_a_u**2*C_a_v**2*Kmfu*Kmu*Vmfv , 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmfv , 3*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmu*rq , 8*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmfv , 2*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmu*rq , C_a_u**2*Kmfu**2*Kmv**2*Vmfv , C_a_u**2*Kmfu**2*Kmv**2*Vmu*rq , 4*C_a_u**2*Kmfu*Kmu*Kmv**2*Vmfv ,- 2*C_a_u**2*Kmfu*Kmu*Kmv**2*Vmu*rq , 2*C_a_u*C_a_v**2*Kmfu**2*Kmu*Vmfv , 2*C_a_u*C_a_v**2*Kmfu*Kmu**2*Vmfv , 4*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmfv , C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmu*rq , 4*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv*Vmfv , 2*C_a_u*Kmfu**2*Kmu*Kmv**2*Vmfv ,- C_a_u*Kmfu**2*Kmu*Kmv**2*Vmu*rq , 2*C_a_u*Kmfu*Kmu**2*Kmv**2*Vmfv , C_a_v**2*Kmfu**2*Kmu**2*Vmfv , 2*C_a_v*Kmfu**2*Kmu**2*Kmv*Vmfv , Kmfu**2*Kmu**2*Kmv**2*Vmfv]
#denom_l  =[C_a_u**4*C_a_v**2,  2*C_a_u**4*C_a_v*Kmv , C_a_u**4*Kmv**2 , 2*C_a_u**3*C_a_v**2*Kmfu , 2*C_a_u**3*C_a_v**2*Kmu , 4*C_a_u**3*C_a_v*Kmfu*Kmv , 4*C_a_u**3*C_a_v*Kmu*Kmv , 2*C_a_u**3*Kmfu*Kmv**2 , 2*C_a_u**3*Kmu*Kmv**2, C_a_u**2*C_a_v**2*Kmfu**2 , 4*C_a_u**2*C_a_v**2*Kmfu*Kmu , C_a_u**2*C_a_v**2*Kmu**2 , 2*C_a_u**2*C_a_v*Kmfu**2*Kmv , 8*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv , 2*C_a_u**2*C_a_v*Kmu**2*Kmv , C_a_u**2*Kmfu**2*Kmv**2 , 4*C_a_u**2*Kmfu*Kmu*Kmv**2 , C_a_u**2*Kmu**2*Kmv**2 , 2*C_a_u*C_a_v**2*Kmfu**2*Kmu , 2*C_a_u*C_a_v**2*Kmfu*Kmu**2 , 4*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv , 4*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv , 2*C_a_u*Kmfu**2*Kmu*Kmv**2 , 2*C_a_u*Kmfu*Kmu**2*Kmv**2 , C_a_v**2*Kmfu**2*Kmu**2 , 2*C_a_v*Kmfu**2*Kmu**2*Kmv , Kmfu**2*Kmu**2*Kmv**2]


### constant belonging to F_v
########################################
########################################

## calculating
substitutionBest = [ (C_a_u, 8.49), (C_a_v, 0.0163), (Kmv, 27.2438), (Kmu, 0.4103), (Vmu, 0.0031717356019947854), (Vmfv,0.0010016037021638822), (Kmfu, 0.11), (rq,0.95)]
substitutionWorst = [ (C_a_u, 0.9), (C_a_v, 0.31), (Kmv, 27.2438), (Kmu, 0.4103), (Vmu, 0.00021354264845612348) , (Vmfv,0.00014867742939012076), (Kmfu, 0.11), (rq,0.95)]
i = 0
tot_sum_numB = 0
sum_needed_numB = 0
tot_sum_numW = 0
sum_needed_numW = 0
actual_num = []
actual_denom = []
for elem in num_l:
	i = i+1
	ex =  elem 
	B = ex.subs(substitutionBest)
	W = ex.subs(substitutionWorst)
	tot_sum_numB +=abs(B) 
	tot_sum_numW += abs(W)
	if( abs(B)>15 or abs(W) > 0.01): 
		sum_needed_numB +=abs(B)
		sum_needed_numW += abs(W)
		actual_num.append(elem)
		print("Best case: ",i ,B)
		print("Worst case: ",i ,W)
		print("=============================================")

print("current fraction: \n", "Best: ", sum_needed_numB/tot_sum_numB, "Worst: ", sum_needed_numW/tot_sum_numW) 

i = 0
tot_sum_denomB = 0
sum_needed_denomB = 0
tot_sum_denomW = 0
sum_needed_denomW = 0
print("=============================================")
print("=============================================")
print("=============================================")
for elem in denom_l:
	i = i+1
	ex =  elem 
	B = ex.subs(substitutionBest)
	W = ex.subs(substitutionWorst)
	tot_sum_denomB += abs(B) 
	tot_sum_denomW += abs(W) 

	if( abs(B)>10000 or abs(W) > 70): 
		sum_needed_denomB +=abs(B) 
		sum_needed_denomW += abs(W) 
		actual_denom.append(elem)
		print("Best case: ",i ,B)
		print("Worst case: ",i ,W)
		print("=============================================")

print("current fraction: \n", "Best: ", sum_needed_denomB/tot_sum_denomB, "Worst: ", sum_needed_denomW/tot_sum_denomW) 
##
print("num : \n", actual_num, "\ndenom: \n", actual_denom)
# C_a = 9 -> 0.9
# C_v = 0.01 -> 2
# Kmv = 27
# Kmu = 0.4
# Kmfu = 0.1
# Vmu = 3*10-4 -> 3*10-2
# Vmfu = 10-4 -> 10-2

## New approach go over each term and calculate its significance for its best case


#print('who knws: \n', sp.simplify(denom))
#print('lets test: ', sp.simplify(T31/c) ) 
#trial = (C_a_u +Kmu)*(C_a_v + Kmv)**2
#print('trial: ', sp.simplify(trial*c))
































######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
#########################################
### INITIALIZATION DONE #################
#########################################
'''

#  R_U: F- vector for  -e - n + 1  :
FU1 = J*C_a_u*Kmv*Vmu*(4*C_a_u*C_a_v*r1 + 2*C_a_u*C_a_v*r2 + 2*C_a_u*C_a_v*r3 + 2*C_a_u*Kmv*r1 + C_a_u*Kmv*r2 + C_a_u*Kmv*r3 + 2*C_a_v*Kmu*r1 + C_a_v*Kmu*r2 + C_a_v*Kmu*r3)/(24*(C_a_u**2*C_a_v**2 + 2*C_a_u**2*C_a_v*Kmv + C_a_u**2*Kmv**2 + 2*C_a_u*C_a_v**2*Kmu + 4*C_a_u*C_a_v*Kmu*Kmv + 2*C_a_u*Kmu*Kmv**2 + C_a_v**2*Kmu**2 + 2*C_a_v*Kmu**2*Kmv + Kmu**2*Kmv**2))




# R_U: F- vector for  e  :
FU2 = J*C_a_u*Kmv*Vmu*(2*C_a_u*C_a_v*r1 + 4*C_a_u*C_a_v*r2 + 2*C_a_u*C_a_v*r3 + C_a_u*Kmv*r1 + 2*C_a_u*Kmv*r2 + C_a_u*Kmv*r3 + C_a_v*Kmu*r1 + 2*C_a_v*Kmu*r2 + C_a_v*Kmu*r3)/(24*(C_a_u**2*C_a_v**2 + 2*C_a_u**2*C_a_v*Kmv + C_a_u**2*Kmv**2 + 2*C_a_u*C_a_v**2*Kmu + 4*C_a_u*C_a_v*Kmu*Kmv + 2*C_a_u*Kmu*Kmv**2 + C_a_v**2*Kmu**2 + 2*C_a_v*Kmu**2*Kmv + Kmu**2*Kmv**2))

# R_U: F- vector for  n  :
FU3 = J*C_a_u*Kmv*Vmu*(2*C_a_u*C_a_v*r1 + 2*C_a_u*C_a_v*r2 + 4*C_a_u*C_a_v*r3 + C_a_u*Kmv*r1 + C_a_u*Kmv*r2 + 2*C_a_u*Kmv*r3 + C_a_v*Kmu*r1 + C_a_v*Kmu*r2 + 2*C_a_v*Kmu*r3)/(24*(C_a_u**2*C_a_v**2 + 2*C_a_u**2*C_a_v*Kmv + C_a_u**2*Kmv**2 + 2*C_a_u*C_a_v**2*Kmu + 4*C_a_u*C_a_v*Kmu*Kmv + 2*C_a_u*Kmu*Kmv**2 + C_a_v**2*Kmu**2 + 2*C_a_v*Kmu**2*Kmv + Kmu**2*Kmv**2))

#const = sp.simplify(sp.simplify(FU1)/ (2*r1 + r2 + r3) )
#print("constant: ", const)
#print("F1: ", sp.simplify(FU1/const) )
#print("F2: ", sp.simplify(FU2/const) )
#print("F3: ", sp.simplify(FU3/const) )



#  R_U: F- vector for  -e - n + 1  :
FV1 = J*(4*C_a_u**4*C_a_v*Kmv*Vmu*r1*rq + 2*C_a_u**4*C_a_v*Kmv*Vmu*r2*rq + 2*C_a_u**4*C_a_v*Kmv*Vmu*r3*rq + 2*C_a_u**4*Kmv**2*Vmu*r1*rq + C_a_u**4*Kmv**2*Vmu*r2*rq + C_a_u**4*Kmv**2*Vmu*r3*rq + 4*C_a_u**3*C_a_v**2*Kmfu*Vmfv*r1 + 2*C_a_u**3*C_a_v**2*Kmfu*Vmfv*r2 + 2*C_a_u**3*C_a_v**2*Kmfu*Vmfv*r3 + 8*C_a_u**3*C_a_v*Kmfu*Kmv*Vmfv*r1 + 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmfv*r2 + 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmfv*r3 + 8*C_a_u**3*C_a_v*Kmfu*Kmv*Vmu*r1*rq + 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmu*r2*rq + 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmu*r3*rq + 2*C_a_u**3*C_a_v*Kmu*Kmv*Vmu*r1*rq + C_a_u**3*C_a_v*Kmu*Kmv*Vmu*r2*rq + C_a_u**3*C_a_v*Kmu*Kmv*Vmu*r3*rq + 4*C_a_u**3*Kmfu*Kmv**2*Vmfv*r1 + 2*C_a_u**3*Kmfu*Kmv**2*Vmfv*r2 + 2*C_a_u**3*Kmfu*Kmv**2*Vmfv*r3 + 4*C_a_u**3*Kmfu*Kmv**2*Vmu*r1*rq + 2*C_a_u**3*Kmfu*Kmv**2*Vmu*r2*rq + 2*C_a_u**3*Kmfu*Kmv**2*Vmu*r3*rq + 2*C_a_u**2*C_a_v**2*Kmfu**2*Vmfv*r1 + C_a_u**2*C_a_v**2*Kmfu**2*Vmfv*r2 + C_a_u**2*C_a_v**2*Kmfu**2*Vmfv*r3 + 8*C_a_u**2*C_a_v**2*Kmfu*Kmu*Vmfv*r1 + 4*C_a_u**2*C_a_v**2*Kmfu*Kmu*Vmfv*r2 + 4*C_a_u**2*C_a_v**2*Kmfu*Kmu*Vmfv*r3 + 4*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmfv*r1 + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmfv*r2 + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmfv*r3 + 4*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmu*r1*rq + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmu*r2*rq + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmu*r3*rq + 16*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmfv*r1 + 8*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmfv*r2 + 8*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmfv*r3 + 4*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmu*r1*rq + 2*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmu*r2*rq + 2*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmu*r3*rq + 2*C_a_u**2*Kmfu**2*Kmv**2*Vmfv*r1 + C_a_u**2*Kmfu**2*Kmv**2*Vmfv*r2 + C_a_u**2*Kmfu**2*Kmv**2*Vmfv*r3 + 2*C_a_u**2*Kmfu**2*Kmv**2*Vmu*r1*rq + C_a_u**2*Kmfu**2*Kmv**2*Vmu*r2*rq + C_a_u**2*Kmfu**2*Kmv**2*Vmu*r3*rq + 8*C_a_u**2*Kmfu*Kmu*Kmv**2*Vmfv*r1 + 4*C_a_u**2*Kmfu*Kmu*Kmv**2*Vmfv*r2 + 4*C_a_u**2*Kmfu*Kmu*Kmv**2*Vmfv*r3 + 4*C_a_u*C_a_v**2*Kmfu**2*Kmu*Vmfv*r1 + 2*C_a_u*C_a_v**2*Kmfu**2*Kmu*Vmfv*r2 + 2*C_a_u*C_a_v**2*Kmfu**2*Kmu*Vmfv*r3 + 4*C_a_u*C_a_v**2*Kmfu*Kmu**2*Vmfv*r1 + 2*C_a_u*C_a_v**2*Kmfu*Kmu**2*Vmfv*r2 + 2*C_a_u*C_a_v**2*Kmfu*Kmu**2*Vmfv*r3 + 8*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmfv*r1 + 4*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmfv*r2 + 4*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmfv*r3 + 2*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmu*r1*rq + C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmu*r2*rq + C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmu*r3*rq + 8*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv*Vmfv*r1 + 4*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv*Vmfv*r2 + 4*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv*Vmfv*r3 + 4*C_a_u*Kmfu**2*Kmu*Kmv**2*Vmfv*r1 + 2*C_a_u*Kmfu**2*Kmu*Kmv**2*Vmfv*r2 + 2*C_a_u*Kmfu**2*Kmu*Kmv**2*Vmfv*r3 + 4*C_a_u*Kmfu*Kmu**2*Kmv**2*Vmfv*r1 + 2*C_a_u*Kmfu*Kmu**2*Kmv**2*Vmfv*r2 + 2*C_a_u*Kmfu*Kmu**2*Kmv**2*Vmfv*r3 + 2*C_a_v**2*Kmfu**2*Kmu**2*Vmfv*r1 + C_a_v**2*Kmfu**2*Kmu**2*Vmfv*r2 + C_a_v**2*Kmfu**2*Kmu**2*Vmfv*r3 + 4*C_a_v*Kmfu**2*Kmu**2*Kmv*Vmfv*r1 + 2*C_a_v*Kmfu**2*Kmu**2*Kmv*Vmfv*r2 + 2*C_a_v*Kmfu**2*Kmu**2*Kmv*Vmfv*r3 + 2*Kmfu**2*Kmu**2*Kmv**2*Vmfv*r1 + Kmfu**2*Kmu**2*Kmv**2*Vmfv*r2 + Kmfu**2*Kmu**2*Kmv**2*Vmfv*r3)/(24*(C_a_u**4*C_a_v**2 + 2*C_a_u**4*C_a_v*Kmv + C_a_u**4*Kmv**2 + 2*C_a_u**3*C_a_v**2*Kmfu + 2*C_a_u**3*C_a_v**2*Kmu + 4*C_a_u**3*C_a_v*Kmfu*Kmv + 4*C_a_u**3*C_a_v*Kmu*Kmv + 2*C_a_u**3*Kmfu*Kmv**2 + 2*C_a_u**3*Kmu*Kmv**2 + C_a_u**2*C_a_v**2*Kmfu**2 + 4*C_a_u**2*C_a_v**2*Kmfu*Kmu + C_a_u**2*C_a_v**2*Kmu**2 + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv + 8*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv + 2*C_a_u**2*C_a_v*Kmu**2*Kmv + C_a_u**2*Kmfu**2*Kmv**2 + 4*C_a_u**2*Kmfu*Kmu*Kmv**2 + C_a_u**2*Kmu**2*Kmv**2 + 2*C_a_u*C_a_v**2*Kmfu**2*Kmu + 2*C_a_u*C_a_v**2*Kmfu*Kmu**2 + 4*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv + 4*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv + 2*C_a_u*Kmfu**2*Kmu*Kmv**2 + 2*C_a_u*Kmfu*Kmu**2*Kmv**2 + C_a_v**2*Kmfu**2*Kmu**2 + 2*C_a_v*Kmfu**2*Kmu**2*Kmv + Kmfu**2*Kmu**2*Kmv**2))



# R_V: F- vector for  e  :
FV2 = J*(2*C_a_u**4*C_a_v*Kmv*Vmu*r1*rq + 4*C_a_u**4*C_a_v*Kmv*Vmu*r2*rq + 2*C_a_u**4*C_a_v*Kmv*Vmu*r3*rq + C_a_u**4*Kmv**2*Vmu*r1*rq + 2*C_a_u**4*Kmv**2*Vmu*r2*rq + C_a_u**4*Kmv**2*Vmu*r3*rq + 2*C_a_u**3*C_a_v**2*Kmfu*Vmfv*r1 + 4*C_a_u**3*C_a_v**2*Kmfu*Vmfv*r2 + 2*C_a_u**3*C_a_v**2*Kmfu*Vmfv*r3 + 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmfv*r1 + 8*C_a_u**3*C_a_v*Kmfu*Kmv*Vmfv*r2 + 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmfv*r3 + 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmu*r1*rq + 8*C_a_u**3*C_a_v*Kmfu*Kmv*Vmu*r2*rq + 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmu*r3*rq + C_a_u**3*C_a_v*Kmu*Kmv*Vmu*r1*rq + 2*C_a_u**3*C_a_v*Kmu*Kmv*Vmu*r2*rq + C_a_u**3*C_a_v*Kmu*Kmv*Vmu*r3*rq + 2*C_a_u**3*Kmfu*Kmv**2*Vmfv*r1 + 4*C_a_u**3*Kmfu*Kmv**2*Vmfv*r2 + 2*C_a_u**3*Kmfu*Kmv**2*Vmfv*r3 + 2*C_a_u**3*Kmfu*Kmv**2*Vmu*r1*rq + 4*C_a_u**3*Kmfu*Kmv**2*Vmu*r2*rq + 2*C_a_u**3*Kmfu*Kmv**2*Vmu*r3*rq + C_a_u**2*C_a_v**2*Kmfu**2*Vmfv*r1 + 2*C_a_u**2*C_a_v**2*Kmfu**2*Vmfv*r2 + C_a_u**2*C_a_v**2*Kmfu**2*Vmfv*r3 + 4*C_a_u**2*C_a_v**2*Kmfu*Kmu*Vmfv*r1 + 8*C_a_u**2*C_a_v**2*Kmfu*Kmu*Vmfv*r2 + 4*C_a_u**2*C_a_v**2*Kmfu*Kmu*Vmfv*r3 + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmfv*r1 + 4*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmfv*r2 + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmfv*r3 + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmu*r1*rq + 4*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmu*r2*rq + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmu*r3*rq + 8*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmfv*r1 + 16*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmfv*r2 + 8*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmfv*r3 + 2*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmu*r1*rq + 4*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmu*r2*rq + 2*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmu*r3*rq + C_a_u**2*Kmfu**2*Kmv**2*Vmfv*r1 + 2*C_a_u**2*Kmfu**2*Kmv**2*Vmfv*r2 + C_a_u**2*Kmfu**2*Kmv**2*Vmfv*r3 + C_a_u**2*Kmfu**2*Kmv**2*Vmu*r1*rq + 2*C_a_u**2*Kmfu**2*Kmv**2*Vmu*r2*rq + C_a_u**2*Kmfu**2*Kmv**2*Vmu*r3*rq + 4*C_a_u**2*Kmfu*Kmu*Kmv**2*Vmfv*r1 + 8*C_a_u**2*Kmfu*Kmu*Kmv**2*Vmfv*r2 + 4*C_a_u**2*Kmfu*Kmu*Kmv**2*Vmfv*r3 + 2*C_a_u*C_a_v**2*Kmfu**2*Kmu*Vmfv*r1 + 4*C_a_u*C_a_v**2*Kmfu**2*Kmu*Vmfv*r2 + 2*C_a_u*C_a_v**2*Kmfu**2*Kmu*Vmfv*r3 + 2*C_a_u*C_a_v**2*Kmfu*Kmu**2*Vmfv*r1 + 4*C_a_u*C_a_v**2*Kmfu*Kmu**2*Vmfv*r2 + 2*C_a_u*C_a_v**2*Kmfu*Kmu**2*Vmfv*r3 + 4*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmfv*r1 + 8*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmfv*r2 + 4*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmfv*r3 + C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmu*r1*rq + 2*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmu*r2*rq + C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmu*r3*rq + 4*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv*Vmfv*r1 + 8*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv*Vmfv*r2 + 4*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv*Vmfv*r3 + 2*C_a_u*Kmfu**2*Kmu*Kmv**2*Vmfv*r1 + 4*C_a_u*Kmfu**2*Kmu*Kmv**2*Vmfv*r2 + 2*C_a_u*Kmfu**2*Kmu*Kmv**2*Vmfv*r3 + 2*C_a_u*Kmfu*Kmu**2*Kmv**2*Vmfv*r1 + 4*C_a_u*Kmfu*Kmu**2*Kmv**2*Vmfv*r2 + 2*C_a_u*Kmfu*Kmu**2*Kmv**2*Vmfv*r3 + C_a_v**2*Kmfu**2*Kmu**2*Vmfv*r1 + 2*C_a_v**2*Kmfu**2*Kmu**2*Vmfv*r2 + C_a_v**2*Kmfu**2*Kmu**2*Vmfv*r3 + 2*C_a_v*Kmfu**2*Kmu**2*Kmv*Vmfv*r1 + 4*C_a_v*Kmfu**2*Kmu**2*Kmv*Vmfv*r2 + 2*C_a_v*Kmfu**2*Kmu**2*Kmv*Vmfv*r3 + Kmfu**2*Kmu**2*Kmv**2*Vmfv*r1 + 2*Kmfu**2*Kmu**2*Kmv**2*Vmfv*r2 + Kmfu**2*Kmu**2*Kmv**2*Vmfv*r3)/(24*(C_a_u**4*C_a_v**2 + 2*C_a_u**4*C_a_v*Kmv + C_a_u**4*Kmv**2 + 2*C_a_u**3*C_a_v**2*Kmfu + 2*C_a_u**3*C_a_v**2*Kmu + 4*C_a_u**3*C_a_v*Kmfu*Kmv + 4*C_a_u**3*C_a_v*Kmu*Kmv + 2*C_a_u**3*Kmfu*Kmv**2 + 2*C_a_u**3*Kmu*Kmv**2 + C_a_u**2*C_a_v**2*Kmfu**2 + 4*C_a_u**2*C_a_v**2*Kmfu*Kmu + C_a_u**2*C_a_v**2*Kmu**2 + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv + 8*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv + 2*C_a_u**2*C_a_v*Kmu**2*Kmv + C_a_u**2*Kmfu**2*Kmv**2 + 4*C_a_u**2*Kmfu*Kmu*Kmv**2 + C_a_u**2*Kmu**2*Kmv**2 + 2*C_a_u*C_a_v**2*Kmfu**2*Kmu + 2*C_a_u*C_a_v**2*Kmfu*Kmu**2 + 4*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv + 4*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv + 2*C_a_u*Kmfu**2*Kmu*Kmv**2 + 2*C_a_u*Kmfu*Kmu**2*Kmv**2 + C_a_v**2*Kmfu**2*Kmu**2 + 2*C_a_v*Kmfu**2*Kmu**2*Kmv + Kmfu**2*Kmu**2*Kmv**2))


# R_V: F- vector for  n  :
FV3 = J*(2*C_a_u**4*C_a_v*Kmv*Vmu*r1*rq + 2*C_a_u**4*C_a_v*Kmv*Vmu*r2*rq + 4*C_a_u**4*C_a_v*Kmv*Vmu*r3*rq + C_a_u**4*Kmv**2*Vmu*r1*rq + C_a_u**4*Kmv**2*Vmu*r2*rq + 2*C_a_u**4*Kmv**2*Vmu*r3*rq + 2*C_a_u**3*C_a_v**2*Kmfu*Vmfv*r1 + 2*C_a_u**3*C_a_v**2*Kmfu*Vmfv*r2 + 4*C_a_u**3*C_a_v**2*Kmfu*Vmfv*r3 + 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmfv*r1 + 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmfv*r2 + 8*C_a_u**3*C_a_v*Kmfu*Kmv*Vmfv*r3 + 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmu*r1*rq + 4*C_a_u**3*C_a_v*Kmfu*Kmv*Vmu*r2*rq + 8*C_a_u**3*C_a_v*Kmfu*Kmv*Vmu*r3*rq + C_a_u**3*C_a_v*Kmu*Kmv*Vmu*r1*rq + C_a_u**3*C_a_v*Kmu*Kmv*Vmu*r2*rq + 2*C_a_u**3*C_a_v*Kmu*Kmv*Vmu*r3*rq + 2*C_a_u**3*Kmfu*Kmv**2*Vmfv*r1 + 2*C_a_u**3*Kmfu*Kmv**2*Vmfv*r2 + 4*C_a_u**3*Kmfu*Kmv**2*Vmfv*r3 + 2*C_a_u**3*Kmfu*Kmv**2*Vmu*r1*rq + 2*C_a_u**3*Kmfu*Kmv**2*Vmu*r2*rq + 4*C_a_u**3*Kmfu*Kmv**2*Vmu*r3*rq + C_a_u**2*C_a_v**2*Kmfu**2*Vmfv*r1 + C_a_u**2*C_a_v**2*Kmfu**2*Vmfv*r2 + 2*C_a_u**2*C_a_v**2*Kmfu**2*Vmfv*r3 + 4*C_a_u**2*C_a_v**2*Kmfu*Kmu*Vmfv*r1 + 4*C_a_u**2*C_a_v**2*Kmfu*Kmu*Vmfv*r2 + 8*C_a_u**2*C_a_v**2*Kmfu*Kmu*Vmfv*r3 + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmfv*r1 + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmfv*r2 + 4*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmfv*r3 + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmu*r1*rq + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmu*r2*rq + 4*C_a_u**2*C_a_v*Kmfu**2*Kmv*Vmu*r3*rq + 8*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmfv*r1 + 8*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmfv*r2 + 16*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmfv*r3 + 2*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmu*r1*rq + 2*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmu*r2*rq + 4*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv*Vmu*r3*rq + C_a_u**2*Kmfu**2*Kmv**2*Vmfv*r1 + C_a_u**2*Kmfu**2*Kmv**2*Vmfv*r2 + 2*C_a_u**2*Kmfu**2*Kmv**2*Vmfv*r3 + C_a_u**2*Kmfu**2*Kmv**2*Vmu*r1*rq + C_a_u**2*Kmfu**2*Kmv**2*Vmu*r2*rq + 2*C_a_u**2*Kmfu**2*Kmv**2*Vmu*r3*rq + 4*C_a_u**2*Kmfu*Kmu*Kmv**2*Vmfv*r1 + 4*C_a_u**2*Kmfu*Kmu*Kmv**2*Vmfv*r2 + 8*C_a_u**2*Kmfu*Kmu*Kmv**2*Vmfv*r3 + 2*C_a_u*C_a_v**2*Kmfu**2*Kmu*Vmfv*r1 + 2*C_a_u*C_a_v**2*Kmfu**2*Kmu*Vmfv*r2 + 4*C_a_u*C_a_v**2*Kmfu**2*Kmu*Vmfv*r3 + 2*C_a_u*C_a_v**2*Kmfu*Kmu**2*Vmfv*r1 + 2*C_a_u*C_a_v**2*Kmfu*Kmu**2*Vmfv*r2 + 4*C_a_u*C_a_v**2*Kmfu*Kmu**2*Vmfv*r3 + 4*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmfv*r1 + 4*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmfv*r2 + 8*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmfv*r3 + C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmu*r1*rq + C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmu*r2*rq + 2*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv*Vmu*r3*rq + 4*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv*Vmfv*r1 + 4*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv*Vmfv*r2 + 8*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv*Vmfv*r3 + 2*C_a_u*Kmfu**2*Kmu*Kmv**2*Vmfv*r1 + 2*C_a_u*Kmfu**2*Kmu*Kmv**2*Vmfv*r2 + 4*C_a_u*Kmfu**2*Kmu*Kmv**2*Vmfv*r3 + 2*C_a_u*Kmfu*Kmu**2*Kmv**2*Vmfv*r1 + 2*C_a_u*Kmfu*Kmu**2*Kmv**2*Vmfv*r2 + 4*C_a_u*Kmfu*Kmu**2*Kmv**2*Vmfv*r3 + C_a_v**2*Kmfu**2*Kmu**2*Vmfv*r1 + C_a_v**2*Kmfu**2*Kmu**2*Vmfv*r2 + 2*C_a_v**2*Kmfu**2*Kmu**2*Vmfv*r3 + 2*C_a_v*Kmfu**2*Kmu**2*Kmv*Vmfv*r1 + 2*C_a_v*Kmfu**2*Kmu**2*Kmv*Vmfv*r2 + 4*C_a_v*Kmfu**2*Kmu**2*Kmv*Vmfv*r3 + Kmfu**2*Kmu**2*Kmv**2*Vmfv*r1 + Kmfu**2*Kmu**2*Kmv**2*Vmfv*r2 + 2*Kmfu**2*Kmu**2*Kmv**2*Vmfv*r3)/(24*(C_a_u**4*C_a_v**2 + 2*C_a_u**4*C_a_v*Kmv + C_a_u**4*Kmv**2 + 2*C_a_u**3*C_a_v**2*Kmfu + 2*C_a_u**3*C_a_v**2*Kmu + 4*C_a_u**3*C_a_v*Kmfu*Kmv + 4*C_a_u**3*C_a_v*Kmu*Kmv + 2*C_a_u**3*Kmfu*Kmv**2 + 2*C_a_u**3*Kmu*Kmv**2 + C_a_u**2*C_a_v**2*Kmfu**2 + 4*C_a_u**2*C_a_v**2*Kmfu*Kmu + C_a_u**2*C_a_v**2*Kmu**2 + 2*C_a_u**2*C_a_v*Kmfu**2*Kmv + 8*C_a_u**2*C_a_v*Kmfu*Kmu*Kmv + 2*C_a_u**2*C_a_v*Kmu**2*Kmv + C_a_u**2*Kmfu**2*Kmv**2 + 4*C_a_u**2*Kmfu*Kmu*Kmv**2 + C_a_u**2*Kmu**2*Kmv**2 + 2*C_a_u*C_a_v**2*Kmfu**2*Kmu + 2*C_a_u*C_a_v**2*Kmfu*Kmu**2 + 4*C_a_u*C_a_v*Kmfu**2*Kmu*Kmv + 4*C_a_u*C_a_v*Kmfu*Kmu**2*Kmv + 2*C_a_u*Kmfu**2*Kmu*Kmv**2 + 2*C_a_u*Kmfu*Kmu**2*Kmv**2 + C_a_v**2*Kmfu**2*Kmu**2 + 2*C_a_v*Kmfu**2*Kmu**2*Kmv + Kmfu**2*Kmu**2*Kmv**2))


const = sp.simplify( sp.simplify(sp.simplify(FV1)/ (2*r1 + r2 + r3) ))
#const = sp.simplify(sp.simplify(FV1)/ sp.simplify(sp.simplify(FV2)))
print("constant: ", const)
print("F1: ", sp.simplify(FV1/const) )
print("F2: ", sp.simplify(FV2/const) )
print("F3: ", sp.simplify(FV3/const) )



##
# R_U: Delta_c1..c3..cM1..cM3 vectors for  -e - n + 1  :

D1FU1 = Kmu*Kmv*Vmu*J*(3*r1 + r2 + r3)/(60*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))
D2FU1 = Kmu*Kmv*Vmu*J*(2*r1 + 2*r2 + r3)/(120*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))
D3FU1 = Kmu*Kmv*Vmu*J*(2*r1 + r2 + 2*r3)/(120*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))
DM1FU1 = -C_a_u*Kmv*Vmu*J*(3*r1 + r2 + r3)/(60*C_a_u*C_a_v**2 + 120*C_a_u*C_a_v*Kmv + 60*C_a_u*Kmv**2 + 60*C_a_v**2*Kmu + 120*C_a_v*Kmu*Kmv + 60*Kmu*Kmv**2)
DM2FU1 = -C_a_u*Kmv*Vmu*J*(2*r1 + 2*r2 + r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)
DM3FU1 = -C_a_u*Kmv*Vmu*J*(2*r1 + r2 + 2*r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)

# R_U: Delta_c1..c3..cM1..cM3 vectors for  e:
D1FU2 = Kmu*Kmv*Vmu*J*(2*r1 + 2*r2 + r3)/(120*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))
D2FU2 = Kmu*Kmv*Vmu*J*(r1 + 3*r2 + r3)/(60*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))
D3FU2 = Kmu*Kmv*Vmu*J*(r1 + 2*r2 + 2*r3)/(120*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))
DM1FU2 = -C_a_u*Kmv*Vmu*J*(2*r1 + 2*r2 + r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)
DM2FU2 = -C_a_u*Kmv*Vmu*J*(r1 + 3*r2 + r3)/(60*C_a_u*C_a_v**2 + 120*C_a_u*C_a_v*Kmv + 60*C_a_u*Kmv**2 + 60*C_a_v**2*Kmu + 120*C_a_v*Kmu*Kmv + 60*Kmu*Kmv**2)
DM3FU2 = -C_a_u*Kmv*Vmu*J*(r1 + 2*r2 + 2*r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)

# R_U: Delta_c1..c3..cM1..cM3 vectors for  n:
D1FU3 = Kmu*Kmv*Vmu*J*(2*r1 + r2 + 2*r3)/(120*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))
D2FU3 = Kmu*Kmv*Vmu*J*(r1 + 2*r2 + 2*r3)/(120*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))
D3FU3 = Kmu*Kmv*Vmu*J*(r1 + r2 + 3*r3)/(60*(C_a_u**2*C_a_v + C_a_u**2*Kmv + 2*C_a_u*C_a_v*Kmu + 2*C_a_u*Kmu*Kmv + C_a_v*Kmu**2 + Kmu**2*Kmv))
DM1FU3 = -C_a_u*Kmv*Vmu*J*(2*r1 + r2 + 2*r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)
DM2FU3 = -C_a_u*Kmv*Vmu*J*(r1 + 2*r2 + 2*r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)
DM3FU3 = -C_a_u*Kmv*Vmu*J*(r1 + r2 + 3*r3)/(60*C_a_u*C_a_v**2 + 120*C_a_u*C_a_v*Kmv + 60*C_a_u*Kmv**2 + 60*C_a_v**2*Kmu + 120*C_a_v*Kmu*Kmv + 60*Kmu*Kmv**2)


##
# R_V: Delta_c1..c3..cM1..cM3 vectors for  -e - n + 1  :

D1FV1 = -J*(3*C_a_u**2*C_a_v*Kmfu*Vmfv*r1 + C_a_u**2*C_a_v*Kmfu*Vmfv*r2 + C_a_u**2*C_a_v*Kmfu*Vmfv*r3 + 3*C_a_u**2*Kmfu*Kmv*Vmfv*r1 + C_a_u**2*Kmfu*Kmv*Vmfv*r2 + C_a_u**2*Kmfu*Kmv*Vmfv*r3 - 3*C_a_u**2*Kmu*Kmv*Vmu*r1*rq - C_a_u**2*Kmu*Kmv*Vmu*r2*rq - C_a_u**2*Kmu*Kmv*Vmu*r3*rq + 6*C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r1 + 2*C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r2 + 2*C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r3 + 6*C_a_u*Kmfu*Kmu*Kmv*Vmfv*r1 + 2*C_a_u*Kmfu*Kmu*Kmv*Vmfv*r2 + 2*C_a_u*Kmfu*Kmu*Kmv*Vmfv*r3 - 6*C_a_u*Kmfu*Kmu*Kmv*Vmu*r1*rq - 2*C_a_u*Kmfu*Kmu*Kmv*Vmu*r2*rq - 2*C_a_u*Kmfu*Kmu*Kmv*Vmu*r3*rq + 3*C_a_v*Kmfu*Kmu**2*Vmfv*r1 + C_a_v*Kmfu*Kmu**2*Vmfv*r2 + C_a_v*Kmfu*Kmu**2*Vmfv*r3 - 3*Kmfu**2*Kmu*Kmv*Vmu*r1*rq - Kmfu**2*Kmu*Kmv*Vmu*r2*rq - Kmfu**2*Kmu*Kmv*Vmu*r3*rq + 3*Kmfu*Kmu**2*Kmv*Vmfv*r1 + Kmfu*Kmu**2*Kmv*Vmfv*r2 + Kmfu*Kmu**2*Kmv*Vmfv*r3)/(60*C_a_u**4*C_a_v + 60*C_a_u**4*Kmv + 120*C_a_u**3*C_a_v*Kmfu + 120*C_a_u**3*C_a_v*Kmu + 120*C_a_u**3*Kmfu*Kmv + 120*C_a_u**3*Kmu*Kmv + 60*C_a_u**2*C_a_v*Kmfu**2 + 240*C_a_u**2*C_a_v*Kmfu*Kmu + 60*C_a_u**2*C_a_v*Kmu**2 + 60*C_a_u**2*Kmfu**2*Kmv + 240*C_a_u**2*Kmfu*Kmu*Kmv + 60*C_a_u**2*Kmu**2*Kmv + 120*C_a_u*C_a_v*Kmfu**2*Kmu + 120*C_a_u*C_a_v*Kmfu*Kmu**2 + 120*C_a_u*Kmfu**2*Kmu*Kmv + 120*C_a_u*Kmfu*Kmu**2*Kmv + 60*C_a_v*Kmfu**2*Kmu**2 + 60*Kmfu**2*Kmu**2*Kmv)
D2FV1 = J*(-C_a_u**2*C_a_v*Kmfu*Vmfv*r1/60 - C_a_u**2*C_a_v*Kmfu*Vmfv*r2/60 - C_a_u**2*C_a_v*Kmfu*Vmfv*r3/120 - C_a_u**2*Kmfu*Kmv*Vmfv*r1/60 - C_a_u**2*Kmfu*Kmv*Vmfv*r2/60 - C_a_u**2*Kmfu*Kmv*Vmfv*r3/120 + C_a_u**2*Kmu*Kmv*Vmu*r1*rq/60 + C_a_u**2*Kmu*Kmv*Vmu*r2*rq/60 + C_a_u**2*Kmu*Kmv*Vmu*r3*rq/120 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r1/30 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r2/30 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r3/60 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r1/30 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r2/30 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r3/60 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r1*rq/30 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r2*rq/30 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r3*rq/60 - C_a_v*Kmfu*Kmu**2*Vmfv*r1/60 - C_a_v*Kmfu*Kmu**2*Vmfv*r2/60 - C_a_v*Kmfu*Kmu**2*Vmfv*r3/120 + Kmfu**2*Kmu*Kmv*Vmu*r1*rq/60 + Kmfu**2*Kmu*Kmv*Vmu*r2*rq/60 + Kmfu**2*Kmu*Kmv*Vmu*r3*rq/120 - Kmfu*Kmu**2*Kmv*Vmfv*r1/60 - Kmfu*Kmu**2*Kmv*Vmfv*r2/60 - Kmfu*Kmu**2*Kmv*Vmfv*r3/120)/(C_a_u**4*C_a_v + C_a_u**4*Kmv + 2*C_a_u**3*C_a_v*Kmfu + 2*C_a_u**3*C_a_v*Kmu + 2*C_a_u**3*Kmfu*Kmv + 2*C_a_u**3*Kmu*Kmv + C_a_u**2*C_a_v*Kmfu**2 + 4*C_a_u**2*C_a_v*Kmfu*Kmu + C_a_u**2*C_a_v*Kmu**2 + C_a_u**2*Kmfu**2*Kmv + 4*C_a_u**2*Kmfu*Kmu*Kmv + C_a_u**2*Kmu**2*Kmv + 2*C_a_u*C_a_v*Kmfu**2*Kmu + 2*C_a_u*C_a_v*Kmfu*Kmu**2 + 2*C_a_u*Kmfu**2*Kmu*Kmv + 2*C_a_u*Kmfu*Kmu**2*Kmv + C_a_v*Kmfu**2*Kmu**2 + Kmfu**2*Kmu**2*Kmv)
D3FV1 = J*(-C_a_u**2*C_a_v*Kmfu*Vmfv*r1/60 - C_a_u**2*C_a_v*Kmfu*Vmfv*r2/120 - C_a_u**2*C_a_v*Kmfu*Vmfv*r3/60 - C_a_u**2*Kmfu*Kmv*Vmfv*r1/60 - C_a_u**2*Kmfu*Kmv*Vmfv*r2/120 - C_a_u**2*Kmfu*Kmv*Vmfv*r3/60 + C_a_u**2*Kmu*Kmv*Vmu*r1*rq/60 + C_a_u**2*Kmu*Kmv*Vmu*r2*rq/120 + C_a_u**2*Kmu*Kmv*Vmu*r3*rq/60 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r1/30 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r2/60 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r3/30 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r1/30 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r2/60 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r3/30 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r1*rq/30 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r2*rq/60 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r3*rq/30 - C_a_v*Kmfu*Kmu**2*Vmfv*r1/60 - C_a_v*Kmfu*Kmu**2*Vmfv*r2/120 - C_a_v*Kmfu*Kmu**2*Vmfv*r3/60 + Kmfu**2*Kmu*Kmv*Vmu*r1*rq/60 + Kmfu**2*Kmu*Kmv*Vmu*r2*rq/120 + Kmfu**2*Kmu*Kmv*Vmu*r3*rq/60 - Kmfu*Kmu**2*Kmv*Vmfv*r1/60 - Kmfu*Kmu**2*Kmv*Vmfv*r2/120 - Kmfu*Kmu**2*Kmv*Vmfv*r3/60)/(C_a_u**4*C_a_v + C_a_u**4*Kmv + 2*C_a_u**3*C_a_v*Kmfu + 2*C_a_u**3*C_a_v*Kmu + 2*C_a_u**3*Kmfu*Kmv + 2*C_a_u**3*Kmu*Kmv + C_a_u**2*C_a_v*Kmfu**2 + 4*C_a_u**2*C_a_v*Kmfu*Kmu + C_a_u**2*C_a_v*Kmu**2 + C_a_u**2*Kmfu**2*Kmv + 4*C_a_u**2*Kmfu*Kmu*Kmv + C_a_u**2*Kmu**2*Kmv + 2*C_a_u*C_a_v*Kmfu**2*Kmu + 2*C_a_u*C_a_v*Kmfu*Kmu**2 + 2*C_a_u*Kmfu**2*Kmu*Kmv + 2*C_a_u*Kmfu*Kmu**2*Kmv + C_a_v*Kmfu**2*Kmu**2 + Kmfu**2*Kmu**2*Kmv)
DM1FV1 = -C_a_u*Kmv*Vmu*rq*J*(3*r1 + r2 + r3)/(60*C_a_u*C_a_v**2 + 120*C_a_u*C_a_v*Kmv + 60*C_a_u*Kmv**2 + 60*C_a_v**2*Kmu + 120*C_a_v*Kmu*Kmv + 60*Kmu*Kmv**2)
DM2FV1 =-C_a_u*Kmv*Vmu*rq*J*(2*r1 + 2*r2 + r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)
DM3FV1 = -C_a_u*Kmv*Vmu*rq*J*(2*r1 + r2 + 2*r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)

# R_V: Delta_c1..c3..cM1..cM3 vectors for  e:

D1FV2 = J*(-C_a_u**2*C_a_v*Kmfu*Vmfv*r1/60 - C_a_u**2*C_a_v*Kmfu*Vmfv*r2/60 - C_a_u**2*C_a_v*Kmfu*Vmfv*r3/120 - C_a_u**2*Kmfu*Kmv*Vmfv*r1/60 - C_a_u**2*Kmfu*Kmv*Vmfv*r2/60 - C_a_u**2*Kmfu*Kmv*Vmfv*r3/120 + C_a_u**2*Kmu*Kmv*Vmu*r1*rq/60 + C_a_u**2*Kmu*Kmv*Vmu*r2*rq/60 + C_a_u**2*Kmu*Kmv*Vmu*r3*rq/120 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r1/30 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r2/30 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r3/60 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r1/30 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r2/30 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r3/60 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r1*rq/30 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r2*rq/30 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r3*rq/60 - C_a_v*Kmfu*Kmu**2*Vmfv*r1/60 - C_a_v*Kmfu*Kmu**2*Vmfv*r2/60 - C_a_v*Kmfu*Kmu**2*Vmfv*r3/120 + Kmfu**2*Kmu*Kmv*Vmu*r1*rq/60 + Kmfu**2*Kmu*Kmv*Vmu*r2*rq/60 + Kmfu**2*Kmu*Kmv*Vmu*r3*rq/120 - Kmfu*Kmu**2*Kmv*Vmfv*r1/60 - Kmfu*Kmu**2*Kmv*Vmfv*r2/60 - Kmfu*Kmu**2*Kmv*Vmfv*r3/120)/(C_a_u**4*C_a_v + C_a_u**4*Kmv + 2*C_a_u**3*C_a_v*Kmfu + 2*C_a_u**3*C_a_v*Kmu + 2*C_a_u**3*Kmfu*Kmv + 2*C_a_u**3*Kmu*Kmv + C_a_u**2*C_a_v*Kmfu**2 + 4*C_a_u**2*C_a_v*Kmfu*Kmu + C_a_u**2*C_a_v*Kmu**2 + C_a_u**2*Kmfu**2*Kmv + 4*C_a_u**2*Kmfu*Kmu*Kmv + C_a_u**2*Kmu**2*Kmv + 2*C_a_u*C_a_v*Kmfu**2*Kmu + 2*C_a_u*C_a_v*Kmfu*Kmu**2 + 2*C_a_u*Kmfu**2*Kmu*Kmv + 2*C_a_u*Kmfu*Kmu**2*Kmv + C_a_v*Kmfu**2*Kmu**2 + Kmfu**2*Kmu**2*Kmv)
D2FV2 = J*(-C_a_u**2*C_a_v*Kmfu*Vmfv*r1 - 3*C_a_u**2*C_a_v*Kmfu*Vmfv*r2 - C_a_u**2*C_a_v*Kmfu*Vmfv*r3 - C_a_u**2*Kmfu*Kmv*Vmfv*r1 - 3*C_a_u**2*Kmfu*Kmv*Vmfv*r2 - C_a_u**2*Kmfu*Kmv*Vmfv*r3 + C_a_u**2*Kmu*Kmv*Vmu*r1*rq + 3*C_a_u**2*Kmu*Kmv*Vmu*r2*rq + C_a_u**2*Kmu*Kmv*Vmu*r3*rq - 2*C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r1 - 6*C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r2 - 2*C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r3 - 2*C_a_u*Kmfu*Kmu*Kmv*Vmfv*r1 - 6*C_a_u*Kmfu*Kmu*Kmv*Vmfv*r2 - 2*C_a_u*Kmfu*Kmu*Kmv*Vmfv*r3 + 2*C_a_u*Kmfu*Kmu*Kmv*Vmu*r1*rq + 6*C_a_u*Kmfu*Kmu*Kmv*Vmu*r2*rq + 2*C_a_u*Kmfu*Kmu*Kmv*Vmu*r3*rq - C_a_v*Kmfu*Kmu**2*Vmfv*r1 - 3*C_a_v*Kmfu*Kmu**2*Vmfv*r2 - C_a_v*Kmfu*Kmu**2*Vmfv*r3 + Kmfu**2*Kmu*Kmv*Vmu*r1*rq + 3*Kmfu**2*Kmu*Kmv*Vmu*r2*rq + Kmfu**2*Kmu*Kmv*Vmu*r3*rq - Kmfu*Kmu**2*Kmv*Vmfv*r1 - 3*Kmfu*Kmu**2*Kmv*Vmfv*r2 - Kmfu*Kmu**2*Kmv*Vmfv*r3)/(60*(C_a_u**4*C_a_v + C_a_u**4*Kmv + 2*C_a_u**3*C_a_v*Kmfu + 2*C_a_u**3*C_a_v*Kmu + 2*C_a_u**3*Kmfu*Kmv + 2*C_a_u**3*Kmu*Kmv + C_a_u**2*C_a_v*Kmfu**2 + 4*C_a_u**2*C_a_v*Kmfu*Kmu + C_a_u**2*C_a_v*Kmu**2 + C_a_u**2*Kmfu**2*Kmv + 4*C_a_u**2*Kmfu*Kmu*Kmv + C_a_u**2*Kmu**2*Kmv + 2*C_a_u*C_a_v*Kmfu**2*Kmu + 2*C_a_u*C_a_v*Kmfu*Kmu**2 + 2*C_a_u*Kmfu**2*Kmu*Kmv + 2*C_a_u*Kmfu*Kmu**2*Kmv + C_a_v*Kmfu**2*Kmu**2 + Kmfu**2*Kmu**2*Kmv))
D3FV2 = J*(-C_a_u**2*C_a_v*Kmfu*Vmfv*r1/120 - C_a_u**2*C_a_v*Kmfu*Vmfv*r2/60 - C_a_u**2*C_a_v*Kmfu*Vmfv*r3/60 - C_a_u**2*Kmfu*Kmv*Vmfv*r1/120 - C_a_u**2*Kmfu*Kmv*Vmfv*r2/60 - C_a_u**2*Kmfu*Kmv*Vmfv*r3/60 + C_a_u**2*Kmu*Kmv*Vmu*r1*rq/120 + C_a_u**2*Kmu*Kmv*Vmu*r2*rq/60 + C_a_u**2*Kmu*Kmv*Vmu*r3*rq/60 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r1/60 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r2/30 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r3/30 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r1/60 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r2/30 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r3/30 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r1*rq/60 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r2*rq/30 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r3*rq/30 - C_a_v*Kmfu*Kmu**2*Vmfv*r1/120 - C_a_v*Kmfu*Kmu**2*Vmfv*r2/60 - C_a_v*Kmfu*Kmu**2*Vmfv*r3/60 + Kmfu**2*Kmu*Kmv*Vmu*r1*rq/120 + Kmfu**2*Kmu*Kmv*Vmu*r2*rq/60 + Kmfu**2*Kmu*Kmv*Vmu*r3*rq/60 - Kmfu*Kmu**2*Kmv*Vmfv*r1/120 - Kmfu*Kmu**2*Kmv*Vmfv*r2/60 - Kmfu*Kmu**2*Kmv*Vmfv*r3/60)/(C_a_u**4*C_a_v + C_a_u**4*Kmv + 2*C_a_u**3*C_a_v*Kmfu + 2*C_a_u**3*C_a_v*Kmu + 2*C_a_u**3*Kmfu*Kmv + 2*C_a_u**3*Kmu*Kmv + C_a_u**2*C_a_v*Kmfu**2 + 4*C_a_u**2*C_a_v*Kmfu*Kmu + C_a_u**2*C_a_v*Kmu**2 + C_a_u**2*Kmfu**2*Kmv + 4*C_a_u**2*Kmfu*Kmu*Kmv + C_a_u**2*Kmu**2*Kmv + 2*C_a_u*C_a_v*Kmfu**2*Kmu + 2*C_a_u*C_a_v*Kmfu*Kmu**2 + 2*C_a_u*Kmfu**2*Kmu*Kmv + 2*C_a_u*Kmfu*Kmu**2*Kmv + C_a_v*Kmfu**2*Kmu**2 + Kmfu**2*Kmu**2*Kmv)
DM1FV2 = -C_a_u*Kmv*Vmu*rq*J*(2*r1 + 2*r2 + r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)
DM2FV2 =-C_a_u*Kmv*Vmu*rq*J*(r1 + 3*r2 + r3)/(60*C_a_u*C_a_v**2 + 120*C_a_u*C_a_v*Kmv + 60*C_a_u*Kmv**2 + 60*C_a_v**2*Kmu + 120*C_a_v*Kmu*Kmv + 60*Kmu*Kmv**2)
DM3FV2 = -C_a_u*Kmv*Vmu*rq*J*(r1 + 2*r2 + 2*r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)
# R_V: Delta_c1..c3..cM1..cM3 vectors for  n:

D1FV3 = J*(-C_a_u**2*C_a_v*Kmfu*Vmfv*r1/60 - C_a_u**2*C_a_v*Kmfu*Vmfv*r2/120 - C_a_u**2*C_a_v*Kmfu*Vmfv*r3/60 - C_a_u**2*Kmfu*Kmv*Vmfv*r1/60 - C_a_u**2*Kmfu*Kmv*Vmfv*r2/120 - C_a_u**2*Kmfu*Kmv*Vmfv*r3/60 + C_a_u**2*Kmu*Kmv*Vmu*r1*rq/60 + C_a_u**2*Kmu*Kmv*Vmu*r2*rq/120 + C_a_u**2*Kmu*Kmv*Vmu*r3*rq/60 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r1/30 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r2/60 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r3/30 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r1/30 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r2/60 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r3/30 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r1*rq/30 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r2*rq/60 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r3*rq/30 - C_a_v*Kmfu*Kmu**2*Vmfv*r1/60 - C_a_v*Kmfu*Kmu**2*Vmfv*r2/120 - C_a_v*Kmfu*Kmu**2*Vmfv*r3/60 + Kmfu**2*Kmu*Kmv*Vmu*r1*rq/60 + Kmfu**2*Kmu*Kmv*Vmu*r2*rq/120 + Kmfu**2*Kmu*Kmv*Vmu*r3*rq/60 - Kmfu*Kmu**2*Kmv*Vmfv*r1/60 - Kmfu*Kmu**2*Kmv*Vmfv*r2/120 - Kmfu*Kmu**2*Kmv*Vmfv*r3/60)/(C_a_u**4*C_a_v + C_a_u**4*Kmv + 2*C_a_u**3*C_a_v*Kmfu + 2*C_a_u**3*C_a_v*Kmu + 2*C_a_u**3*Kmfu*Kmv + 2*C_a_u**3*Kmu*Kmv + C_a_u**2*C_a_v*Kmfu**2 + 4*C_a_u**2*C_a_v*Kmfu*Kmu + C_a_u**2*C_a_v*Kmu**2 + C_a_u**2*Kmfu**2*Kmv + 4*C_a_u**2*Kmfu*Kmu*Kmv + C_a_u**2*Kmu**2*Kmv + 2*C_a_u*C_a_v*Kmfu**2*Kmu + 2*C_a_u*C_a_v*Kmfu*Kmu**2 + 2*C_a_u*Kmfu**2*Kmu*Kmv + 2*C_a_u*Kmfu*Kmu**2*Kmv + C_a_v*Kmfu**2*Kmu**2 + Kmfu**2*Kmu**2*Kmv)
D2FV3 = J*(-C_a_u**2*C_a_v*Kmfu*Vmfv*r1/120 - C_a_u**2*C_a_v*Kmfu*Vmfv*r2/60 - C_a_u**2*C_a_v*Kmfu*Vmfv*r3/60 - C_a_u**2*Kmfu*Kmv*Vmfv*r1/120 - C_a_u**2*Kmfu*Kmv*Vmfv*r2/60 - C_a_u**2*Kmfu*Kmv*Vmfv*r3/60 + C_a_u**2*Kmu*Kmv*Vmu*r1*rq/120 + C_a_u**2*Kmu*Kmv*Vmu*r2*rq/60 + C_a_u**2*Kmu*Kmv*Vmu*r3*rq/60 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r1/60 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r2/30 - C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r3/30 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r1/60 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r2/30 - C_a_u*Kmfu*Kmu*Kmv*Vmfv*r3/30 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r1*rq/60 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r2*rq/30 + C_a_u*Kmfu*Kmu*Kmv*Vmu*r3*rq/30 - C_a_v*Kmfu*Kmu**2*Vmfv*r1/120 - C_a_v*Kmfu*Kmu**2*Vmfv*r2/60 - C_a_v*Kmfu*Kmu**2*Vmfv*r3/60 + Kmfu**2*Kmu*Kmv*Vmu*r1*rq/120 + Kmfu**2*Kmu*Kmv*Vmu*r2*rq/60 + Kmfu**2*Kmu*Kmv*Vmu*r3*rq/60 - Kmfu*Kmu**2*Kmv*Vmfv*r1/120 - Kmfu*Kmu**2*Kmv*Vmfv*r2/60 - Kmfu*Kmu**2*Kmv*Vmfv*r3/60)/(C_a_u**4*C_a_v + C_a_u**4*Kmv + 2*C_a_u**3*C_a_v*Kmfu + 2*C_a_u**3*C_a_v*Kmu + 2*C_a_u**3*Kmfu*Kmv + 2*C_a_u**3*Kmu*Kmv + C_a_u**2*C_a_v*Kmfu**2 + 4*C_a_u**2*C_a_v*Kmfu*Kmu + C_a_u**2*C_a_v*Kmu**2 + C_a_u**2*Kmfu**2*Kmv + 4*C_a_u**2*Kmfu*Kmu*Kmv + C_a_u**2*Kmu**2*Kmv + 2*C_a_u*C_a_v*Kmfu**2*Kmu + 2*C_a_u*C_a_v*Kmfu*Kmu**2 + 2*C_a_u*Kmfu**2*Kmu*Kmv + 2*C_a_u*Kmfu*Kmu**2*Kmv + C_a_v*Kmfu**2*Kmu**2 + Kmfu**2*Kmu**2*Kmv)
D3FV3 = -J*(C_a_u**2*C_a_v*Kmfu*Vmfv*r1 + C_a_u**2*C_a_v*Kmfu*Vmfv*r2 + 3*C_a_u**2*C_a_v*Kmfu*Vmfv*r3 + C_a_u**2*Kmfu*Kmv*Vmfv*r1 + C_a_u**2*Kmfu*Kmv*Vmfv*r2 + 3*C_a_u**2*Kmfu*Kmv*Vmfv*r3 - C_a_u**2*Kmu*Kmv*Vmu*r1*rq - C_a_u**2*Kmu*Kmv*Vmu*r2*rq - 3*C_a_u**2*Kmu*Kmv*Vmu*r3*rq + 2*C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r1 + 2*C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r2 + 6*C_a_u*C_a_v*Kmfu*Kmu*Vmfv*r3 + 2*C_a_u*Kmfu*Kmu*Kmv*Vmfv*r1 + 2*C_a_u*Kmfu*Kmu*Kmv*Vmfv*r2 + 6*C_a_u*Kmfu*Kmu*Kmv*Vmfv*r3 - 2*C_a_u*Kmfu*Kmu*Kmv*Vmu*r1*rq - 2*C_a_u*Kmfu*Kmu*Kmv*Vmu*r2*rq - 6*C_a_u*Kmfu*Kmu*Kmv*Vmu*r3*rq + C_a_v*Kmfu*Kmu**2*Vmfv*r1 + C_a_v*Kmfu*Kmu**2*Vmfv*r2 + 3*C_a_v*Kmfu*Kmu**2*Vmfv*r3 - Kmfu**2*Kmu*Kmv*Vmu*r1*rq - Kmfu**2*Kmu*Kmv*Vmu*r2*rq - 3*Kmfu**2*Kmu*Kmv*Vmu*r3*rq + Kmfu*Kmu**2*Kmv*Vmfv*r1 + Kmfu*Kmu**2*Kmv*Vmfv*r2 + 3*Kmfu*Kmu**2*Kmv*Vmfv*r3)/(60*C_a_u**4*C_a_v + 60*C_a_u**4*Kmv + 120*C_a_u**3*C_a_v*Kmfu + 120*C_a_u**3*C_a_v*Kmu + 120*C_a_u**3*Kmfu*Kmv + 120*C_a_u**3*Kmu*Kmv + 60*C_a_u**2*C_a_v*Kmfu**2 + 240*C_a_u**2*C_a_v*Kmfu*Kmu + 60*C_a_u**2*C_a_v*Kmu**2 + 60*C_a_u**2*Kmfu**2*Kmv + 240*C_a_u**2*Kmfu*Kmu*Kmv + 60*C_a_u**2*Kmu**2*Kmv + 120*C_a_u*C_a_v*Kmfu**2*Kmu + 120*C_a_u*C_a_v*Kmfu*Kmu**2 + 120*C_a_u*Kmfu**2*Kmu*Kmv + 120*C_a_u*Kmfu*Kmu**2*Kmv + 60*C_a_v*Kmfu**2*Kmu**2 + 60*Kmfu**2*Kmu**2*Kmv)
DM1FV3 = -C_a_u*Kmv*Vmu*rq*J*(2*r1 + r2 + 2*r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)
DM2FV3 =-C_a_u*Kmv*Vmu*rq*J*(r1 + 2*r2 + 2*r3)/(120*C_a_u*C_a_v**2 + 240*C_a_u*C_a_v*Kmv + 120*C_a_u*Kmv**2 + 120*C_a_v**2*Kmu + 240*C_a_v*Kmu*Kmv + 120*Kmu*Kmv**2)
DM3FV3 = -C_a_u*Kmv*Vmu*rq*J*(r1 + r2 + 3*r3)/(60*C_a_u*C_a_v**2 + 120*C_a_u*C_a_v*Kmv + 60*C_a_u*Kmv**2 + 60*C_a_v**2*Kmu + 120*C_a_v*Kmu*Kmv + 60*Kmu*Kmv**2)

#constant = sp.simplify(sp.simplify(D3FV1)/( (2*r1 + r2 + 2*r3)/1))
#res1 =sp.simplify(sp.simplify(D1FV3)/constant)
#res2 =sp.simplify(sp.simplify(D2FV3)/constant)
#res3 =sp.simplify(sp.simplify(D3FV3)/constant)
#print( "SU1 ", res1)
#print( "SU2 ", res2)
#print( "SU3 ", res3)




substitutionU = [ (C_a_u, 101300*0.208/(8.32*293)), (C_a_v, 101300*0.0/(8.32*293)), (Kmv, 27.2438), (Kmu, 0.4103), (Vmu, 2.39e-4*math.exp( (80200/8.32) * (1/293 - 1/293) ) ) ]

substitutionV = [ (C_a_u, 101300*0.208/(8.32*293)), (C_a_v, 101300*0.0/(8.32*293)), (Kmv, 27.2438), (Kmu, 0.4103), (Vmu, 2.39e-4*math.exp( (80200/8.32) * (1/293 - 1/293) ) ), (Vmfv,1.61e-4*math.exp((56700/8.32)*(1/293 - 1/293))) ,(Kmfu, 0.1149), (rq, 0.97)]
########################################
###### evaluating and printing #########
########################################

'''

'''
print("F-vector: U part")
print("1-e-n: ", FU1.subs(substitutionU) )
print("e: ", FU2.subs(substitutionU) )
print("n: ", FU3.subs(substitutionU) )

print("F-vector: V part")
print("1-e-n: ", FV1.subs(substitutionV) )
print("e: ", FV2.subs(substitutionV) )
print("n: ", FV3.subs(substitutionV) )


print("Stiffness matrix: U part")


print("1-e-n: ")

print(D1FU1.subs(substitutionU) )
print(D2FU1.subs(substitutionU) )
print(D3FU1.subs(substitutionU) )
print(DM1FU1.subs(substitutionU) )
print(DM2FU1.subs(substitutionU) )
print(DM3FU1.subs(substitutionU) )

print("e: ", D1FV1.subs(substitutionU) )
print(D1FU2.subs(substitutionU) )
print(D2FU2.subs(substitutionU) )
print(D3FU2.subs(substitutionU) )
print(DM1FU2.subs(substitutionU) )
print(DM2FU2.subs(substitutionU) )
print(DM3FU2.subs(substitutionU) )

print("n: ", D1FV1.subs(substitutionU) )
print(D1FU3.subs(substitutionU) )
print(D2FU3.subs(substitutionU) )
print(D3FU3.subs(substitutionU) )
print(DM1FU3.subs(substitutionU) )
print(DM2FU3.subs(substitutionU) )
print(DM3FU3.subs(substitutionU) )

print("Stiffness matrix: V part")


print("1-e-n: ")

print(D1FV1.subs(substitutionV) )
print(D2FV1.subs(substitutionV) )
print(D3FV1.subs(substitutionV) )
print(DM1FV1.subs(substitutionV) )
print(DM2FV1.subs(substitutionV) )
print(DM3FV1.subs(substitutionV) )

print("e: ", D1FV1.subs(substitutionU) )
print(D1FV2.subs(substitutionV) )
print(D2FV2.subs(substitutionV) )
print(D3FV2.subs(substitutionV) )
print(DM1FV2.subs(substitutionV) )
print(DM2FV2.subs(substitutionV) )
print(DM3FV2.subs(substitutionV) )

print("n: ", D1FV1.subs(substitutionU) )
print(D1FV3.subs(substitutionV) )
print(D2FV3.subs(substitutionV) )
print(D3FV3.subs(substitutionV) )
print(DM1FV3.subs(substitutionV) )
print(DM2FV3.subs(substitutionV) )
print(DM3FV3.subs(substitutionV) )

'''
