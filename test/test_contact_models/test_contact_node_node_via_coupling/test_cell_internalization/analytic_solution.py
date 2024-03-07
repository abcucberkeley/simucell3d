import autograd.numpy as np   
from autograd import grad, jacobian
import scipy.optimize as optimize
from os import _exit


def H(r, a): 
    return 0.5 * (a**2 + r**2 + a * np.sqrt(a**2 + r**2))



def V(r, a):
    return (2./3.) * np.pi * (np.power((r**2) + (a**2), 1.5) + (a**3) + (3./2.) * (r**2) * a)






#-------------------------------------------------------------------------------------------------------------
def lagrangian(x):

    #Extract the variables
    r, a1, a2, ac, p1, p2 = x

    #Compute the volumes of the cells
    vc = V(r, ac)
    v1 = V(r, a1) + vc
    v2 = V(r, a2) - vc

    #Normalize the volumes by dividing the measured cell volume 2
    norm_v1 = v1 / measured_v2
    norm_v2 = v2 / measured_v2

    #The lagrangian with the 2 lagranian multipliers that enforce the volume constraints
    L = sigma * H(r, a1) + H(r, a2) +  alpha * H(r, ac) - p1 * (norm_v1 - beta**3) - p2 * (norm_v2  - 1.)
    return L
#-------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    #Automatically compute the gradients of the lagrangian and the jacobian of the lagrangian
    lagrangian_gradient = grad(lagrangian)
    lagrangian_jacobian = jacobian(lagrangian_gradient)

    #Compute the initial guess
    #psi_1 =  1.3843449168377018
    #psi_2 =  1.0623528186459092
    #psi_c =  1.7155917382548804
    #R1 =  1.378688313737518e-05
    #R2 =  1.550977639773187e-05
    #Rc =  1.3723307019310015e-05


    psi_1 =  1.4233297431807483
    psi_2 =  1.1485940336047804
    psi_c =  2.2605029847633573
    R1 =  1.364579089581736e-05
    R2 =  1.4796441802831265e-05
    Rc =  1.7479588910986176e-05

    r =  R1 * np.sin(psi_1)
    a1 = R1 * np.cos(psi_1)
    a2 = R2 * np.cos(psi_2)
    ac = Rc * np.cos(psi_c)

    #measured_v1 = 8.282e-15
    #measured_v2 = 8.781e-15

    measured_v1 = 8.374e-15
    measured_v2 = 8.766e-15

    #gamma_1 = 1.9e-3
    #gamma_2 = 1.0e-3 
    #gamma_c = 5.0e-4

    gamma_1 = 1.7e-3
    gamma_2 = 1.0e-3 
    gamma_c = 5.0e-4

    sigma = gamma_1 / gamma_2
    alpha = gamma_c / gamma_2
    beta = (measured_v1 / measured_v2)**(1./3.)

    #measured_pressure_1 = 275.5
    #measured_pressure_2 = 128.9

    measured_pressure_1 = 249.4
    measured_pressure_2 = 134.9

    characteristic_length = np.power(3. * measured_v2 / (4. * np.pi), 1./3.)
    normalized_p1 = measured_pressure_1 * characteristic_length
    normalized_p2 = measured_pressure_2 * characteristic_length

    #measured_cell_1_apical_area  =  7.05141e-10  
    #measured_cell_1_lateral_area =  9.7567e-10

    measured_cell_1_apical_area =  1.32669e-9
    measured_cell_1_lateral_area = 7.05141e-10

    x0 = np.array([r, a1, a2, ac, normalized_p1, normalized_p2])

    res = optimize.root(lagrangian_gradient, x0,  method='hybr', jac=lagrangian_jacobian, tol=1e-3, options={'maxfev': int(5e5)})
    root = res.x

    print("x0  ", ["{:.2e}".format(x) for x in x0],   "L=", lagrangian(x0), "grad_l", ["{:.2e}".format(x) for x in lagrangian_gradient(x0)] )
    print("root", ["{:.2e}".format(x) for x in root], "L=", lagrangian(root),"grad_l", ["{:.2e}".format(x) for x in lagrangian_gradient(root)] )

    #Compute the estimated apical area of cell 1
    theory_r, theory_a1, theory_a2, theory_ac, theory_p1, theory_p2 = root
    theory_apical_area_1  = 4. * np.pi * H(theory_r, theory_a1)
    theory_lateral_area_1 = 4. * np.pi * H(theory_r, theory_ac)

    print()
    print(measured_cell_1_apical_area, measured_cell_1_lateral_area, measured_cell_1_apical_area / (measured_cell_1_apical_area + measured_cell_1_lateral_area))
    print(theory_apical_area_1, theory_lateral_area_1, theory_apical_area_1 / (theory_apical_area_1+ theory_lateral_area_1) )





    #print(4. * np.pi * gamma_2 * theory_p1)



