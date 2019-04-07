
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import minimize

#Instrument variables 
q = 1.602e-19 #C
u = 1.661e-27 # kg
Vl = 1715 # V
L1 = 0.0145 #m
L2 = 0.0395 #m
L3 = 0.0079 #m
Ld = 0.9102 #m

#Ion trajectories 
def p(x):
    Va=x[0] 
    Vex=5000 - Va
    tp=x[1]
    q1=x[2]
    q2=x[3]
    q3=x[4]
    print(x)
    a=[]
    for m in range(325,610,10):           
        m = m*u
        for v0 in range(10,1010,50):    

    # Repeller-Extractor region (vd1, td1 = extractor)
            v1 = np.sqrt(v0*v0 + 2*q*Vex/m)
            t1 = ((m*v0*L1)/(q*Vex))*(-1+np.sqrt(1+(2*q*Vex)/(m*v0*v0)))
        
     # Extractor-Control region (vd2, td2 = control, dzd = pulse position, vpd = velocity at pulse position)
            Lx = v1*(tp-t1)+((q*Va)/((2*m*L2))*(tp-t1)**2)
            vpv = np.sqrt((2/m)*(m*v0*v0/2 + q*Vex + q*Va*Lx/L2))
            root_guess = 3
        
            func = lambda root: (Lx-L2)+(vpv*root)+(q/(m*L2))*((((Va+q1)*(root)*(root)/2)+((q2)*(root)*(root)*(root))/6)+((q3)*(root)*(root)*(root)*(root))/12)
            root = fsolve(func, root_guess)
            t2 = root + tp
            if root >= 0:
                v2 = vpv + (q/(m*L2))*(((Va+q1)*root)+((q2*(root)*(root))/2)+((q3*(root)*(root)*(root))/3))
            else: 
                v2 = np.sqrt(v0*v0 + (2*q*(Vex +Va))/m)
        
            # Einzel region
            t4 = t2 + 2*(((v2*m*L3)/(q*Vl)))*(1+np.sqrt(1-((2*q*Vl)/(m*v2*v2))))
  
            # Ground-Detector
            t5 = t4 + Ld/v2        
            tdiff = t5 - quad_equations(m,500)
            a = np.append(a,tdiff)
    a= np.absolute(a)
    tofsum= sum(a)       
    print('tof sum=', tofsum)
    return sum(a)

#Trajectory calculation 2 - used for time difference calculation 
def quad_equations(m, v0): 
    
    # Repeller-Extractor region (vd1, td1 = extractor)
    v1 = np.sqrt(v0*v0 + 2*q*Vex/m)
    t1 = ((m*v0*L1)/(q*Vex))*(-1+np.sqrt(1+(2*q*Vex)/(m*v0*v0)))

    # Extractor-Control region (vd2, td2 = control, dzd = pulse position, vpd = velocity at pulse position)
    Lx = v1*(tp-t1)+((q*Va)/((2*m*L2))*(tp-t1)**2)
    vpv = np.sqrt((2/m)*(m*v0*v0/2 + q*Vex + q*Va*Lx/L2))
    root_guess = 3
    func = lambda root: (Lx-L2)+(vpv*root)+(q/(m*L2))*((((Va+q1)*(root)*(root)/2)+((q2)*(root)*(root)*(root))/6)+((q3)*(root)*(root)*(root)*(root))/12)
    root = fsolve(func, root_guess)
    t2 = root + tp
    if root >= 0:
        v2 = vpv + (q/(m*L2))*(((Va+q1)*root)+((q2*(root)*(root))/2)+((q3*(root)*(root)*(root))/3))
    else: 
        v2 = np.sqrt(v0*v0 + (2*q*(Vex +Va))/m)

    # Einzel region
    t4 = t2 + 2*(((v2*m*L3)/(q*Vl)))*(1+np.sqrt(1-((2*q*Vl)/(m*v2*v2))))
    
    # Ground-Detector
    tf = t4 + Ld/v2        

    return tf     
    
#Optimisation variables 
Va,tp,q1,q2,q3=4.55000020e+03, 2.59458086e-06, 5.81519474e+02, 4.95802506e+08, 2.01288344e+13
Vex=5000 - Va

ini = np.array([Va,tp,q1,q2,q3])
res = minimize(p, ini, method='Nelder-Mead', options={'xtol': 1e-15, 'disp': True,'maxiter': 10000})
print('output=',res.x)
