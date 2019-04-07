
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
#from scipy.optimize import minimize

q = 1.602e-19 #C
u = 1.661e-27 # kg

Vl = 1715 # V
L1 = 0.0145 #m
L2 = 0.0395 #m
L3 = 0.0079 #m
Ld = 0.9102 #m


#
#Vp0, Vmax,tpq,quad1,quad2= 5.69840909e+02, 1.02752244e+03, 2.75728232e-06, -9.52807888e+05, 6.11177932e+06

#tp = 2.85e-6


#Va,tp,q1,q2,q3= 8.91364804e+00, 2.40966217e-06, 6.0599639765e+04, 6.0e+11, 4.807091e+18

#Vex,Va,tp,q1,q2,q3=600,  3.82517787e+03,  2.8552042e-06,  8.1372595e+02,  9.911192e+08, 1.33964169e+14
#Va,tp,q1,q2,q3= 8.91364804e+00, 2.40966217e-06, 6.0599639765e+04, 6.0e+11, 4.807091e+18

#
#Vex,Va,tp =500,  4500,  2.70e-06,



Va,tp,q1,q2,q3=4.55111026e+03, 2.59781828e-06, 5.91763922e+02,4.72579070e+08,3.49607999e+13
#4.545012980e+03, 2.71408663e-06, 6.451283408e+02, 3.17067678e+08, 1.30376710e+14
Vex=5000 - Va


def quad_equations(m, v0):  
    # Repeller-Extractor region (vd1, td1 = extractor)
    v1 = np.sqrt(v0*v0 + 2*q*Vex/m)
    t1 = ((m*v0*L1)/(q*Vex))*(-1+np.sqrt(1+(2*q*Vex)/(m*v0*v0)))

    # Extractor-Control region (vd2, td2 = control, dzd = pulse position, vpd = velocity at pulse position)
    Lx = v1*(tp-t1)+((q*Va)/((2*m*L2))*(tp-t1)**2)

    vpv = np.sqrt((2/m)*(m*v0*v0/2 + q*Vex + q*Va*Lx/L2))
#    print(vpv)
    root_guess = 0.00003

    func = lambda root: (Lx-L2)+(vpv*root)+(q/(m*L2))*((((Va+q1)*(root)*(root)/2)+((q2)*(root)*(root)*(root))/6)+((q3)*(root)*(root)*(root)*(root))/12)
    root = fsolve(func, root_guess)
#    print(root)
    t2 = root + tp
    
    if root >= 0:

        v2 = vpv + (q/(m*L2))*(((Va+q1)*root)+((q2*(root)*(root))/2)+((q3*(root)*(root)*(root))/3))
    else: 
        
        v2 = np.sqrt(v0*v0 + (2*q*(Vex +Va))/m)

    # Einzel region
    
#    print(v2)
    t4 = t2 + 2*(((v2*m*L3)/(q*Vl)))*(1+np.sqrt(1-((2*q*Vl)/(m*v2*v2))))
#    print('here',tv4)
   
    
    
    # Ground-Detector
    tf = t4 + Ld/v2        

    return tf
quad_equations = np.vectorize(quad_equations)



XQ = np.arange(325,610,10)
YQ = np.arange(10,1010,50)
XQ, YQ = np.meshgrid(XQ, YQ)
ZQ = quad_equations(XQ*u,YQ)-quad_equations(XQ*u, 500) # value from quad equation subtract average value 
ttd = np.sum(np.absolute(ZQ)) #total variation in flight time 
print("total time =",ttd)



fig,ax5= plt.subplots()
mapp = ax5.imshow(ZQ, origin='lower', cmap='Spectral', vmin=-1e-9, vmax=1e-9)
#ax5.set_xticks(np.arange(0,101,20)) 
#ax5.set_xticklabels(np.arange(300,601,60))
#ax5.set_yticks(np.arange(0,101,5)) 
#ax5.set_yticklabels(np.arange(0,1010,250))
cbaxes = fig.add_axes([0.93, 0.125, 0.025, 0.75])
fig.colorbar(mapp,cax=cbaxes, ticks=[-1e-9, 0, 1e-9],orientation='vertical')
