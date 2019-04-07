
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

#Instrument variables 
q = 1.602e-19 #C
u = 1.661e-27 # kg
Vex=5000 - Va
Va=4.55111026e+03
tp=2.59781828e-06
q1=5.91763922e+02
q2=4.72579070e+08
q3=3.49607999e+13
Vl = 1715 # V
L1 = 0.0145 #m
L2 = 0.0395 #m
L3 = 0.0079 #m
Ld = 0.9102 #m

#Ion trajectories calculation
def quad_equations(m, v0):  
    # Repeller-Extractor region (vd1, td1 = extractor)
    v1 = np.sqrt(v0*v0 + 2*q*Vex/m)
    t1 = ((m*v0*L1)/(q*Vex))*(-1+np.sqrt(1+(2*q*Vex)/(m*v0*v0)))

    # Extractor-Control region (vd2, td2 = control, dzd = pulse position, vpd = velocity at pulse position)
    Lx = v1*(tp-t1)+((q*Va)/((2*m*L2))*(tp-t1)**2)
    vpv = np.sqrt((2/m)*(m*v0*v0/2 + q*Vex + q*Va*Lx/L2))
    root_guess = 0.00003
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
