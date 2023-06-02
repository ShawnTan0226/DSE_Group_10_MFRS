#Libraries
import numpy as np
import matplotlib.pyplot as plt
import math



#Mean aerodynamic chord
def MAC_part(cr, ct, sweep, b):
    #MAC
    y = 2/3 * (cr + ct - (cr*ct)/(cr+ct))
    #offset from that parts root on quarter chord
    off_x = (b/2/(ct-cr))*(y-ct)
    off_y = cr*0.25-ct*0.25 + np.tan(sweep)*off_x
    return y, off_x, off_y

def listgenerator(cri,cti,cro,cto,sweepi,sweepo,bi,bo):
    MACi, off_x_i, off_y_i = MAC_part(cri,cti,sweepi,bi)
    MACo, off_x_o, off_y_o = MAC_part(cro, cto, sweepo, bo)

    MAC_list = np.array([MACi,MACo])
    S_list = np.array(([(cri+cti)/2*bi,(cro+cto)/2*bo]))
    x_list = np.array(([off_x_i,bi+off_x_o]))
    return MAC_list,S_list,x_list

def MAC_aircraft(cri,cti,cro,cto,sweepi,sweepo,bi,bo):#Make sure to use numpy array
    MAC_list,S_list,x_list=listgenerator(cri, cti, cro, cto, sweepi, sweepo, bi, bo)
    MAC = np.sum((MAC_list*S_list))/np.sum(S_list)
    x_quarter = np.sum(x_list)/np.sum(S_list)
    return MAC, x_quarter

'''INPUT'''
cri=0
cti=0
cro=0
cto=0
sweepi=0
sweepo=0
bi=0
bo=0
'''CALCULATION'''
print(MAC_aircraft(cri,cti,cro,cto,sweepi,sweepo,bi,bo))