#Libraries
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.interpolate as sp

from Plane import Plane



class control:
    def __init__(self,plane,tail,cf_c,delta,CL_max,Cmac):
        self.plane=plane
        self.tail=tail
        self.cf_c=cf_c
        self.delta=np.deg2rad(delta)
        self.CL_max=CL_max
        self.Cmac=Cmac


        self.calc_x_cg_MAC()
        self.calc_dCm_req()
        print('dcm----',self.dCm_req)
        self.calc_eta_i()
    
    def moment_diff(self,eta_i):
        return self.dCm_req+self.calc_Delta_Cm(eta_i,0.95)

    def calc_eta_i(self):
        x=self.newtonRaphson(self.moment_diff,0.5,0.0001,100,0.0001,0.5)[2]
        self.eta_i=x
        

    def calc_dCm_req(self):
        self.dCm_req=self.CL_max*(self.plane.x_cg-self.plane.x_quarter)/self.plane.MAC+np.sum(self.Cmac)

    def calc_delta_Cm_outer(self):
        dcm=self.calc_delta_CL(self.delta,self.plane.b[1]/self.plane.b_tot+0.05,0.95)
        return dcm
        

    def calc_Delta_Cm(self,eta_i,eta_o):
        self.Delta_CL_ref=self.calc_delta_CL(self.delta,0,1)
        self.dcm_dCMref=0
        Delta_Cm=(self.x_cg_MAC-0.25)*self.calc_delta_CL(self.delta,eta_i,eta_o)+self.K_Lambda(eta_i,eta_o)*self.plane.A/1.5*self.Delta_CL_ref*np.tan(self.plane.sweep_eq)+self.K_Lambda(eta_i,eta_o)*self.Delta_CL_ref*self.dcm_dCMref
        return Delta_Cm




    def calc_x_cg_MAC(self):
        self.x_cg_MAC=(self.plane.x_cg-self.plane.x_quarter-self.plane.MAC*0.25)/self.plane.MAC

    def calc_delta_cm(self,delta):
        self.delta_cm=delta*self.f_cl_delta(self.cf_c)*self.f_k(delta)

    def K_Lambda(self,eta_i,eta_o):
        K_Lambda=[0,0.48,0.74,0.9,1,1.07,1.13,1.16,1.19,1.22,1.22]
        eta=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
        f=sp.interp1d(eta,K_Lambda)
        return f(eta_o)-f(eta_i)

    def K_b(self,eta_i,eta_o):
        K_b=[0,0.16,0.3,0.45,0.55,0.67,0.77,0.86,0.93,0.97,1]
        eta=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]

        f=sp.interp1d(eta,K_b)
        return f(eta_o)-f(eta_i)
    
    def calc_delta_cl(self,delta):

        self.delta_cl=delta*self.f_cl_delta(self.cf_c)*self.f_k(delta)

    def calc_delta_CL(self,delta,eta_i,eta_o):
        self.calc_delta_cl(delta)
        delta_CL=self.delta_cl*self.K_b(eta_i,eta_o)*self.f_ad_CL_ad_cl(self.cf_c)
        return delta_CL

    def f_ad_CL_ad_cl(self,cfc):
        ad_cl_list=[0,0.4,0.55,0.68,0.77,0.84,0.89,0.96,1]
        cfc_list=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
        f_adcl=sp.interp1d(cfc_list,ad_cl_list)
        ad_cl=f_adcl(cfc)
        print('Aspect ratio: ',self.plane.A,'ad_cl: ',ad_cl)
        self.ad_CL_ad_cl=1#float(input('Enter ad_CL/ad_cl(Roskam 6, pg 261, fig 8.53): '))
        return self.ad_CL_ad_cl

    def f_cl_delta(self,cfc):
        cl_delta=[1.75,2.5,3.25,3.75,4.25,4.6,5,5.3,5.7,6]
        cf_c=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]
        f=sp.interp1d(cf_c,cl_delta)
        return f(cfc)

    def f_k(self,delta):
        k=[1,1,1,0.9,0.7,0.6,0.55,0.52,0.48,0.46,0.45,0.43,0.42]
        deflection=[0,5,10,15,20,25,30,35,40,45,50,55,60]
        f_k=sp.interp1d(deflection,k)
        return f_k(delta)
    
    def gradient(self, f, x, step):
        return (f(x + step) - f(x - step)) / (step * 2)

    def newtonRaphson(self, f, x0, e, N, h, relax):
        #print('\n\n*** NEWTON RAPHSON METHOD IMPLEMENTATION ***')
        i = 0
        step = 1
        flag = 1
        condition = True
        while condition:
            # if g(f,x0,h) == 0.0:
            #     print('Divide by zero error!')
            #     break
            #print('x0---', x0)
            #print('value---', f(x0))
            #print('grad---', self.gradient(f, x0, h))
            x1 = x0 * relax + (x0 - f(x0) / (self.gradient(f, x0, h))) * (1 - relax)
            # print('Iteration-%d, x1 = %0.6f and f(x1) = %0.6f' % (step, x1, f(x1)))
            x0 = x1
            step = step + 1
            newvalue = f(x1)
            #print(newvalue)
            # if g(f,buildingno,x0,h)<0:
            #     x1=x1/relax

            if abs(np.max(newvalue)) < e:
                condition = False
            if step > N:
                print('\nNot Convergent.')
                flag = 2
                condition = False
            i += 1
            print('x1---', x1)

        if flag == 1:
            #print('\nRequired root is: %0.8f', x1)
            return x0, i, x1
        else:
            #print('\nNot Convergent.')
            return 1000, i, "No solution found"

        