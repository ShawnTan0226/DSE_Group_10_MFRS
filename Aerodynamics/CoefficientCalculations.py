import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import interp1d
from Plane import Plane

# Twist, dCmdEps, CmO_airfoil should be structured as follows
# Twist = [0, twist of 1st part, twist of 2nd part, ...]
# dCmdEps = [0, dCmdEps of 1st part, dCmdEps of 2nd part, ...] Source: Roskam 6
# CmO_airfoil = [[CmO_root, CmO_tip] for 1st part, [[CmO_root, CmO_tip] for 2nd part, ...] Source: Roskam 6
#The horizontal stabiliser are at the wing tips and have a sweep of zeroµ

#TODO:check all the inputs units
class AerodynamicProperties:
    def __init__(self,plane, dCmdEps, twist, CmO_airfoil, h=5000 ,V=110): #Altitude in ft
        self.h=h
        self.V=V
        self.plane=plane

        self.dCmdEps = np.array(dCmdEps)
        self.twist = np.array(twist)

        self.CmO_root_list = np.array([])
        self.CmO_tip_list = np.array([])
        for i in range(len(CmO_airfoil)):
            self.CmO_root_list = np.concatenate((self.CmO_root_list,[CmO_airfoil[i][0]]))
            self.CmO_tip_list = np.concatenate((self.CmO_tip_list,[CmO_airfoil[i][1]]))

        directory_path = "./Xflr dat"
        file_list = []
        data={}

        for filename in os.listdir(directory_path):
            file_list.append(filename)

        for file in file_list:
            filename=directory_path+'/'+file
            df = pd.read_csv(filename,header=None)
            df=np.array(df)
            data[file[:4]]=df

        self.data=data[self.plane.planename]
        self.calc_C_L_alpha()

    def calc_C_L_alpha(self):
        start=np.where(self.data[:,0]==-5)[0][0]
        end=np.where(self.data[:,0]==5)[0][0]
        coeff = np.polyfit(self.data[start:end,0],self.data[start:end,2],1)
        self.C_L_alpha=coeff[0]
        return self.C_L_alpha

    ### ISA CALCULATIONS ###
    def atmos(self):
        self.T = 288.15 - 0.0065 * self.h
        self.rho = 1.225 * (self.T / 288.15) ** (-1 + 9.81 / (287.058 * 0.0065))
        self.nu = 0.0000169 #Must be changed manually

    def aerodynamic_properties(self):
        self.atmos()
        self.Re = self.rho * self.V * self.plane.MAC / self.nu
        self.Re_list = self.rho * self.V * self.plane.MAC_list / self.nu
        self.a = np.sqrt(1.4 * 287.058 * self.T)
        self.M = self.V / self.a

    ### AERODYNAMIC COEFFICIENT ###
    def define_C_f(self, laminar_frac, part):
        C_f_laminar = 1.328 / np.sqrt(self.Re_list[part])
        C_f_turbulent = 0.455 / ((np.log10(self.Re_list[part])) ** 2.58 * (1 + 0.144 * self.M ** 2))
        C_f_total = laminar_frac * C_f_laminar + (1 - laminar_frac) * C_f_turbulent
        print('C_f_total', C_f_total)
        return C_f_total

    def define_C_D_part_wing(self, laminar_frac, part):
        # modern blended winglet IF= 1-1.01, so IF can be neglected
        C_f = self.define_C_f(laminar_frac, part)
        FF = (1 + 0.6 / self.plane.max_thickness_location * self.plane.max_thickness + 100 * (self.plane.max_thickness) ** (4)) * (
                    1.34 * self.M ** 0.18 * (np.cos(self.plane.sweep[part])) ** 0.28)
        print('FF', FF)
        S_wet = 2 * self.plane.S_list[part] * 1.07
        C_D_part = FF * C_f * S_wet / self.plane.S
        return C_D_part

    def define_C_D_0(self, laminar_frac):
        self.C_D_0 = 0
        for i in range(len(self.plane.taper)):
            self.C_D_0 += self.define_C_D_part_wing(laminar_frac, i)
        self.C_D_0 += 0.1 * self.C_D_0
        return self.C_D_0

    def calc_Cmac(self):

        self.Cmac = ((self.plane.A_list * (np.cos(self.plane.sweep) ** 2)) /
                     (self.plane.A_list + 2 * np.cos(self.plane.sweep))) * (self.CmO_root_list + self.CmO_tip_list) / 2\
                    + self.dCmdEps[1:]*self.twist[1:]
        print("Cmac", self.Cmac)

    ### DYNAMIC STABILITY COEFFICIENTS ###

    def CLCDalphadot(self):
        self.CLdot=0
        self.CDdot=0

    def calc_C_Z_0(self):
        self.C_Z_0=0

    def calc_C_X_0(self):
        self.C_X_0=-self.C_l-self.T_c*(self.alpha0+self.ip)


    #Speed derivatives
    def calc_C_X_u(self):
        self.C_X_u=-3*self.C_D*(1-self.C_D_T_c)

    def calc_C_Z_u(self,ip):
        self.C_Z_u=-2*self.C_L+self.C_D*(-(self.alpha0+ip)+3*self.C_D_T_c)

    def calc_C_m_u(self):
        self.C_m_u=-3*self.C_D*self.C_m_T_c

    #angle of attack derivatives
    def calc_C_X_alpha(self):
        self.C_X_alpha=-self.C_L*(1-2*self.C_L_alpha/(self.pi*self.plane.A*self.e))

    def calc_C_Z_alpha(self):
        self.C_Z_alpha=-self.C_L_alpha-self.C_D

    def calc_C_m_alpha(self):
        self.C_m_alpha=-self.C_L_alpha*(self.x_ac/self.MAC)

    #pitch rate derivatives
    def calc_C_X_q(self):
        self.C_X_q=0
    
    def calc_C_Z_q(self):
        self.C_Z_q=-self.C_L_alpha*(self.x_ac-self.x_cg)/(self.MAC)

    def calc_C_m_q(self):
        self.C_m_q=-self.C_L_alpha*(self.x_ac-self.x_cg)**2/(self.MAC**2)

    #Roll rate derivatives
    def calc_C_Y_p(self):
        #Should this be neglected?
        self.C_Y_p=(self.z_vcos(self.aoa)-self.l_v*np.sin(self.aoa)-self.z_v)/self.plane.b_tot*self.C_Y_b_v

    def calc_C_l_p(self):
        y=np.arange(0,self.plane.b_tot/2,self.plane.b_tot/2000)
        Delta_alpha=y
        Delta_C_l=self.normalisedC_L*Delta_alpha
        self.C_l_p=-2*np.trapz(Delta_C_l*y,y)
    
    def calc_C_n_p(self):
        self.C_n_p_w=0
        self.C_n_p_v=(self.l_v*np.cos(self.aoa) +self.z_v*np.sin(self.aoa))(self.z_vcos(self.aoa)-self.l_v*np.sin(self.aoa)-self.z_v)/self.plane.b_tot*self.C_Y_b_v
        self.C_n_p=self.C_n_p_w+self.C_n_p_v

    #Yaw rate derivatives
    def calc_C_Y_r(self):
        self.C_Y_r=-2*self.C_Y_beta_v*(self.x_ac_v-self.x_cg)/(self.MAC)

    def calc_C_l_r(self):
        self.C_l_r=0
    
    def calc_C_n_r(self):
        self.C_n_r=0

    #def Cmdot(self): #Make CLalpha_h equal to zero if there is no tail, tail is a straight conventional wing at wingtips
    #    self.eta_h =
    #    self.x_ac_h = self.plane.offset[-1]+ self.plane.c
    #    self.Cmdot = -2*self.CLalpha_h * self.eta_h *(self.plane.offset[-1]+ self.C


    def calc_C_Y_beta(self):
        #wing
        self.C_Y_b_w = -0.00573 * self.dihedral  # dihedral in deg

        print("Are you training ")
        cz=float(input())

        #single vertical tail
        constant = 0.724 + 3.06*((self.Sv/self.plane.S)/(1+np.cos(self.plane.sweep_eq))) + 0.009 * self.plane.A #
        print("Considering the vertical stabiliser airfoil and placement what are the following values ? (Roskam 6 - p. 386)")
        AvfAv=float(input('Avf/Av (if you dont know assume 1'))
        Av_hfAvf = float(input('Avf/Av (if you dont know assume 1.2)'))
        Kvh = float(input('Kvh (if you dont know assume 1.2)'))
        self.A_v_eff = AvfAv * self.Av * (1+ Kvh*((Av_hfAvf)-1))
        self.C_L_alpha_v = 2*np.pi*self.A_v_eff / (2+((self.A_v_eff**2/(self.Cl_alpha_v/(2*np.pi))**2)*(1+(np.tan(self.sweep_half_v))**2)+4)**0.5) #omiting beta correction eq. 8.22

        # Assuming bv/2ri>3.5-> kv=1
        self.C_Y_b_v_single = -self.CL_alpha_v * constant *self.Sv/self.plane.S ##Roskam 6 eq.10.28, kv=1 - Roskam 6 fig. 10.12

        #Twin vertical tail
        self.C_Y_b_v_twin = 2*self.C_Y_b_v_single #Simplify instead of eq.10.32 Roskam 6, omits effect of vertical stabiliser on each other

        #Final C_Y_beta
        self.C_Y_b_single =  self.C_Y_b_w + self.C_Y_b_v_single
        self.C_Y_b_twin = self.C_Y_b_w + self.C_Y_b_v_twin

    def calc_C_l_beta(self):
        #Wing-fuselage



    def help(self,option):
        print("Help for the coefficients, options:\n -Vertical Tail \n -Wing \n -Cmac")
        if option == 'Vertical Tail':
            print("The vertical tail sizing takes the follwing input\n -")
        if option == 'Wing':
            print("The wing sizing takes the follwing input\n -dihedral \n -")


'''Inputs'''
#Aerodynamic coefficients




dCmdEps = [0,0.4,0.4]
twist = [0,0.1,0.1]
CmO_airfoil = [[0.1,0.05],[0.2,0.1]] #Get from xfoil
plane=Plane(9.72,[0.3,0.267],[38,38],[8.82,22])
coefficients = AerodynamicProperties(plane, dCmdEps, twist, CmO_airfoil)
print(coefficients.Cmac())
