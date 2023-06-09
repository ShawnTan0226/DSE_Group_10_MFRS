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
    def __init__(self,plane, tail,Steady_state, CmO_airfoil, h=5000 ,V=110): #Altitude in ft
        self.plane=plane
        self.tail=tail
        self.Steady_state=Steady_state

        self.MAC=self.plane.MAC
        self.b=self.plane.b_tot
        self.x_ac=self.plane.x_quarter
        self.x_cg=self.plane.x_cg
        self.dihedral=self.plane.dihedral

        self.x_tail=self.tail.x_tail
        self.z_tail=self.tail.z_tail
        self.S_tail=self.tail.S_tail
        self.Cl_alpha_v=self.tail.C_L_alpha_v
        self.sweep_half_v=self.tail.sweep_half_v
        self.Av=self.tail.A_v
        self.Sv=self.tail.S_v

        self.l_v=self.x_tail-self.plane.x_cg
        self.z_v=self.z_tail-self.plane.z_cg
        
        self.h=Steady_state.h
        self.V=Steady_state.V
        self.C_L=Steady_state.C_L
        self.C_D=Steady_state.C_D
        self.aoa=Steady_state.aoa
        self.T_c=Steady_state.T_c
        self.horisteady=Steady_state.horisteady


        self.CmO_root_list = np.array([])
        self.CmO_tip_list = np.array([])

        self.atmos()
        self.aerodynamic_properties()
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
            data[file]=df

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
    def define_C_f_wing(self, laminar_frac, part):
        C_f_laminar = 1.328 / np.sqrt(self.Re_list[part])
        C_f_turbulent = 0.455 / ((np.log10(self.Re_list[part])) ** 2.58 * (1 + 0.144 * self.M ** 2))
        C_f_total = laminar_frac * C_f_laminar + (1 - laminar_frac) * C_f_turbulent
        print('C_f_total', C_f_total)
        return C_f_total

    def define_C_D_part_wing(self, laminar_frac, part):
        # modern blended winglet IF= 1-1.01, so IF can be neglected
        C_f = self.define_C_f_wing(laminar_frac, part)
        FF = (1 + 0.6 / self.plane.max_thickness_location * self.plane.max_thickness + 100 * (self.plane.max_thickness) ** (4)) * (
                    1.34 * self.M ** 0.18 * (np.cos(self.plane.sweep[part])) ** 0.28)
        print('FF', FF)
        S_wet = 2 * self.plane.S_list[part] * 1.07
        C_D_part = FF * C_f * S_wet / self.plane.S
        return C_D_part
    
    def define_C_f_other(self, laminar_frac, length):
        Re = self.rho * self.V * length / self.nu
        C_f_laminar = 1.328 / np.sqrt(Re)
        C_f_turbulent = 0.455 / ((np.log10(Re)) ** 2.58 * (1 + 0.144 * self.M ** 2))
        C_f_total = laminar_frac * C_f_laminar + (1 - laminar_frac) * C_f_turbulent
        print('C_f_total', C_f_total)
        return C_f_total
    
    def define_C_D_part_nacelle(self, laminar_frac):
        cabinlength=17.40
        l=1.5/9.6*cabinlength
        d=0.4/9.6*cabinlength
        C_f = self.define_C_f_other(laminar_frac, l)
        f=l/d
        FF = 1+0.35/f
        S_wet = np.pi*d**2/4+np.pi*d*l
        IF=0.036*(self.plane.MAC*d/self.plane.S)*(0.2)**2
        C_D_part = FF * C_f * S_wet / self.plane.S+IF
        return C_D_part

    def define_C_D_0(self, laminar_frac):
        self.C_D_0 = 0
        for i in range(len(self.plane.taper)):
            self.C_D_0 += self.define_C_D_part_wing(laminar_frac, i)
        self.C_D_0 += self.define_C_D_part_nacelle(laminar_frac)
        self.C_D_0 += 0.1 * self.C_D_0
        return self.C_D_0

    def calc_Cmac(self):

        self.Cmac = ((self.plane.A_list * (np.cos(self.plane.sweep) ** 2)) /
                     (self.plane.A_list + 2 * np.cos(self.plane.sweep))) * (self.CmO_root_list + self.CmO_tip_list) / 2\
                    + self.dCmdEps[1:]*self.twist[1:]
        print("Cmac", self.Cmac)

    ### DYNAMIC STABILITY COEFFICIENTS ###

    def CLCDalphadot(self):
        self.C_L_alphadot=0
        self.C_D_alphadot=0

    def calc_C_Z_0(self):
        self.C_X_0=-self.C_l-self.T_c*(self.aoa+self.plane.ip)

    def calc_C_X_0(self):
        #Steady flight
        if self.horisteady:
            self.C_Z_0=0
        else:
            self.C_Z_0=-self.C_D+self.T_c
        


    #----------------Speed derivatives----------------
    def calc_C_X_u(self):
        self.C_X_u=-3*self.C_D*(1-self.C_D_T_c)

    def calc_C_Z_u(self):
        self.C_Z_u=-2*self.C_L+self.C_D*(-(self.alpha0+self.plane.ip)+3*self.C_D_T_c)

    def calc_C_m_u(self):
        self.C_m_u=-3*self.C_D*self.C_m_T_c

    #----------------angle of attack derivatives----------------
    def calc_C_X_alpha(self):
        self.C_X_alpha=-self.C_L*(1-2*self.C_L_alpha/(np.pi*self.plane.A*self.e))

    def calc_C_Z_alpha(self):
        self.C_Z_alpha=-self.C_L_alpha-self.C_D

    def calc_C_m_alpha(self):
        self.C_m_alpha=-self.C_L_alpha*(self.x_ac/self.MAC)

    #----------------pitch rate derivatives----------------
    def calc_C_X_q(self):
        self.C_X_q=0
    
    def calc_C_Z_q(self):
        self.C_Z_q=-self.C_L_alpha*(self.x_ac-self.x_cg)/(self.MAC)

    def calc_C_m_q(self):
        self.C_m_q=-self.C_L_alpha*(self.x_ac-self.x_cg)**2/(self.MAC**2) 

    #----------------Side slip derivatives----------------
    def calc_C_Y_beta(self):
        #wing
        self.C_Y_b_w = -0.00573 * np.rad2deg(self.dihedral)  # dihedral in deg

        #single vertical tail
        constant = 0.724 + 3.06*((self.Sv/self.plane.S)/(1+np.cos(self.plane.sweep_eq))) + 0.009 * self.plane.A #
        print("Considering the vertical stabiliser airfoil and placement what are the following values ? (Roskam 6 - p. 386)")
        AvfAv=float(input('Avf/Av (if you dont know assume 1'))
        Av_hfAvf = float(input('Avf/Av (if you dont know assume 1.2)'))
        Kvh = float(input('Kvh (if you dont know assume 1.2)'))
        self.A_v_eff = AvfAv * self.Av * (1+ Kvh*((Av_hfAvf)-1))
        self.C_L_alpha_v = 2*np.pi*self.A_v_eff / (2+((self.A_v_eff**2/(self.Cl_alpha_v/(2*np.pi))**2)*(1+(np.tan(self.sweep_half_v))**2)+4)**0.5) #omiting beta correction eq. 8.22

        if self.tail.tailnumber == 1:
            tn=1
        else:
            tn=2
        # Assuming bv/2ri>3.5-> kv=1
        self.C_Y_b_v = -tn*self.C_L_alpha_v * constant *self.Sv/self.plane.S ##Roskam 6 eq.10.28, kv=1 - Roskam 6 fig. 10.12
        # Simplify instead of eq.10.32 Roskam 6, omits effect of vertical stabiliser on each other

        #Final C_Y_beta
        self.C_Y_b =  self.C_Y_b_w + self.C_Y_b_v

    def calc_C_l_beta(self): #eq: 10.34 Roskam 6

        print("Considering: sweep at half chord {}, Aspect ratio {} and Taper ratio {}? (Roskam 6 - p. 393)".format(np.rad2deg(self.plane.sweep_eq_half),self.plane.A,self.plane.taper_eq))#Use of sweep equivalent
        clbetacLwf = float(input("What is Cl_beta/CL from fig. 10.20"))

        print("Considering: Aspect ratio {} and Taper ratio {}".format(self.plane.A,self.plane.taper_eq))
        clbetacLA = float(input("What is Cl_beta/CL from fig. 10.23"))

        print("Considering: sweep at half chord {}, Aspect ratio {} and Taper ratio {}? (Roskam 6 - p. 395)".format(np.rad2deg(self.plane.sweep_eq_half),self.plane.A,self.plane.taper_eq))
        clbetadihed = float(input("What is Cl_beta/dihedral from fig. 10.24"))

        print("Considering:Aspect ratio {} and Taper ratio {} (Roskam 6 - p. 396)".format(self.plane.A,self.plane.taper_eq))
        x = float(input("What is the dCL/(epsilon*tan(sweep)) from fig 10.26"))

        #Wing-fuselage
        self.C_l_beta_wf = 57.3*(self.C_L*(clbetacLwf*1*1 + clbetacLA) + self.dihedral*clbetadihed + x * self.plane.twist[-1] * np.tan(self.plane.sweep_eq))#Assumed no compressibility correction, fuselage correction =1 since no fuselage

        #Vertical Tail
        self.C_l_beta_v = self.C_Y_b_v*(self.z_v*np.cos(self.aoa) - self.l_v*np.sin(self.aoa))/self.plane.b[-1] #Roskam 6 eq. 10.34

    def calc_C_n_beta(self): #eq: 10.35 Roskam 6
        self.C_n_beta=-self.C_Y_b_v*(self.l_v*np.cos(self.aoa) + self.z_v*np.sin(self.aoa))/self.plane.b[-1]

    def calc_C_Y_betadot(self):
        self.C_Y_betadot = 0 #Usually neglected for high AR

    #----------------Roll rate derivatives----------------
    def calc_C_Y_p(self): 
        #Should this be neglected?
        self.C_Y_p=2*self.C_Y_b_v*(self.z_v*np.cos(self.aoa)-self.l_v*np.sin(self.aoa)-self.z_v)/self.b

    def calc_C_l_p(self):
        Beta_Comp=(1-self.M**2)**0.5
        Gamma_b=np.arctan(np.tan(self.plane.sweep_eq)/Beta_Comp)
        k=self.c_l_alpha*Beta_Comp/(2*np.pi) 
        BA_k=Beta_Comp*self.plane.A/k
        print("Gamma_b =",Gamma_b,"BA_k =",BA_k,"Taper =",self.plane.taper_eq,"What is BetaClp/k? (Roskam 6 - p. 418)")
        BetaClp_k=float(input('BetaClp/k ='))
        print("A =",self.plane.A,"Quarter chord sweep(Lambda_1/4) =",self.plane.sweep_eq,"What is Clp/CL^2? (Roskam 6 - p. 420)")
        Clp_CL2=float(input('Clp/CL^2 ='))
        Deltaclp_drag=Clp_CL2*self.C_L**2-0.125*self.C_D_0

        c_l_p_w=BetaClp_k*k/Beta_Comp*1+Deltaclp_drag
        c_l_p_v=2/self.plane.b_tot**2((self.z_v*np.cos(self.aoa)-self.l_v*np.sin(self.aoa))(self.z_v*np.cos(self.aoa)-self.l_v*np.sin(self.aoa)-self.z_v))*self.C_Y_b_v
        self.C_l_p=c_l_p_w+c_l_p_v
    
    def calc_C_n_p(self):
        print("A =",self.plane.A,"Taper ratio =",self.plane.taper_eq,"What is Cnp/eps? (Roskam 6 - p. 420)")
        Cnp_eps=float(input('Cnp/eps ='))
        Sweep=self.plane.sweep_eq
        A=self.plane.A
        xbar_cbar=0.25
        C_n_p_CL=-1/6*(A+6*(A+np.cos(Sweep))*(xbar_cbar*np.tan(Sweep)/A+np.tan(Sweep)**2/12))/(A+4*np.cos(Sweep))
        self.C_n_p_w=C_n_p_CL*self.C_L+Cnp_eps*self.plane.twist[-1]
        self.C_n_p_v=-2/self.b**2((self.l_v*np.cos(self.aoa)+self.z_v*np.sin(self.aoa))(self.z_v*np.cos(self.aoa)-self.l_v*np.sin(self.aoa)-self.z_v))*self.C_Y_b_v
        self.C_n_p=self.C_n_p_w+self.C_n_p_v

    #----------------Yaw rate derivatives----------------
    def calc_C_Y_r(self):
        # self.C_Y_r=-2*self.C_Y_beta_v*(self.x_ac_v-self.x_cg)/(self.MAC)
        #NOTE: incorporate x_ac_v and z_ac_b and z_cg in the plane object code
        self.lv = self.x_ac_v - self.x_cg
        self.zv = self.z_ac_v - self.z_cg

        # self.calc_C_Y_beta()

        self.C_Y_r = -2 * self.C_Y_beta_v * (self.lv * np.cos(self.aoa) + self.zv * np.sin(self.aoa)) / self.plane.b[-1]


    def calc_C_l_r(self):
        self.C_l_r = self.C_l_r_w + self.C_l_r_v

        self.B = (1 - (self.M**2) * np.cos(self.plane.sweep_eq**2))**0.5
        self.C_l_r_slope_low = float(input("From figure 10.41 (Page 430) in Roskam 6, provide slope value (OTHERWISE: assume 3.3 for A=6, Taper =0.26 and sweep=38deg"))

        self.C_l_r_slope = ((1 + (self.plane.A * (1 - self.B**2)) / (2 * self.B * (self.plane.A * self.B + 2 * np.cos(self.plane.sweep_eq)))) +
                           (((self.plane.A * self.B + 2*np.cos(self.plane.sweep_eq)) / (self.plane.A * self.B + 4*np.cos(self.plane.sweep_eq))) * ((np.tan(self.plane.sweep_eq)**2)/8)) * self.C_l_r_slope_low) / \
                           (1 + ((self.plane.A + 2*np.cos(self.plane.sweep_eq))/ (self.plane.A + 4*np.cos(self.plane.sweep_eq))) * ((np.tan(self.plane.sweep_eq)**2) / 8))

        self.C_l_r_gamma = 0.083 * (np.pi * self.plane.A * np.sin(self.plane.sweep_eq)) / (self.plane.A + 4*np.cos(self.plane.sweep_eq))

        self.C_l_r_epsilon = float(input("From figure 10.42 (page 430) in Roskam 6, provide slope value (OTHERWISE: assume 0.125 for A=6, Taper=0.26 and sweep=38"))

        self.C_l_r_dflaps = float(input("From figure 10.43 (page 431) in Roskam 6, provide slope value (OTHERWISE: assume 0 since dflaps is 0"))
        self.dflaps = 0 #Assume no flaps 
        self.aoa_dflaps = 0 #assume 0 because flaps is 0 but technically


        self.C_l_r_w = (self.C_L) * self.C_l_r_slope + self.C_l_r_gamma * self.dihedral + self.C_l_r_epsilon * self.twist[-1] + self.C_l_r_dflaps * self.aoa_dflaps * self.dflaps


        self.lv = self.x_ac_v - self.plane.x_cg
        self.zv = self.z_ac_v - self.z_cg
        self.C_l_r_v = -(2/(self.plane.b[-1]**2)) * (self.lv * np.cos(self.aoa) + self.zv * np.sin(self.aoa)) * \
                       (self.zv * np.cos(self.aoa) - self.lv * np.sin(self.aoa)) * self.C_Y_b_v

        self.C_l_r = self.C_l_r_w + self.C_l_r_v


    
    def calc_C_n_r(self):
        self.lv = self.x_ac_v - self.plane.x_cg
        self.zv = self.z_ac_v - self.z_cg

        self.C_n_r_cl = float(input("From figure 10.44 (page 433) in Roskam 6, provide slope value"))
        self.C_n_r_cd0 = float(input("From figure 10.45 (page 434) in Roskam 6, provide slope value"))

        self.C_n_r_w = self.C_n_r_cl * self.C_L**2 + self.C_n_r_cd0 * self.C_D_0

        self.C_n_r_v = (2/(self.plane.b[-1]**2)) * ((self.lv * np.cos(self.aoa) + self.zv * np.sin(self.aoa))**2) * self.C_Y_b_v

        self.C_n_r = self.C_n_r_w + self.C_n_r_v

    #def Cmdot(self): #Make CLalpha_h equal to zero if there is no tail, tail is a straight conventional wing at wingtips
    #    self.eta_h =
    #    self.x_ac_h = self.plane.offset[-1]+ self.plane.c
    #    self.Cmdot = -2*self.CLalpha_h * self.eta_h *(self.plane.offset[-1]+ self.C


    

            
    


    def help(self,option):
        print("Help for the coefficients, options:\n -Vertical Tail \n -Wing \n -Cmac")
        if option == 'Vertical Tail':
            print("The vertical tail sizing takes the follwing input\n -")
        if option == 'Wing':
            print("The wing sizing takes the follwing input\n -dihedral \n -")


'''Inputs'''
#Aerodynamic coefficients




# dCmdEps = [0,0.4,0.4]
# twist = [0,0.1,0.1]
# CmO_airfoil = [[0.1,0.05],[0.2,0.1]] #Get from xfoil
# plane=Plane(9.72,[0.3,0.267],[38,38],[8.82,22])
# coefficients = AerodynamicProperties(plane, dCmdEps, twist, CmO_airfoil)
# print(coefficients.Cmac())
