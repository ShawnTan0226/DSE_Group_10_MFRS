import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import interp1d

# Twist, dCmdEps, CmO_airfoil should be structured as follows
# Twist = [0, twist of 1st part, twist of 2nd part, ...]
# dCmdEps = [0, dCmdEps of 1st part, dCmdEps of 2nd part, ...] Source: Roskam 6
# CmO_airfoil = [[CmO_root, CmO_tip] for 1st part, [[CmO_root, CmO_tip] for 2nd part, ...] Source: Roskam 6
class AerodynamicProperties:
    def __init__(self,plane, dCmdEps, twist, CmO_airfoil, h=5000 ,V=110): #Altitude in ft
        self.h=h
        self.V=V
        self.plane=plane

        self.dCmdEps = dCmdEps
        self.twist = twist

        self.CmO_root_list = np.array([])
        self.CmO_tip_list = np.array([])
        for i in range(len(CmO_airfoil)):
            self.CmO_root_list = np.concatenate((self.CmO_root_list,[CmO_airfoil[i][0]]))
            self.CmO_tip_list = np.concatenate((self.CmO_tip_list,[CmO_airfoil[i][1]]))

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

    def Cmac(self):
        self.Cmac = ((self.plane.A_list * (np.cos(self.plane.sweep) ** 2)) /
                     (self.plane.A_list + 2 * np.cos(self.plane.sweep))) * (self.CmO_root_list + self.CmO_tip_list) / 2\
                    + self.dCmdEps[1:]*self.twist[1:]
        print("Cmac", self.Cmac)