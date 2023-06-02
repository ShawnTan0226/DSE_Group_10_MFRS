import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import interp1d

class AerodynamicProperties:
    def __init__(self,planegeometry,h=5000,V=110):
        self.h=h
        self.V=V

    ### ISA CALCULATIONS ###
    def atmos(self):
        self.T = 288.15 - 0.0065 * self.h
        self.rho = 1.225 * (self.T / 288.15) ** (-1 + 9.81 / (287.058 * 0.0065))
        self.nu = 0.0000169
    def aerodynamic_properties(self):
        self.atmos()
        self.Re = self.rho * self.V * self.MAC / self.nu
        self.Re_list = self.rho * self.V * self.MAC_list / self.nu
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
        FF = (1 + 0.6 / self.max_thickness_location * self.max_thickness + 100 * (self.max_thickness) ** (4)) * (
                    1.34 * self.M ** 0.18 * (np.cos(self.sweep[part])) ** 0.28)
        print('FF', FF)
        S_wet = 2 * self.S_list[part] * 1.07
        C_D_part = FF * C_f * S_wet / self.S
        return C_D_part
    def define_C_D_0(self, laminar_frac):
        self.C_D_0 = 0
        for i in range(len(self.taper)):
            self.C_D_0 += self.define_C_D_part_wing(laminar_frac, i)
        self.C_D_0 += 0.1 * self.C_D_0
        return self.C_D_0
    def Cm0(self):
        self.deltaCm0epsilon = -0.06
        self.epsilon = 3  # radians
        self.CmO_root = -0.6
        self.CmO_tip = -0.6
        self.Cmac = ((self.A * (np.cos(self.sweep) ** 2)) / (self.A + 2 * np.cos(self.sweep))) * (
                    self.CmO_root + self.CmO_tip) / 2 +