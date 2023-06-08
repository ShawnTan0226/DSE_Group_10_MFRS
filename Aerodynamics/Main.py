from Plane import Plane
from CoefficientCalculations import AerodynamicProperties

import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.join(current_dir, "..", "Systems_and_Structures")
sys.path.append(parent_dir)
from LG_design import LandingGear

from BatterySize import Planform_calculation

LandingGear(10,10,10,10,10,10,10,10)
# MakePlane = Planform_calculation(".\Airfoil_dat\MH 91  14.98%.dat",".\Airfoil_dat\MH 91  14.98%.dat",19902,2409,7.5,0.25,0.4)
# MakePlane.makeplane()
# MakePlane.plane.xflrvelues()
# MakePlane.plane.define_C_D_0(0.2)
# MakePlane.plane.plot_plane()
# print(MakePlane.plane.MAC)
# print(MakePlane.plane.C_D_0)

# MakePlane = Planform_calculation(".\Airfoil_dat\MH 91  14.98%.dat",".\Airfoil_dat\MH 91  14.98%.dat",19902,2409,20.5,0.25,0.4)
# MakePlane.makeplane()
# MakePlane.plane.xflrvelues()
# MakePlane.plane.define_C_D_0(0.2)
# MakePlane.plane.plot_plane()
# print(MakePlane.plane.MAC)
# print(MakePlane.plane.C_D_0)
Batterysize=Planform_calculation(".\Airfoil_dat\MH 91  14.98%.dat",".\Airfoil_dat\MH 91  14.98%.dat",19902,2409,7.5,0.25,0.4)
plane=Batterysize.makeplane()
Coeff=AerodynamicProperties(plane,1,1,1)

plane=Plane(10,[0.4,0.5],[30, 30],[10,20])

plane.MAC=4.04
plane.S=81
Coeff=AerodynamicProperties(plane,1,1,1)
Coeff.define_C_D_part_nacelle(0.2)
plane=Plane(10,[0.4,0.5],[30, 30],[10,20],[0,5])
print(plane.y_quarter,plane.x_quarter,plane.MAC,plane.b,plane.MAC_list,plane.y_list)
plane.plot_plane()
plane.calculate_COG(pylon_cg=0.85, lg_cg=0.65, vertical_tail_cg=0.875, engine_mass=1224, engine_cg=0.85, battery_mass=10600, battery_cg=0.281, payload_mass=2903, payload_cg=0.281, system_mass=300, system_cg=0.281, MTOW=19900)
print(f"x_cg: {plane.x_cg}, x_relative_cg: {plane.x_relative_cg}")
print(f"total x_ac: {plane.x_quarter}, x_ac_rel: {plane.x_quarter/plane.length}")
print(f"x_abs_dif: {plane.x_quarter - plane.x_cg}")

print('x-list', plane.x_list)
print('MAC_list',plane.MAC_list)
print('offset',plane.offset)

print('MAC_eq',plane.MAC_eq)
print('MAC_calc',plane.MAC)
print('S',plane.S)
print('S_eq',plane.S_eq)
# print('S_eq2',plane.S_eq2)
print('cr',plane.cr_eq)
print('ct',plane.ct_eq)
print("A", plane.A)
print('b', plane.b, plane.b_tot)
# # test.plot_plane()
# test.xflrvelues()
# # test.drawbox(0.5)
# test.drawtail(0.2)
# print(test.MAC_aircraft())

#print(f"a= {plane.c}")
#print(plane.determine_CoG(engine_cg=0.85, battery_cg=0.27, payload_cg=0.27))

'''Inputs'''
# #Aerodynamic coefficients
# dCmdEps = [0,0.4,0.4]
# twist = [0,0.1,0.1]
# CmO_airfoil = [[0.1,0.05],[0.2,0.1]]
#
# coefficients = AerodynamicProperties(plane, dCmdEps, twist, CmO_airfoil)
# print(coefficients.Cmac())