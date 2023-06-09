from Plane import Plane, Trim,Tail
from CoefficientCalculations import AerodynamicProperties

import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.join(current_dir, "..", "Systems_and_Structures")
sys.path.append(parent_dir)
from LG_design import LandingGear
# from ..Systems_and_structures.LG_design import LandingGear

from BatterySize import Planform_calculation

import matplotlib.pyplot as plt
import numpy as np

MTOW=6521
m_payload=2596.51
m_struc=1695.5
m_eng=400
m_batprop=1829.21
m_system=20
x_cg_pylon=0.85
x_cg_lg=0.65
x_cg_vtail=0.875
x_cg_payload=0.281
x_cg_eng=0.85
x_cg_batprop=0.281
x_cg_system=0.281
Wingloading=1412
Volume=6.5

#Tail class inputs
eta =0.25
SrS =0.36 # based on cessna citation 500 - Roskam
T_engine = 5000 #[N]
d_engine = 2 #[m] distance of 1 engine from centerline


Batterysize=Planform_calculation(".\Airfoil_dat\MH 91  14.98%.dat",".\Airfoil_dat\MH 91  14.98%.dat",MTOW,Wingloading,Volume,0.25,0.4)
plane=Batterysize.makeplane()
x_cg=plane.calculate_COG(x_cg_pylon,x_cg_lg,x_cg_vtail,m_eng,x_cg_eng,m_batprop,x_cg_batprop,m_payload,x_cg_payload,m_system,x_cg_system,MTOW)
tail = Tail(plane,eta,Srs,T_engine,d_engine,x_cg)
plane.calculate_MOI(0,0,0)
LG=LandingGear(x_cg,plane.b_tot,plane.sweep[0],MTOW,plane.c[1],plane.c[0],plane.MAC,pusher=True)
plt.scatter([LG.track_width_MLG/2,0,0],[LG.pos_x_MLG,-LG.pos_x_NLG,x_cg])
plane.plot_plane()


Batterysize=Planform_calculation(".\Airfoil_dat\MH 91  14.98%.dat",".\Airfoil_dat\MH 91  14.98%.dat",MTOW,Wingloading,Volume,0.25,LG.track_width_MLG/(plane.b_tot))
plane=Batterysize.makeplane()
x_cg=plane.x_quarter
tail = Tail(plane,eta,Srs,T_engine,d_engine,x_cg)
LG=LandingGear(x_cg,plane.b_tot,plane.sweep[0],MTOW,plane.c[1],plane.c[0],plane.MAC,pusher=True)
plt.scatter([LG.track_width_MLG/2,0,0],[LG.pos_x_MLG,-LG.pos_x_NLG,x_cg])
plane.plot_plane()
plane.xflrvelues()

trim=Trim(0.5, 0.02, 0.02,110,5)




Coeff=AerodynamicProperties(plane,tail,trim,)

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
