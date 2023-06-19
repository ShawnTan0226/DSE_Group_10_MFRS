from Plane import Plane, Trim, Tail
from CoefficientCalculations import AerodynamicProperties
from Control import control
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.join(current_dir, "..", "Systems_and_Structures")
sys.path.append(parent_dir)
from LG_design import LandingGear
# from ..Systems_and_structures.LG_design import LandingGear

from Stability import Stab
from BatterySize import Planform_calculation

import matplotlib.pyplot as plt
import numpy as np

MTOW = 5514.97
m_struc=1378.55
m_eng=203
m_batprop=1316.92
m_payload=2596.51
m_batt=m_payload+m_batprop
m_pl=0
m_system=0
Wingloading=1210.377
V_prop=3
V_pl=3.5
Volume=V_prop+V_pl


#Tail class inputs
eta = 1
SrS =0.36 # based on cessna citation 500 - Roskam
T_engine = 5600 #[N]
d_engine = 3 #[m] distance of 1 engine from centerline
l_engine = 2


plot_plane_confiig=True
plot_payload_config=True
record=False

Batterysize=Planform_calculation(".\Airfoil_dat\MH 91  14.98%.dat",".\Airfoil_dat\MH 91  14.98%.dat",MTOW,Wingloading,V_prop,V_pl,0.25,0.4)
plane=Batterysize.makeplane()
x_cg=plane.x_quarter
tail = Tail(plane,eta,SrS,T_engine,l_engine, d_engine,x_cg)
tail.tail_sizing_2()
#tail.funct_f_b_wt(2.215)

# Define the function to plot

# Generate x values
x = np.linspace(-10, 20)
# y = tail.funct_f_b_wt(x)
# plt.plot(x, y)
# plt.show()

# Generate y values by applying the function to each x value
# Create the plot
LG=LandingGear(x_cg,plane.b_tot,plane.sweep[0],MTOW,plane.c[1],plane.c[0],plane.MAC,pusher=True)

pushers=True

for i in range(10):
    Batterysize=Planform_calculation(".\Airfoil_dat\MH 91  14.98%.dat",".\Airfoil_dat\MH 91  14.98%.dat",MTOW,Wingloading,V_prop,V_pl,0.25,LG.track_width_MLG/(plane.b_tot),sweep_inner=np.deg2rad(30),sweep_outer=np.deg2rad(28))
    plane=Batterysize.makeplane()
    x_cg_batt=Batterysize.battery_placement()
    plane.x_rect_batt=Batterysize.x_rect_batt
    plane.cg_list=Batterysize.cg_list

    #x_cg_pylon=0.85*plane.b_tot/2
    #x_cg_payload=0.281*plane.b_tot/2
    #x_cg_eng=0.85*plane.b_tot/2
    #x_cg_system=0.281*plane.b_tot/2
    
    
    x_cg_pylon=plane.coords_bot[1]
    x_cg_payload=0.281*plane.b_tot/2
    x_cg_eng=plane.coords_bot[1]+0.65
    x_cg_system=0.281*plane.b_tot/2

    LG=LandingGear(x_cg,plane.b_tot,plane.sweep[0],MTOW,plane.c[1],plane.c[0],plane.MAC,pusher=pushers)
    tail = Tail(plane,eta,SrS,T_engine,l_engine,d_engine,x_cg)
    tail.tail_sizing_2()

    x_cg=plane.calculate_COG(x_cg_pylon,LG.pos_x_MLG,tail.x_tail,m_eng,x_cg_eng,m_batt,x_cg_batt,m_pl,x_cg_payload,m_system,x_cg_system,MTOW)

tail.tail_dimensions(7)


# print('Bruhhhhh',tail.S_v_b1,tail.x_offset_engine,tail.S_v_wt1)

LG=LandingGear(x_cg,plane.b_tot,plane.sweep[0],MTOW,plane.c[1],plane.c[0],plane.MAC,pusher=pushers)


# plane.plot_plane(True)

plane.calculate_MOI(0,0,0)

trim=Trim(0.5, 0.02, 0.02,110,4.75)





Coeff=AerodynamicProperties(plane,tail,trim,[[0.025,0.025],[0.025,0.025]])


Cs=control(plane,tail,0.4,-30,1.2,Coeff.Cmac,Coeff.C_L_alpha,Coeff.C_l_alpha)
print('dcm----',Cs.dCm_req)
print(Cs.calc_delta_Cm_outer())
print(Cs.eta_i)
print('new_CL_max',Cs.CL_max_new)

coord_y_cs=np.array([((plane.b[1]+1.3*2)/plane.b[2]+0.05)*plane.b[2]/2,Cs.eta_i*plane.b[2]/2,0.95*plane.b[2]/2])
chord_cs=2*(plane.c[2]-plane.c[1])/(plane.b[2]-plane.b[1])*(coord_y_cs-plane.b[1]/2)+plane.c[1]
coord_x_cs=2*(plane.coords_bot[2]-plane.coords_bot[1])/(plane.b[2]-plane.b[1])*(coord_y_cs-plane.b[1]/2)+plane.coords_bot[1]
coord_top_x_cs=coord_x_cs-chord_cs*0.4
coord_full_y_cs=np.concatenate((coord_y_cs,coord_y_cs[::-1]))
coord_full_x_cs=np.concatenate((coord_top_x_cs,coord_x_cs[::-1]))

coord_y_cs_plot=np.copy(coord_y_cs)
coord_y_cs_plot=np.insert(coord_y_cs_plot,1,coord_y_cs[1])
coord_y_cs_plot=np.insert(coord_y_cs_plot,1,coord_y_cs[1])

coord_x_cs_plot=np.copy(coord_top_x_cs)
coord_x_cs_plot=np.insert(coord_x_cs_plot,1,coord_x_cs[1])
coord_x_cs_plot=np.insert(coord_x_cs_plot,1,coord_top_x_cs[1])

coord_full_y_cs_plot=np.concatenate((coord_y_cs_plot,coord_y_cs[::-1],[coord_y_cs_plot[0]]))
coord_full_x_cs_plot=np.concatenate((coord_x_cs_plot,coord_x_cs[::-1],[coord_x_cs_plot[0]]))


if plot_plane_confiig:
    if plot_payload_config:
        Wing_batt_color=(0/255, 118/255, 194/255)
        Wing_batt_a=1
        Body_batt_color=(0/255, 184/255, 200/255)
        Body_batt_a=1
        engine_color=(224/255, 60/255, 49/255)
        engine_a=1
        computer_color=(237/255, 104/255, 66/255)
        computer_a=1
        other_color=(220/255, 220/255, 220/255)
        other_a=1
        plane.draw_battery_placement(Wing_batt_color,Body_batt_color,engine_color,computer_color,other_color,Wing_batt_a,Body_batt_a,engine_a,computer_a,other_a,False)
    else:
        plane.plot_plane(False)
    plt.scatter([-LG.track_width_MLG/2,LG.track_width_MLG/2,0],[LG.pos_x_MLG,LG.pos_x_MLG,LG.pos_x_NLG],color="blue",label="LG",zorder=8)
    plt.scatter([0],[x_cg],color="red",label="CG",zorder=8)
    plt.scatter([plane.y_quarter],[plane.x_quarter],color="black",label="ac",zorder=8)
    print(tail.cr_v_b,plane.coords_bot[0])
    plt.plot([0,0],[plane.coords_bot[0]-tail.cr_v_b,plane.coords_bot[0]],color="brown",label="tail")

    plt.plot(coord_full_y_cs_plot,coord_full_x_cs_plot,color="black",linewidth=0.75)
    plt.plot(-coord_full_y_cs_plot,coord_full_x_cs_plot,color="black",linewidth=0.75)
    plt.fill(coord_full_y_cs,coord_full_x_cs,color="orange",label="control surface",alpha=1)
    plt.fill(-coord_full_y_cs,coord_full_x_cs,color="orange",alpha=1)

    plt.legend()
    plt.axis('equal')
    print('LG height',LG.height_MLG)
    print('x_cg',x_cg)
    plt.show()
    plane.xflrvelues()

Stability=Stab(plane,Coeff,MTOW)
Stability.get_asymm_eigen()
Stability.get_symm_eigen()
Stability.halftimes()
Stability.damping()
print(Stability.eigenvalues_symm,Stability.eigenvalues_asymm)
print('A---',Stability.A_asymm,'\nB---',Stability.B_asymm,'\nC---',Stability.C_asymm,'\nD---',Stability.D_asymm,'\nE---',Stability.E_asymm)
print(Stability.Routh_discriminant())
print('damping asymm: ',Stability.damping_asymm,'damping symm',Stability.damping_symm)
print('halftime asymm: ',Stability.halftime_asymm,'halftime symm',Stability.halftime_symm)


if record:
    print(Cs.eta_i)
    print(LG.track_width_MLG,LG.height_MLG,LG.pos_x_MLG,LG.pos_x_NLG)
    plane.record_planform()
    Stability.record_stability()
    tail.record_tail()

'''
plane=Plane(10,[0.4,0.5],[30, 30],[10,20])

plane.MAC=4.04
plane.S=81
Coeff=AerodynamicProperties(plane,1,1,1)
Coeff.define_C_D_part_nacelle(0.2)
plane=Plane(10,[0.4,0.5],[30, 30],[10,20],[0,5])
print(plane.y_quarter,plane.x_quarter,plane.MAC,plane.b,plane.MAC_list,plane.y_list)
plane.plot_plane()
plane.calculate_COG(pylon_cg=0.85, lg_cg=0.65, vertical_tail_cg=0.875, engine_mass=1224, engine_cg=0.85, battery_mass=10600, battery_cg=0.281, payload_mass=2903, payload_cg=0.281, system_mass=300, system_cg=0.281, MTOW=19900)
'''