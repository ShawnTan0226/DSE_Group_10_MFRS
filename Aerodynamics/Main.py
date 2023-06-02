from Plane import Plane
from BatterySize import Planform_calculation

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



plane=Plane(9.72,[0.3,0.267],[38,38],[8.82,22])
# # test.plot_plane()
# test.xflrvelues()
# # test.drawbox(0.5)
# test.drawtail(0.2)
# print(test.MAC_aircraft())

#print(f"a= {plane.c}")
#print(plane.determine_CoG(engine_cg=0.85, battery_cg=0.27, payload_cg=0.27))

'''Inputs'''
#Aerodynamic coefficients
dCmdEps = [0,0.4,0.4]
twist = [0,0.1,0.1]
CmO_airfoil = [[0.1,0.05],[0.2,0.1]]

coefficients = AerodynamicProperties(plane, dCmdEps, twist, CmO_airfoil)
print(coefficients.Cmac())