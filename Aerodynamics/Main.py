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



test=Plane(9.72,[0.3,0.267],[38,38],[8.82,22])
test.plot_plane()
print(test.COG(pylon_cg=0.85, lg_cg=0.65, vertical_tail_cg=0.875, engine_mass=1244, engine_cg=0.85, battery_mass=10600, battery_cg=0.281, payload_mass=2903, payload_cg=0.28, system_mass=300, system_cg=0.281, MTOW=19900))
# # test.plot_plane()
# test.xflrvelues()
# # test.drawbox(0.5)
# test.drawtail(0.2)
# print(test.MAC_aircraft())


