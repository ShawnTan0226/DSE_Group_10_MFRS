import numpy as np
import matplotlib.pyplot as plt
import Int_Force_diagrams as ifd

n = 10000 # resolution of data along half span
Sweep05 = np.deg2rad(35.36)

# Weights in the aircraft
BatteryMass = 2903 + 10600 # kg from 2nd weight estimation
StructuralMass = 5174 # kg from 2nd weight estimation
EngineMass = 1224 # kg for 7 engines, evenly spread¨
EngineMass = EngineMass / 7
ThrustTO = 6000 # N thrust of one engine at Lift Off
ThrustCruise = 4288.32 # N thrust of one engine at cruise
# Take off force at hook
acceleration = 2.84
mass_total = BatteryMass + StructuralMass + EngineMass
ForceHookTO = acceleration * mass_total

#CruiseDistribution, TODistribution, y_points = ifd.importdat(n)
#Lcruise, Ltakeoff, Tcruise, W, M_cruise, M_TO, Dcruise, Dtakeoff = ifd.forces(CruiseDistribution, CruiseDistribution, y_points, n)
#Vx, Vy, Mx, My, Mz = ifd.InternalLoads(Lcruise, Tcruise, W, Dcruise, M_TO, n, y_points, CruiseDistribution['y-span'][-1])
#Ixx, Iyy, ts = ifd.PrelimSizing(Vx, Vy, Mx, My, Mz, CruiseDistribution)

def TestCase1():
    # values with calculated on paper
    b = 10 # halfspan
    W = 1000 * 9.81 # half weight in N
    n = 10000
    sweep = np.deg2rad(45)
    Zaxis = np.arange(0,10,b/n)
    lf = 1 # loadfactor
    Ltot = W
    Ldis = Ltot/b
    W1 = W * 0.8
    W2 = W * 0.2
    Ddis = Ldis/10
    Mdis = 200
    L = np.ones(n) * Ldis
    D = np.ones(n) * Ddis
    M = np.ones(n) * Mdis
    W = np.zeros(n)
    W[0] = W1 * b/n
    W[int(round(n/2))] = W2 / (b/n)
    T = np.zeros(n)
    T[int(round(n/2))] = Ltot/10 / (b/n)
    chord = 1 # meter
    chord = np.ones(n) * chord

    # simulation with the distribution
    Vx, Vy, Mx, My, Mz = ifd.InternalLoads(L, T, W, D, M, n, Zaxis, b, sweep)
    Ixx, Iyy, ts = ifd.PrelimSizing(Vx, Vy, Mx, My, Mz, chord)

    print('Test results')
    print('------------------')
    print('Expected Vymax = -7848 N. Simulated Vymax = ', Vy[0])
    print('Expected Vy(z=5m) = -2943 N. Simulated Vy(z=5m) = ', Vy[int(round(n/2))-1])
    print('Expected Mxmax = 39240 Nm. Simulated Mxmax = ', Mx[0])
    print('Expected Mx(z=5m) = 12262 Nm. Simulated Mx(z=5m) = ', Mx[int(round(n / 2)) - 1])
    print('------------------')
    print('Expected Vx(z=0m) = 0 N. Simulated Vx(z=0m) = ', Vx[0])
    print('Expected Vx(z=5m) = +-490.5 N. Simulated Vx(z=5m) = ', Vx[int(round(n/2))-1])
    print('Expected My(z=0m) = 0 Nm. Simulated My(z=0m) = ', My[0])
    print('Expected My(z=5m) = 1226.2 Nm. Simulated My(z=5m) = ', My[int(round(n / 2)) - 1])
    print('------------------')
    print('Internal torque with sweep = 45°')
    print('Expected Mz(z=0m) = 37234 Nm. Simulated Vx(z=0m) = ', Mz[0])
    print('Expected Mz(z=5m) = 11262.5 Nm. Simulated Vx(z=5m) = ', Mz[int(round(n / 2)) - 1])

    #plt.plot(Zaxis, Vx, label='Vx')
    #plt.plot(Zaxis, Vy, label='Vy')
    #plt.plot(Zaxis, Mx, label='Mx')
    #plt.plot(Zaxis, My, label='My')
    #plt.plot(Zaxis, Mz, label='Mz')

    plt.plot(Zaxis, Ixx, label='Ixx')
    plt.plot(Zaxis, Iyy, label='Iyy')
    #plt.plot(Zaxis, ts, label='thickness of spar web')

    plt.legend()
    plt.show()


TestCase1()