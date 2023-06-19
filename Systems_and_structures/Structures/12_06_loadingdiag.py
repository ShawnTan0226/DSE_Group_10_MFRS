import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import integrate

n = 10 # resolution of data along half span

Sweep05 = np.deg2rad(35.36)

# Weights in the aircraft
BatteryMass = 1829 + 2596 # kg from 2nd weight estimation
StructuralMass = 1695 # kg from 2nd weight estimation
EngineMass = 0#400 # kg for 7 engines, evenly spread
EngineMass = EngineMass / 2
ThrustTO = 0 # N thrust of one engine at Lift Off
ThrustCruise = 3011 # N thrust of one engine at cruise
# Take off force at hook
acceleration = 0

mass_total = BatteryMass + StructuralMass + EngineMass
ForceHookTO = 0 # acceleration * mass_total - ThrustTO * 2
#print(ForceHookTO)

##### Import planform and aerodynamic data #####

def importdat(n):
    CruiseDistribution = '4.75deg cruise lift dist.csv'
    data_cruise = pd.read_csv(CruiseDistribution, index_col=None)
    data_cruise = data_cruise.to_dict(orient='list')
    y_old =  data_cruise['y-span']

    # Interpolating to n data points with a new y axis

    y_points = np.linspace(data_cruise['y-span'][0], data_cruise['y-span'][-1], n)
    y_points = np.linspace(data_cruise['y-span'][0], data_cruise['y-span'][-1], n)
    for i in data_cruise.keys():
        data_cruise[i] = np.interp(y_points, y_old, data_cruise[i])


    return data_cruise, y_points

##### Define the force distributions along the beam (half span) #####

def forces(CruiseDistribution, y_points, n):

    # Lift and drag at cruise :
    Cl = CruiseDistribution['Cl']
    Cdi = CruiseDistribution['ICd'] + 0.012# total cd cdi + cd0
    c = CruiseDistribution['Chord']
    V = 110 # m/s
    rho = 0.7361 # kg/m^3
    dy = (CruiseDistribution['y-span'][-1] - CruiseDistribution['y-span'][0]) / n
    Lcruise = Cl * 0.5 * rho * V**2 * c # N/m
    Dcruise = Cdi * 0.5 * rho * V**2 * c # N/m

    print(np.average(c),n)
    # Weight distribution along span
    # distributed weights of batteries, payload and structure
    Wave = (StructuralMass + BatteryMass)/2 * 9.81 / CruiseDistribution['y-span'][-1] # N/m
    W = Wave * CruiseDistribution['Chord'] / ((CruiseDistribution['Chord'][0] + CruiseDistribution['Chord'][-1]) / 2) # N/m mass spread normalised with chord distribution.
    #W[0] += EngineMass * 9.81 / (CruiseDistribution['y-span'][-1] / n)
    #W[-1] += EngineMass*2 * 9.81 / (CruiseDistribution['y-span'][-1] / n)
    #W[int(round(n/3))] += EngineMass * 9.81 / (CruiseDistribution['y-span'][-1] / n)
    W[int(round(n*(1.47/CruiseDistribution['y-span'][-1])))] += EngineMass * 9.81 / (CruiseDistribution['y-span'][-1] / n)


    T = np.zeros(n)
    T[int(round(n*(1.47/CruiseDistribution['y-span'][-1])))] += ThrustCruise
    TC = T / (CruiseDistribution['y-span'][-1] / n)

    # Moment distribution :
    V = 110 # m/s
    rho = 0.7361 # kg/m^3
    M_cruise = CruiseDistribution['CmAirf@chord/4'] * 1/2 * rho * V**2 * CruiseDistribution['Chord']

    return Lcruise, TC, W, M_cruise, Dcruise,

def TestForces(Lcruise, Tcruise, W):
    dy = (CruiseDistribution['y-span'][-1] - CruiseDistribution['y-span'][0]) / n
    c = CruiseDistribution['Chord']

def InternalLoads(L, T, W, D, M, n, y_points, halfspan,sweep):
    b = halfspan * 2
    Dtot = T - D # drag and thrust act on the x axis
    Vx = integrate.cumtrapz(np.flip(Dtot * b / (2 * n)))[::-1]
    Vy = integrate.cumtrapz(np.flip((-L + W) * b / (2 * n)))[::-1]
    Vx = np.append(Vx,[0])
    Vy = np.append(Vy,[0])

    My = -integrate.cumtrapz(np.flip(Vx * b / (2 * n)))[::-1]
    Mx = -integrate.cumtrapz(np.flip(Vy * b / (2 * n)))[::-1]
    My = np.append(My,[0])
    Mx = np.append(Mx,[0])

    Ml = []
    Mw =[]
    yp = y_points
    count = 2

    # due to lift and weight moment arm
    for i in yp[1:]:
        Mli = IntegrateTorqueFromLift(count, yp, -L, sweep)
        Mwi = IntegrateTorqueFromLift(count, yp, W, sweep)
        Ml.append(Mli)
        Mw.append(Mwi)
        count += 1

    # due to aerodynamic moment:
    Mzm = integrate.cumtrapz(np.flip(M), np.flip(yp))
    Mz = (np.array(Ml) + np.array(Mw) + np.array(Mzm))[::-1]
    Mz = np.append(Mz, [0])

    return Vx, Vy, Mx, My, Mz

def IntegrateTorqueFromLift(c, axis, data, sweep):
    data = np.flip(data)
    data = data[:c]
    axis = np.flip(axis)
    axis = axis[:c]
    data = data * np.tan(sweep) * (axis - axis[c-1])

    M = np.trapz(data, axis)
    return M

def PrelimSizing(Vx, Vy, Mx, My, Mz, chords):
    # Material properties
    sigma = 324000000 # Pa yield strength 2024 Al
    tau = 283000000 # Pa yield shear strength 2024 Al
    tc = 0.15 # thickness over chord ratio of airfoil
    t = chords * tc

    Ixx = abs(Mx) * t / (2 * sigma)
    Iyy = abs(My) * chords / (2 * sigma)
    # assuming two area points at the thickness extremities:
    A = Ixx * 4 / t ** 2
    # determining the necessary spar thickness for shear stresses due to internal shear force:
    Q = 2 * A * 0.5 * t
    ts = abs(Vy) * Q / (Iyy * tau) # spar thickness needed to withstand shear stress due to bending

    # thickness of wingbox for torque
    twb = abs(Mz)/(2 * chords**2 * (0.4*0.15) * tau)

    #plt.plot(y_points, A)
    #plt.plot(y_points, Ixx)
    #plt.show()
    return Ixx, Iyy, ts, twb

nl = -1 # loadfactor
CruiseDistribution, y_points = importdat(n)
Lcruise, TC, W, M_cruise, Dcruise = forces(CruiseDistribution, y_points, n)
Vx, Vy, Mx, My, Mz = InternalLoads(nl*Lcruise, TC, W, abs(nl)*Dcruise, nl*M_cruise, n, y_points, CruiseDistribution['y-span'][-1], Sweep05)
#Ixx, Iyy, ts, twb = PrelimSizing(Vx, Vy, Mx, My, Mz, CruiseDistribution['Chord'])

def truss(Force, n, ydistr):
    truss_nodes = [0,0.226917131,0.453834262,0.680751393,0.907668524,1.134585655,1.361502786,1.588419917,1.815337048,2.042254179,2.26917131,2.496088441,2.723005572,2.949922703,3.176839834,3.403756965,3.630674096,3.857591227,4.084508358,4.311425489,4.53834262,4.765259751,4.992176882,5.219094013,5.446011144,5.672928275,5.899845406,6.126762537,6.353679668,6.580596799,6.80751393,7.034431061,7.261348192,7.488265323,7.715182454,7.942099585,8.169016716,8.24]
    force_nodes = np.interp(truss_nodes, ydistr, Force)
    len_nodes = (np.pad(truss_nodes, (0, 1), 'constant') - np.pad(truss_nodes, (1, 0), 'constant'))[1:]
    force_nodes = force_nodes * len_nodes / 4
    force_nodes[-1] = Force[-1] * (truss_nodes[-1]-truss_nodes[-2])/2 # holds force on only one side
    force_nodes[0] = Force[0] * (truss_nodes[1]-truss_nodes[0]) # holds force from 2 sides
    print(repr(force_nodes))
    plt.scatter(truss_nodes, 1000* len_nodes)
    return truss_nodes, force_nodes

ytruss, forcetruss = truss(-W * 3.08 *1.3, 20, y_points)

def plots():
    plt.plot(y_points, Lcruise - W, label='Lift - Weight')
    plt.scatter(ytruss, forcetruss, label='Force applied at nodes')
    #plt.plot(y_points, -W, label='Ixx')
    #plt.plot(y_points, Mx, label='Mx')
    #plt.plot(y_points, , label='tot')
    #plt.plot(y_points, My, label='My')
    #plt.plot(y_points, Mz, label='Mz')
    plt.xlabel("yb [m]")
    plt.ylabel("-F_z [N], -w_z [N]")
    plt.grid()
    plt.legend()
    plt.show()
    np.savetxt("output_load_nodes.csv", forcetruss, delimiter=",")

plots()