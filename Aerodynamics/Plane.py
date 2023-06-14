import numpy as np
import math as math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import interp1d
#TODO: Implement multiple airfoil
class Plane:
    def __init__(self,Cri,taper,sweep,b,Volume,twist=[0,0],dihedral=0,planename='Main.csv',ip=0,h=5000,V=110,airfoil=".\Airfoil_dat\MH 91  14.98%.dat", number_of_tail=1):
        #Plane object has n sections
        self.c=np.array([Cri]) #array of chord [m] form: [Middle c,c1,c2,c3,...,cn]
        self.taper = np.array(taper) #array of taper ratio form: [taper1,taper2,taper3,...,tapern]
        self.sweep = np.array(np.deg2rad(sweep)) #array of sweep angle [rad] form: [sweep1,sweep2,sweep3,...,sweepn]
        self.twist = np.array(np.deg2rad(twist))
        self.dihedral = np.deg2rad(dihedral)
        self.b = np.concatenate(([0],b)) #array of span [m] form: [0,b1,b2,b3,...,bn]
        self.S_list=np.array([]) #array of surface area of each section [m^2] form: [S1,S2,S3,...,Sn]
        self.battery_volume=Volume #volume of the plane [m^3]
        
        self.b_tot=self.b[-1] #total span [m]
        self.coords=np.array([])
        self.V=V
        self.h=h

        self.planename=planename
        self.ip=ip
        self.cg_list=0

        self.MAC_list=np.array([])
        self.x_list=np.array([])
        self.y_list = np.array([])

        self.draw()
        self.MAC_aircraft()
        self.equivalent_wing()
        self.define_airfoil(airfoil)


        self.tailnumber = number_of_tail

        # self.aerodynamic_properties()

    def help(self):
        print('Visual of top view of plane with plot_plane() \n ',
              'Get the xflr inputs with xflrvalues() \n ',
              'Plot plane with payload area with drawbox(opacity) \n ',
              'Plot plane with tail area with drawtail(opacity) \n ',
              'Get the MAC and x_quarter with MAC_aircraft() \n ',
              'Get the C_D_0 with define_C_D_0(laminar_frac) \n ',)


    ### WING GEOMETRY ###
    def draw(self):
        self.offset=np.array([0])
        count=0
        for i in self.taper:
            self.c=np.concatenate((self.c, [i*self.c[-1]]))
            nextoffset=self.offset[-1]+np.tan(self.sweep[count])*(self.b[count+1]-self.b[count])/2+0.25*self.c[-2]-0.25*self.c[-1]
            self.offset=np.concatenate((self.offset, [nextoffset]))
            self.S_list=np.concatenate((self.S_list,[(self.c[-1]+self.c[-2])/2*(self.b[count+1]-self.b[count])]))
            count+=1


        self.S=np.sum(self.S_list)
        self.A = self.b[-1]**2/self.S
        self.A_list = (self.b[1:]-self.b[:-1])**2/self.S_list

        self.coords_bot=self.offset+self.c
        self.coords=np.concatenate((self.coords,self.offset,self.coords_bot[::-1]))
        self.coords=np.concatenate((self.coords,self.coords[::-1]))
        negative=np.concatenate((-self.b,-self.b[::-1]))[::-1]


        self.bfull=np.concatenate((self.b/2,self.b[::-1]/2,negative/2))

    def plot_plane(self):
        plt.plot(self.bfull,self.coords,color='black')
        plt.fill(self.bfull,self.coords, color='gray', alpha=0.5)
        # plt.plot(np.concatenate((self.y_list,self.y_list[::-1],[self.y_list[0]])),np.concatenate((self.x_list-0.25*self.MAC_list,self.x_list[::-1]+0.75*self.MAC_list[::-1],[self.x_list[0]-0.25*self.MAC_list[0]])),color='red')
        # plt.plot([self.y_quarter,self.y_quarter],[self.x_quarter-0.25*self.MAC,self.x_quarter+0.75*self.MAC])
        # plt.scatter(self.y_list,self.x_list)
        # plt.plot([0,self.b[-1]/2,self.b[-1]/2,0],[self.x_cr_eq,self.x_ct_eq,self.x_ct_eq+self.ct_eq,self.x_cr_eq+self.cr_eq],color='blue',linestyle='--')
        plt.gca().invert_yaxis()
        plt.show()
    def add_text_to_file(self,file_path, text):
        with open(file_path, 'a') as file:
            file.write(text)
    
    def record_planform(self):
        text='S: '+str(self.S)+'\nb: '+str(self.b)+'\nc: '+str(self.c)+'\ntaper: '+str(self.taper)+'\nsweep: '+str(self.sweep)+'\nA: '+str(self.A)+'\nMAC: '+str(self.MAC)+'\n\n'
        print(text)
        self.add_text_to_file('./Record/Planform record.txt', text)

    def get_planform(self):
        print('S: ',self.S)
        print('b: ',self.b)
        print('c: ',self.c)
        print('taper: ',self.taper)
        print('sweep: ',self.sweep)


    def xflrvelues(self):
        xflr={}
        xflr['c']=self.c
        xflr['b']=self.b
        xflr['offset']=self.offset
        print('c: ',self.c)
        print('b: ',self.b)
        print('offset: ',self.offset)
        return xflr
    def drawbox(self,opacity):
        frontbox=self.offset+0.15*self.c
        backbox=self.offset+0.65*self.c
        x_store=np.concatenate((frontbox,backbox[::-1]))
        x_store=np.concatenate((x_store,x_store[::-1]))
        y=np.concatenate((self.b,self.b[::-1]))
        y=np.concatenate((y,-y[::-1]))

        x_red=np.concatenate((self.offset,frontbox[::-1]))
        x_red=np.concatenate((x_red,x_red[::-1]))

        x_red2=np.concatenate((self.offset+self.c,backbox[::-1]))
        x_red2=np.concatenate((x_red2,x_red2[::-1]))

        plt.plot(self.bfull,self.coords, color='black')
        plt.gca().invert_yaxis()
        plt.fill(y,x_store, color='blue', alpha=opacity, label='Battery')
        plt.fill(y,x_red, color='orange', alpha=opacity)
        plt.fill(y,x_red2, color='orange', alpha=opacity)
        plt.show()

    def drawtail(self,opacity):
        x_front=self.offset[-2:]
        x_back=self.coords_bot[-2:]
        x=np.concatenate((x_front,x_back[::-1]))
        y=np.concatenate((self.b[-2:],self.b[-2:][::-1]))
        negy=-y

        plt.plot(self.bfull,self.coords, color='black')
        plt.plot([self.b[-2],self.b[-2]],[self.offset[-2],self.coords_bot[-2]],  color='black')
        plt.plot([-self.b[-2],-self.b[-2]],[self.offset[-2],self.coords_bot[-2]],  color='black')

        plt.gca().invert_yaxis()
        plt.fill(self.bfull,self.coords, color='blue', alpha=opacity)
        plt.fill(y,x, color='orange', alpha=1)
        plt.fill(negy,x, color='orange', alpha=1)
        plt.show()
        
    def MAC_part(self,cr, ct, sweep, b):
        #MAC
        y = 2/3 * (cr + ct + ct**2/cr)/(1+ct/cr)
        #offset from that parts root on quarter chord
        off_y = -(b/2)*(y-cr)/(cr-ct) #b is the span of that part of the wing
        off_x = off_y*np.tan(sweep)
        return y, off_x, off_y

    def listgenerator(self):
        for i in range(len(self.taper)):
            part= self.MAC_part(self.c[i],self.c[i+1],self.sweep[i],self.b[i+1]-self.b[i])
            self.MAC_list=np.concatenate((self.MAC_list,[part[0]]))
            self.x_list=np.concatenate((self.x_list,[part[1]]+self.offset[i]+0.25*self.c[i]))
            self.y_list = np.concatenate((self.y_list, [part[2]]+self.b[i]/2))

    def MAC_aircraft(self):#Make sure to use numpy array
        self.listgenerator()
        self.MAC = np.sum(self.MAC_list*self.S_list)/np.sum(self.S_list)
        self.x_quarter = np.sum(self.x_list*self.S_list)/np.sum(self.S_list) #X is from leading point
        self.y_quarter = np.sum(self.y_list*self.S_list)/np.sum(self.S_list)

    def equivalent_wing(self):
        self.sweep_eq = np.arctan((self.x_list[1] - self.x_list[0])/(self.y_list[1]-self.y_list[0])) #rad
        self.cr_eq = (self.MAC_list[1]-self.MAC_list[0])/(self.y_list[1]-self.y_list[0])*(-self.y_list[0])+self.MAC_list[0]
        self.ct_eq = (self.MAC_list[1]-self.MAC_list[0])/(self.y_list[1]-self.y_list[0])*(self.b[2]/2-self.y_list[0])+self.MAC_list[0]
        self.MAC_eq = self.MAC_part(self.cr_eq,self.ct_eq,self.sweep_eq,self.b[2])
        self.S_eq = (self.cr_eq+self.ct_eq)/2*self.b[2]
        self.x_cr_eq = self.x_list[0]-np.tan(self.sweep_eq)*self.y_list[0]-0.25*self.cr_eq
        self.x_ct_eq = self.x_list[0]+np.tan(self.sweep_eq)*(self.b[-1]/2-self.y_list[0])-0.25*self.ct_eq

        self.sweep_eq_half = np.arctan(((self.x_ct_eq+0.5*self.ct_eq)-(self.x_cr_eq+0.5*self.cr_eq))/(self.b[2]/2))
        self.taper_eq = self.ct_eq/self.cr_eq


    def define_airfoil(self,file_path):
        # Initialize empty arrays for positive and negative values
        positive_column1 = []
        positive_column2 = []
        negative_column1 = []
        negative_column2 = []

        positive_column11 = []
        positive_column21 = []
        negative_column11 = []
        negative_column21 = []

        # Read the .dat file
        with open(file_path, 'r') as file:
            for line in file:
                # Remove leading/trailing whitespaces and split the line by spaces
                data = line.strip().split()
                if len(data) >= 2:  # Ensure the line has at least two columns
                    try:
                        value1 = float(data[0])
                        value2 = float(data[1])
                        if 0.15 <= value1 <= 0.55:
                            if value2 >= 0:
                                positive_column1.append(value1)
                                positive_column2.append(value2)
                            else:
                                negative_column1.append(value1)
                                negative_column2.append(value2)
                            
                    except ValueError:
                        continue
                    try:
                        value1 = float(data[0])
                        value2 = float(data[1])
                        if value2 >= 0:
                            positive_column11.append(value1)
                            positive_column21.append(value2)
                        else:
                            negative_column11.append(value1)
                            negative_column21.append(value2)
                            
                    except ValueError:
                        continue
        # Compute the surface using numpy
        postive_surface = -np.trapz(positive_column2, positive_column1)
        negative_surface = -np.trapz(negative_column2, negative_column1)
        Area = (negative_surface+postive_surface)

        self.eq_positive=interp1d(positive_column11[::-1],positive_column21[::-1],bounds_error=False,fill_value='extrapolate')
        self.eq_negative=interp1d(negative_column11,negative_column21,bounds_error=False,fill_value='extrapolate')
        thickness=self.eq_positive(np.arange(0,1,0.001))-self.eq_negative(np.arange(0,1,0.001))
        self.max_thickness = np.max(thickness)
        self.max_thickness_location = np.argmax(thickness)*0.001
        return Area


    def calculate_COG(self, pylon_cg, lg_cg, vertical_tail_cg, engine_mass, engine_cg, battery_mass, battery_cg, payload_mass, payload_cg, system_mass, system_cg, MTOW=19900):
        """Steps to determine cog wing structure :
        1. Partition wing into sections with each section having the corresponding c
        2. Estimate x according to sweep (determine offset)
        3. Set up equation (include the offset (due to sweep) already for the x of each section)"""
        #Step1
        self.lg_cg = lg_cg
        self.pylon_cg = pylon_cg
        self.vertical_tail_cg = vertical_tail_cg

        self.battery_density=battery_mass/self.battery_volume

        chord_body_section = np.linspace(self.c[0], self.c[1], 100)
        chord_y_body = np.linspace(0, self.b[1]/2, 100)

        chord_wing_sections = np.linspace(self.c[1], self.c[2], 100)
        chord_y_wing = np.linspace(0, self.b[2]/2 - self.b[1]/2, 100)

        #Step2 (estimates x according to sweep)
        x_body = 0.25 * self.c[0] + np.tan(self.sweep[0]) * chord_y_body + 0.1 * chord_body_section
        x_wing = (0.25 * self.c[0] + np.tan(self.sweep[1]) * self.b[1]/2) + np.tan(self.sweep[1]) * chord_y_wing + 0.1 * chord_wing_sections

        body_integrand = x_body * chord_body_section**2
        body_integration = np.trapz(body_integrand, chord_y_body)

        wing_integrand = x_wing * chord_wing_sections**2
        wing_integration = np.trapz(wing_integrand, chord_y_wing)

        #step 3
        body_wing_cg = (2*body_integration + wing_integration) / (np.trapz(chord_wing_sections**2, chord_y_wing) + 2 * np.trapz(chord_body_section**2, chord_y_body))
        """"----------------wing structure section---------------------"""


        #----------------total length of drone according to sweep----------
        total_length = max((0.25 * self.c[0] + np.tan(self.sweep[0]) * 0.5*self.b[1] + np.tan(self.sweep[1]) * (0.5*self.b[2] - 0.5*self.b[1]) + 0.75 * self.c[2]), self.c[0])

        

        #-----------COG DETERMINATION with all other subsystems------------------------------
        self.MTOW=MTOW
        total_structural_mass = 0.26 * MTOW

        #Defining each subsystem

        #--HERE: system structure is a part of body--
        pylon_mf = 0.035
        # pylon_cg = 0.85 #VARIABLE
        pylon_mass = pylon_mf * total_structural_mass #May be Variable

        lg_mf = 0.108
        # lg_cg = 0.65 #VARIABLE
        lg_mass = lg_mf * total_structural_mass #may be VARIABLE

        vertical_tail_mf = 0.024
        # vertical_tail_cg = 0.875 #VARIABLE
        vertical_tail_mass = vertical_tail_mf * total_structural_mass #may be variable

        #HERE: system mass is part of the bodywing mass
        body_wing_mass = total_structural_mass - pylon_mass - lg_mass - vertical_tail_mass - system_mass
        body_wing_cg_relative = body_wing_cg 

        #----------INPUT VARIABLES---------------
        # engine_mass = 1244.168
        # engine_cg = 0.85

        # battery_mass = 10600.13 #battery for drone
        # battery_cg = 0.2814

        # payload_mass = 2903.2 #battery for ac
        # payload_cg = 0.2814

        # system_mass = ??
        # system_cg = 0.2814
        #---------------------------------------

        #COG EQUATION
        x_relative_cg = (body_wing_cg_relative * body_wing_mass + pylon_mass * pylon_cg + lg_mass*lg_cg +
            vertical_tail_mass * vertical_tail_cg + engine_mass*engine_cg + battery_mass * battery_cg +
            payload_mass*payload_cg + system_mass * system_cg) / (MTOW)

        x_cg = x_relative_cg

        self.x_cg = x_cg
        self.x_relative_cg = x_relative_cg
        self.length = total_length

        return x_cg
    
    def calculate_MOI(self,z_cg_lg,z_cg_pylon,z_cg_vtail):
        self.z_cg_lg = z_cg_lg
        self.z_cg_pylon = z_cg_pylon
        self.z_cg_vtail = z_cg_vtail

        chord_body_section = np.linspace(self.c[0], self.c[1], 100)
        chord_y_body = np.linspace(0, self.b[1]/2, 100)

        chord_wing_section = np.linspace(self.c[1], self.c[2], 100)
        chord_y_wing = np.linspace(0, self.b[2]/2 - self.b[1]/2, 100)

        #Step2 (estimates x according to sweep)
        x_body = 0.25 * self.c[0] + np.tan(self.sweep[0]) * chord_y_body + 0.1 * chord_body_section-self.x_cg
        x_wing = (0.25 * self.c[0] + np.tan(self.sweep[1]) * self.b[1]/2) + np.tan(self.sweep[1]) * chord_y_wing + 0.1 * chord_wing_section



        self.Structure_mass=0.26*self.MTOW
        Wing_mass=0.29*self.Structure_mass
        Body_mass=0.54*self.Structure_mass
        Landing_gear_mass=0.108*self.Structure_mass
        Pylon_mass=0.035*self.Structure_mass
        Vertical_tail_mass=0.024*self.Structure_mass

        chord_y_wing=chord_y_wing+self.b[1]/2

        V_tot_structure_body=np.trapz(chord_body_section**2, chord_y_body)
        V_tot_structure_wing=np.trapz(chord_wing_section**2, chord_y_wing)
        structure_density_wing=(Wing_mass)/V_tot_structure_wing
        structure_density_body=(Body_mass)/V_tot_structure_body

        xx_function_wing=structure_density_wing*chord_wing_section**2*chord_y_wing**2
        xx_function_body=structure_density_body*chord_body_section**2*chord_y_body**2

        yy_function_wing=structure_density_wing*chord_wing_section**2*x_wing**2
        yy_function_body=structure_density_body*chord_body_section**2*x_body**2

        zz_function_wing=structure_density_wing*chord_wing_section**2*(x_wing**2+chord_y_wing**2)
        zz_function_body=structure_density_body*chord_body_section**2*(x_body**2+chord_y_body**2)

        I_xx_batt=self.battery_density*self.cg_list[1][1]*self.y_list[1]**2
        I_yy_batt=0
        I_zz_batt=0
        for i in range(len(self.cg_list[1])):
            I_yy_batt+=self.battery_density*self.cg_list[1][i]*(self.cg_list[0][i]-self.x_cg)**2
            if len(self.cg_list[1])==2:
                y_cg=[0,self.y_list[1]]
                I_zz_batt+=self.battery_density*self.cg_list[1][i]*((self.cg_list[0][i]-self.x_cg)**2+(y_cg[i])**2)
            else:
                y_cg=[0,self.y_list[1],0]
                I_zz_batt+=self.battery_density*self.cg_list[1][i]*((self.cg_list[0][i]-self.x_cg)**2+(y_cg[i])**2)


        self.I_xx=np.trapz(xx_function_wing, chord_y_wing)+np.trapz(xx_function_body, chord_y_body)+Landing_gear_mass*z_cg_lg**2+Vertical_tail_mass*z_cg_vtail**2+Pylon_mass*z_cg_pylon**2+I_xx_batt
        self.I_yy=np.trapz(yy_function_wing, chord_y_wing)+np.trapz(yy_function_body, chord_y_body)+Landing_gear_mass*((self.lg_cg-self.x_cg)**2+z_cg_lg**2)+Vertical_tail_mass*((self.vertical_tail_cg-self.x_cg)**2+self.z_cg_pylon**2)+Pylon_mass*((self.pylon_cg-self.x_cg)**2+self.z_cg_pylon**2)+I_yy_batt
        self.I_zz=np.trapz(zz_function_wing, chord_y_wing)+np.trapz(zz_function_body, chord_y_body)+Landing_gear_mass*(self.lg_cg-self.x_cg)**2+Vertical_tail_mass*(self.vertical_tail_cg-self.x_cg)**2+Pylon_mass*(self.pylon_cg-self.x_cg)**2+I_zz_batt
        self.I_xz=Landing_gear_mass*(self.lg_cg-self.x_cg)*z_cg_lg+Vertical_tail_mass*(self.vertical_tail_cg-self.x_cg)*z_cg_vtail+Pylon_mass*(self.pylon_cg-self.x_cg)*z_cg_pylon



    #def vert_stab(self):
    #    #As the moment arm is limiting on BWB without tail, the twin vert stab will be placed in
    #    # such a way that the moment arm is maximised
    #    # C_v_t =
    #    # #Single
    #    # self.asd
    #    x=0
class Tail:
    def __init__(self, plane, eta, SrS, T_engine, l_engine, d_engine, x_cg, taper_v=0.4, A_v=2.5,
                 thickness_v=0.12, def_rudder_emergency = 10, beta_max=30, Cl_alpha=2*np.pi,
                 sweep_half_v=0, V_stall = 50, density = 1.225, iteration = 1):
        self.A_v=A_v
        self.Taper_v=taper_v
        self.V_s = V_stall
        self.density = density


        self.eta = eta #How much of the span does the rudder take
        self.SrS= SrS #Part of surface used for the rudder (Roskam 2 for first iteration)
        self.thickness_v = thickness_v
        self.taper_v = taper_v

        self.df_deg = def_rudder_emergency
        self.df_rad = np.deg2rad(def_rudder_emergency)

        self.cl_alpha_theory = Cl_alpha
        self.sweep_half_v = sweep_half_v

        self.beta_max_rad = np.deg2rad(beta_max)#Make sure to add safety factor
        self.beta_engine_fail = self.beta_max_rad*0.2

        self.plane = plane #plane pobject
        self.x_cg = x_cg
        self.coords_bot=plane.coords_bot
        self.b_i=plane.b[1]

        self.T_engine = T_engine #Max thrust of one engine
        self.dy_engine = self.plane.b[1]/2 #distance of engine from centerline is on the rib section
        self.l_engine = l_engine #length of ducted fan
        self.d_engine = d_engine #diameter of ducted fan

        self.iteration = iteration

    def gradient(self, f, x, step):
        return (f(x + step) - f(x - step)) / (step * 2)

    def newtonRaphson_tail(self, f, x0, e, N, h, relax):
        #print('\n\n*** NEWTON RAPHSON METHOD IMPLEMENTATION ***')
        i = 0
        step = 1
        flag = 1
        condition = True
        while condition:
            # if g(f,x0,h) == 0.0:
            #     print('Divide by zero error!')
            #     break
            #print('x0---', x0)
            #print('value---', f(x0))
            #print('grad---', self.gradient(f, x0, h))
            x1 = x0 * relax + (x0 - f(x0) / (self.gradient(f, x0, h))) * (1 - relax)
            # print('Iteration-%d, x1 = %0.6f and f(x1) = %0.6f' % (step, x1, f(x1)))
            x0 = x1
            step = step + 1
            newvalue = f(x1)
            #print(newvalue)
            # if g(f,buildingno,x0,h)<0:
            #     x1=x1/relax

            if abs(np.max(newvalue)) < e:
                condition = False
            if step > N:
                print('\nNot Convergent.')
                flag = 2
                condition = False
            if x1 <= 0:
                condition = False
                x1=1000
            i += 1
            #print('x1---', x1)

        if flag == 1:
            print('\nRequired root is: %0.8f', x1)
            return x0, i, x1
        else:
            #print('\nNot Convergent.')
            return 1000, i, "No solution found"

    def calc_flap_span_factor(self):
        print("Considering that the flap span is {} ?".format(self.eta))
        self.Kb= 1 #float(input('What is Kb (if you dont know assume 1, Roskam 6 p.260)'))

    def calc_delta_cl(self):
        # Calculates the flap chord/ average chord form Sr/S
        self.cfc = self.SrS / self.eta  # Assumes average chord and sweep of zero

        print("Considering that the cf/c is {} and flap deflection {} deg".format(self.cfc, self.df_deg))
        self.k_theory =1 # float(input('What is k (Roskam 6 - p.228 - fig. 8.13)'))

        print("Considering that the cf/c is {} and t/c is {}".format(self.cfc, self.thickness_v))
        self.cl_delta_theory = 5 # float(input('What is cl_delta_theory (Roskam 6 - p.228 - fig. 8.14)'))

        # We assume Cl_alpha_theory = Cl_alpha

        self.delta_cl = self.df_rad * self.cl_delta_theory * self.k_theory

    def calc_CL_alpha_v(self):
        #Assumed ellipitical lift distribution
        self.CL_alpha_v = 2*np.pi/(1+2*np.pi/np.pi/self.A_v)

    def calc_rudder_effectiveness(self):
        print("Considering that the cf/c is {} and A_v is {}".format(self.cfc, self.A_v))
        print("A_v is only 2 if you choose for a body configuraton, if not, change A_v in class input to the A_v of the wingtip VS and check your findings")
        self.rudder_effectiveness = 1 # float(input('What is the flap effectiveness? ( Roskam 6 - p.261 - fig. 8.53)'))

    def calc_deltacl_rudder(self):
        self.calc_CL_alpha_v()
        self.calc_flap_span_factor()
        self.calc_delta_cl()
        self.calc_rudder_effectiveness()

        self.delta_cL_rudder = self.Kb*self.delta_cl*self.CL_alpha_v/self.cl_alpha_theory*self.rudder_effectiveness

### --- Vertical stabiliser either on the wingtips or 2 on the body --- (iteration: 0)

    def calc_xv(self,S_v): #Body
        cr = (S_v/2/self.A_v)**0.5*2/(self.taper_v+1)
        self.min_dy = np.tan(self.beta_max_rad)*cr

        MAC = 2/3 * cr*(1+self.taper_v+self.taper_v**2)/(1+self.taper_v)
        b = (self.A_v/S_v)**0.5

        self.zv = -(b)*(MAC-cr)/(cr-self.taper_v*cr)

        self.sweep_v = np.arctan((0.25*cr-0.25*MAC)/(self.zv))
        self.x_v = np.tan(self.sweep_v)*self.zv

    def calc_lv(self,S_v):#Body
        self.calc_xv(S_v)
        self.x_v_end_body = (self.coords_bot[1] - self.coords_bot[0]) / (self.b_i) * self.min_dy + self.coords_bot[0]
        # print(S_v)
        cr = (S_v/2/self.A_v)**0.5*2/(self.taper_v+1)
        self.lv =self.x_v_end_body-self.x_cg-cr+self.x_v

    def f_b(self,S_v): #Body
        self.calc_lv(S_v)  # keep in the loop
        f = - self.T_engine * self.dy_engine + 2* 0.5 * (self.CL_alpha_v * self.beta_engine_fail + self.delta_cL_rudder
                                                        ) * self.density * self.V_s ** 2 * S_v * self.lv
        return f

    def tail_sizing_1(self):
        self.calc_deltacl_rudder()

        #Body
        self.S_v_b = self.newtonRaphson_tail(self.f_b,9,0.01,10000,0.0001,0.5)[2]
        if self.S_v_b=="No solution found":
            self.S_v_b = 10000

        #Wingtips
        cr_t = self.plane.c[-1]
        self.MAC_wt = 2/3 * (cr_t+self.taper_v*cr_t+self.taper_v**2*cr_t)/(1+self.taper_v)
        self.l_wt = -0.75 * cr_t + self.coords_bot[-1] - self.x_cg #Assume ac is at quarter root chord
        self.S_v_wt = self.T_engine * self.dy_engine/(2*0.5 * (self.CL_alpha_v *self.beta_engine_fail + self.delta_cL_rudder
                                                    ) * self.density * self.V_s ** 2 * self.l_wt)



        #Compare best option
        if self.S_v_b<self.S_v_wt:
            self.S_v = 2*np.copy(self.S_v_b) #Surface for both tailplanes
            self.x_tail = self.lv+self.x_cg
            self.z_tail = self.zv
            self.A_v = self.A_v
            print("Vertical stabiliser on Body section with dy = {} ".format(self.min_dy))
            self.x_v_cg = self.x_tail

        else:
            self.S_v = 2*np.copy(self.S_v_wt) #Surface for both wingtips
            self.b = (self.S_v_wt/self.MAC_wt)
            self.x_tail = self.l_wt + self.x_cg
            self.z_tail = -(self.b)*(self.MAC_wt-cr_t)/(cr_t-self.taper_v*cr_t)
            self.A_v = self.b**2/self.S_v
            print("Vertical stabiliser on Wingtips, A is {}".format(self.A_v))

###  --- Vertical stabiliser on wingtips + one vertical stabiliser on the body --- (iteration: 1)

    def calc_xv_b_wt(self,S_v): #Body+wing tips
        cr = (S_v/self.A_v)**0.5*2/(self.taper_v+1)
        self.cr_test = cr
        #print('cr',cr)

        MAC = 2/3 * cr*(1+self.taper_v+self.taper_v**2)/(1+self.taper_v)
        self.MAC_test = MAC
        self.b_v_b = (self.A_v/S_v)**0.5

        self.z_v_b = -(self.b_v_b)*(MAC-cr)/(cr-self.taper_v*cr)

        self.sweep_v = np.arctan((0.25*cr-0.25*MAC)/(self.z_v_b))

        self.x_v_b_wt = np.tan(self.sweep_v)*self.z_v_b+0.25*cr

    def calc_lv_b_wt(self,S_v):
        self.calc_xv_b_wt(S_v)
        x_v_end_body = self.coords_bot[0]
        #print('coords', self.coords_bot)
        #print('S_v',S_v)
        cr = (S_v / self.A_v) ** 0.5 * 2 / (self.taper_v + 1)
        self.l_v_b_wt = x_v_end_body - self.x_cg - cr + self.x_v_b_wt
        #print('length',x_v_end_body)
        #print('cg',self.x_cg)
        #print('cr', cr)
        #print('x_v_b_wt',self.x_v_b_wt)

    def funct_f_b_wt(self,S_v):
        self.calc_lv_b_wt(S_v)  # keep in the loop
        f = - self.T_engine * self.dy_engine + 2* 0.5 * (self.CL_alpha_v * self.beta_engine_fail + self.delta_cL_rudder) * self.density * self.V_s ** 2* self.l_v_wt1 * self.S_v_wt1+  0.5 * (self.CL_alpha_v * self.beta_engine_fail + self.delta_cL_rudder) * self.density * self.V_s ** 2 * S_v * self.l_v_b_wt
        #print('Engine',- self.T_engine * self.dy_engine)
        #print('Wt',2* 0.5 * (self.CL_alpha_v * self.beta_engine_fail + self.delta_cL_rudder) * self.density * self.V_s ** 2* self.l_v_wt1 * self.S_v_wt1)
        #print('Body',  0.5 * (self.CL_alpha_v * self.beta_engine_fail + self.delta_cL_rudder) * self.density * self.V_s ** 2 * S_v * self.l_v_b_wt)
        #print('M_wt',2* 0.5 * (self.CL_alpha_v * self.beta_engine_fail + self.delta_cL_rudder) * self.density * self.V_s ** 2* self.l_v_wt1)
        #print('M_b',  (self.CL_alpha_v * self.beta_engine_fail + self.delta_cL_rudder)*0.5* self.density * self.V_s ** 2 * S_v * self.l_v_b_wt, self.l_v_b_wt)
        return f



    #def check_engine_wake(self):
    def tail_sizing_2(self): #Wingtips + one body configuration
        self.calc_deltacl_rudder()
        # Wingtips
        cr_t = self.plane.c[-1]
        self.MAC_wt1 = 2 / 3 * (cr_t + self.taper_v * cr_t + self.taper_v ** 2 * cr_t) / (1 + self.taper_v)
        self.S_v_wt1 = self.A_v * ((cr_t+self.taper_v*cr_t)/2)**2  # Surface for one of the wingtips
        self.b_v_wt1 = (self.S_v_wt1*self.A_v)**0.5
        self.x_tail_wt1 = -0.75 * cr_t + self.coords_bot[-1]
        self.l_v_wt1 = self.x_tail_wt1 - self.x_cg
        self.z_tail_wt1 = -(self.b_v_wt1)*(self.MAC_wt1-cr_t)/(cr_t-self.taper_v*cr_t)


        #Body
        self.S_v_b_wt1 = self.newtonRaphson_tail(self.funct_f_b_wt, 9, 0.005, 100, 0.00001, 0.75)[2]  # surface of body
        self.x_tail_b1 = self.l_v_b_wt + self.x_cg
        self.z_tail_b1 = self.z_v_b
        self.b_v_b1 = self.b_v_b
        self.x_offset_engine = -(self.dy_engine / 2 / np.tan(self.beta_max_rad) - self.coords_bot[0] + self.coords_bot[
            1] - self.l_engine/2)

        if self.x_offset_engine <= 0:
            self.x_offset_engine = 0

        self.S_v = 2 * np.copy(self.S_v_wt1)+np.copy(self.S_v_b_wt1)  # Surface for both wingtips
        #self.b = (self.S_v_wt / self.MAC_wt)
        #self.x_tail = self.l_wt + self.x_cg
        #self.z_tail = -(self.b) * (self.MAC_wt - cr_t) / (cr_t - self.taper_v * cr_t)
        #self.A_v = self.b ** 2 / self.S_v


    def Tail_positioning(self):
        if self.coords_bot[-1]>self.coords_bot[0] and self.plane.b[1]<self.min_dy:
            self.x_v_end=self.coords_bot[-1]
        else:
            self.x_v_end=(self.coords_bot[1]-self.coords_bot[0])/(self.b_i)*self.min_dy+self.coords_bot[0]

# cog
# wing Span
# sweep
# MTOW
# Tip chord
# MAC
# Dihedral
# Enginer placement



class Trim:
    def __init__(self,CL,CD,T_c,V,aoa,theta=0,q=0,beta=0,phi=0,p=0,r=0):
        self.CL=CL
        self.CD=CD
        self.T_c=T_c
        self.V=V
        self.aoa=np.deg2rad(aoa)
        self.theta=theta
        self.q=q
        self.beta=beta
        self.phi=phi
        self.p=p
        self.r=r


# taper_list = [0.267,0.34]
# test=Plane(5.8,taper_list,[38,38],[22,25])
# print(test.S)
# test.plot_plane()
# test.xflrvelues()
# # test.drawbox(0.5)
# test.drawtail(0.2)