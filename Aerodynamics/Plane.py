import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import interp1d

class Plane:
    def __init__(self,Cri,taper,sweep,b,h=5000,V=110,airfoil=".\Airfoil_dat\MH 91  14.98%.dat"):
        #Plane object has n sections
        self.c=np.array([Cri]) #array of chord [m] form: [Middle c,c1,c2,c3,...,cn]
        self.taper = np.array(taper) #array of taper ratio form: [taper1,taper2,taper3,...,tapern]
        self.sweep = np.array(np.deg2rad(sweep)) #array of sweep angle [rad] form: [sweep1,sweep2,sweep3,...,sweepn]
        self.b = np.concatenate(([0],b)) #array of span [m] form: [0,b1,b2,b3,...,bn]
        self.S_list=np.array([]) #array of surface area of each section [m^2] form: [S1,S2,S3,...,Sn]
        
        self.coords=np.array([]) 
        self.V=V
        self.h=h

        self.MAC_list=np.array([])
        self.x_list=np.array([])
        self.y_list = np.array([])

        self.draw()
        self.MAC_aircraft()
        self.define_airfoil(airfoil)
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
            print('sweep',self.sweep[count])
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


        self.bfull=np.concatenate((self.b,self.b[::-1],negative))

    def plot_plane(self):
        plt.plot(self.bfull,self.coords,color='black')
        plt.fill(self.bfull,self.coords, color='gray', alpha=0.5)
        plt.gca().invert_yaxis()
        plt.show()

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
        y = 2/3 * (cr + ct - (cr*ct)/(cr+ct))
        #offset from that parts root on quarter chord
        off_y = -(b/2)*(y-cr)/(cr-ct) #b is the span of that part of the wing
        off_x = off_y*np.tan(sweep)
        return y, off_x, off_y

    def listgenerator(self):
        for i in range(len(self.taper)):
            part= self.MAC_part(self.c[i],self.c[i+1],self.sweep[i],self.b[i+1]-self.b[i])
            self.MAC_list=np.concatenate((self.MAC_list,[part[0]]))
            self.x_list=np.concatenate((self.x_list,[part[1]]+self.offset[i]+0.25*self.c[i]))
            self.y_list = np.concatenate((self.y_list, [part[2]]+self.b[i]))

    def MAC_aircraft(self):#Make sure to use numpy array
        self.listgenerator()
        self.MAC = np.sum(self.MAC_list*self.S_list)/np.sum(self.S_list)
        self.x_quarter = np.sum(self.x_list*self.S_list)/np.sum(self.S_list)
        self.y_quarter = np.sum(self.y_list*self.S_list)/np.sum(self.S_list)

    
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


    def COG(self, pylon_cg, lg_cg, vertical_tail_cg, engine_mass, engine_cg, battery_mass, battery_cg, payload_mass, payload_cg, system_mass, system_cg, MTOW=19900):
        """Steps to determine cog wing structure :
        1. Partition wing into sections with each section having the corresponding c
        2. Add weights to body and wing section
        3. Set up equation (include the offset (due to sweep) already for the x of each section)"""
        #Step1
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
        wing_cg = (2*body_integration + wing_integration) / (np.trapz(chord_wing_sections**2, chord_y_wing) + 2 * np.trapz(chord_body_section**2, chord_y_body))
        """"----------------wing structure section---------------------"""


        #----------------total length of drone according to sweep----------
        total_length = 0.25 * self.c[0] + np.tan(self.sweep[0]) * 0.5*self.b[1] + np.tan(self.sweep[1]) * (0.5*self.b[2] - 0.5*self.b[1]) + 0.75 * self.c[2]

        #-----------COG DETERMINATION with all other subsystems-----
        total_structural_mass = 0.26 * MTOW


        #Defining each subsystem
        wing_mf = 0.290
        # wing_cg = 0.525
        wing_cg_relative = wing_cg / total_length
        wing_mass = wing_mf * total_structural_mass

        body_mf = 0.544
        body_cg = 0.375
        body_mass = body_mf * total_structural_mass

        pylon_mf = 0.035
        # pylon_cg = 0.85 #VARIABLE
        pylon_mass = pylon_mf * total_structural_mass #May be Variable

        lg_mf = 0.108
        # lg_cg = 0.65 #VARIABLE
        lg_mass = lg_mf * total_structural_mass #may be VARIABLE

        vertical_tail_mf = 0.024
        # vertical_tail_cg = 0.875 #VARIABLE
        vertical_tail_mass = vertical_tail_mf * total_structural_mass #may be variable

        # engine_mass = 1244.168
        # engine_cg = 0.85

        # battery_mass = 10600.13 #battery for drone
        # battery_cg = 0.2814

        # payload_mass = 2903.2 #battery for ac
        # payload_cg = 0.2814

        # system_mass = ??
        # system_cg = 0.2814

        #COG EQUATION
        total_relative_cg = (wing_mass * wing_cg_relative + body_mass * body_cg + pylon_mass * pylon_cg + lg_mass*lg_cg +
            vertical_tail_mass * vertical_tail_cg + engine_mass*engine_cg + battery_mass * battery_cg +
            payload_mass*payload_cg + system_mass * system_cg) / (MTOW)

        total_cg = total_relative_cg * total_length
        print(f"cg_rel: {total_relative_cg}")
        print(f"cg_abs: {total_cg}")

        return total_cg






# taper_list = [0.267,0.34]
# test=Plane(5.8,taper_list,[38,38],[22,25])
# print(test.S)
# test.plot_plane()
# test.xflrvelues()
# # test.drawbox(0.5)
# test.drawtail(0.2)