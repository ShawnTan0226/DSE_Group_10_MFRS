''' This code is made to compute the battery size in the wings'''

import numpy as np
import matplotlib.pyplot as plt
import math
from Plane import Plane
from scipy.interpolate import interp1d
'''Inputs'''

# g=9.81 #Gravitational acceleration [m/s^2]
# MTOW= 19902 #Maximum Take Off Weight [kg]
# Wing_loading= 2409 #Wing Loading [N/m^2]
# S=MTOW*g/Wing_loading #Total Surface
# AR = 6 #Aspect ratio
# b = np.sqrt(S*AR) # outer wing wingspan [m]
# V_bat  = 20.5 # Battery Volume [m^3]
# V_body = 0.25*V_bat # Battery Volume [m^3]
# V_tot = V_bat + V_body #Total Volume [m^3]



# taper_outer=0.267354977
# sweep_inner=np.deg2rad(38)
# sweep_outer=np.deg2rad(38)
# b_inner=0.4*b #np.arange(0.1,0.65,0.05)*b
# b_outer=b-b_inner

# ''' Arifoil Properties '''

    
# # Read the .dat file
# file_path_inner = ".\Airfoil_dat\MH 91  14.98%.dat"
# file_path_outer = ".\Airfoil_dat\MH 91  14.98%.dat"

# def airfoilvolume(file_path):
#     # Initialize empty arrays for positive and negative values
#     positive_column1 = []
#     positive_column2 = []
#     negative_column1 = []
#     negative_column2 = []

#     # Read the .dat file
#     with open(file_path, 'r') as file:
#         for line in file:
#             # Remove leading/trailing whitespaces and split the line by spaces
#             data = line.strip().split()
#             if len(data) >= 2:  # Ensure the line has at least two columns
#                 try:
#                     value1 = float(data[0])
#                     value2 = float(data[1])
#                     if 0.15 <= value1 <= 0.55:
#                         if value2 >= 0:
#                             positive_column1.append(value1)
#                             positive_column2.append(value2)
#                         else:
#                             negative_column1.append(value1)
#                             negative_column2.append(value2)
#                 except ValueError:
#                     continue
#     # Print the positive and negative arrays
#     print("Positive Column 1:", positive_column1)
#     print("Positive Column 2:", positive_column2)
#     print("Negative Column 1:", negative_column1)
#     print("Negative Column 2:", negative_column2)


#     # Compute the surface using numpy
#     postive_surface = -np.trapz(positive_column2, positive_column1)
#     negative_surface = -np.trapz(negative_column2, negative_column1)

#     # Plot the surface
#     # plt.plot(positive_column1, positive_column2, label='Positive')
#     # plt.plot(negative_column1, negative_column2, label='Negative')
#     # plt.xlabel('Column 1 (X-axis)')
#     # plt.ylabel('Column 2 (Y-axis)')
#     # plt.title('Surface Plot')
#     # plt.grid(True)
#     # ax = plt.gca()
#     # ax.set_aspect('equal', adjustable='box')
#     # plt.draw()
#     # plt.show()

#     print(negative_surface)
#     print(postive_surface)
#     Area = (negative_surface+postive_surface)
#     return Area

# Area_inner = airfoilvolume(file_path_inner)
# Area_outer = airfoilvolume(file_path_outer)
# print('Expected volume available for bateries : {} m^2 per chord of 1m'.format(Area_outer))

# def f(x): #In here x is the inner taper ratio
#     Cri = 2 * S / ((taper_outer*x+x)*b_outer+(x+1)*b_inner)
#     y=-V_tot + 2*Area_inner * (x**2*b_inner/2+0.5*(1-x)*x*b_inner+1/6*(1-x)**2*b_inner)*Cri**2+ 2*Area_outer * (taper_outer**2*b_outer/2+0.5*(1-taper_outer)*taper_outer*b_outer+1/6*(1-taper_outer)**2*b_outer)*x**2*Cri**2
#     return y
# print(f(0.1))


# def gradient(f, x, step):
#     return (f(x + step) - f(x - step)) / (step * 2)


# def newtonRaphson(f, x0, e, N, h, relax):
#     print('\n\n*** NEWTON RAPHSON METHOD IMPLEMENTATION ***')
#     i = 0
#     step = 1
#     flag = 1
#     condition = True
#     while condition:
#         # if g(f,x0,h) == 0.0:
#         #     print('Divide by zero error!')
#         #     break
#         print('x0---', x0)
#         print('value---', f(x0))
#         print('grad---', gradient(f, x0, h))
#         x1 = x0 * relax + (x0 - f(x0) / (gradient(f, x0, h))) * (1 - relax)
#         # print('Iteration-%d, x1 = %0.6f and f(x1) = %0.6f' % (step, x1, f(x1)))
#         x0 = x1
#         step = step + 1
#         newvalue = f(x1)
#         print(newvalue)
#         # if g(f,buildingno,x0,h)<0:
#         #     x1=x1/relax

#         if abs(np.max(newvalue)) < e:
#             condition = False
#         if step > N:
#             print('\nNot Convergent.')
#             flag = 2
#             condition = False
#         i += 1
#         print('x1---', x1)

#     if flag == 1:
#         print('\nRequired root is: %0.8f' , x1)
#         return x0, i, x1
#     else:
#         print('\nNot Convergent.')
#         return 1000, i


# x1 = newtonRaphson(f,0.4,0.001,1000, 0.01, 0.5)[2]
# Cri =  2 * S / ((taper_outer*x1+x1)*b_outer+(x1+1)*b_inner)

# plt.plot(b_inner,Cri, label='Inner wing')
# plt.plot(b_inner,x1, label='Inner wing')
# # plt.show()

# print("b_outer: ",b_outer)
# print("b_inner: ",b_inner)
# print("b: ",b)
# print("taper_inner: ",x1)
# print("taper_outer: ",taper_outer)
# print("cr_inner: ",Cri)
# print("cr_outer: ",x1*Cri)
# print("ct_outer: ",x1*Cri*taper_outer)
# print("Offset inner: ",Cri*0.25-x1*Cri*0.25 + np.tan(sweep_inner)*b_inner/2)
# print("Offset outer: ",Cri*0.25-x1*Cri*taper_outer*0.25 + np.tan(sweep_outer)*b_outer/2+ np.tan(sweep_inner)*b_inner/2)



class Planform_calculation:
    def __init__(self,file_path_i,file_path_o,MTOW,Wing_loading,V_bat_pl,V_bat_prop,V_frac,b_frac,AR=6,sweep_inner=np.deg2rad(38),sweep_outer=np.deg2rad(38),taper_outer=0.267354977):
        self.g=9.81 #Gravitational acceleration [m/s^2]
        self.MTOW= MTOW #Maximum Take Off Weight [kg]
        self.Wing_loading= Wing_loading #Wing Loading [N/m^2]
        self.S=MTOW*self.g/self.Wing_loading #Total Surface
        self.AR = AR #Aspect ratio
        b = np.sqrt(self.S*AR) # outer wing wingspan [m]
        self.V_bat_pl = V_bat_pl # Battery Volume [m^3]
        self.V_bat_prop = V_bat_prop # Battery Volume [m^3]
        self.V_bat  = V_bat_pl+V_bat_prop # Battery Volume [m^3]
        self.V_body = V_frac*self.V_bat # Battery Volume [m^3]
        self.V_tot = self.V_bat + self.V_body #Total Volume [m^3]



        self.taper_outer=taper_outer
        self.sweep_inner=sweep_inner
        self.sweep_outer=sweep_outer
        self.b_inner=b_frac*b #np.arange(0.1,0.65,0.05)*b
        self.b_outer=b-self.b_inner
        self.file_path=file_path_i

        self.Area_inner=self.airfoilvolume(file_path_i)
        self.Area_outer=self.airfoilvolume(file_path_o)
        self.flying_wing=False


    def airfoilvolume(self,file_path):
        # Initialize empty arrays for positive and negative values
        positive_column1 = []
        positive_column2 = []
        negative_column1 = []
        negative_column2 = []

        positive_column11 = []
        positive_column21 = []
        negative_column11 = []
        negative_column21 = []
        self.coords_negative=np.array([[0,0]]) # Initialize empty arrays for positive and negative values
        self.coords_positive=np.array([[0,0]])

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

        eq_positive=interp1d(positive_column11[::-1],positive_column21[::-1],bounds_error=False,fill_value='extrapolate')
        eq_negative=interp1d(negative_column11,negative_column21,bounds_error=False,fill_value='extrapolate')
        thickness=eq_positive(np.arange(0,1,0.001))-eq_negative(np.arange(0,1,0.001))
        self.max_thickness = np.max(thickness)
        self.max_thickness_location = np.argmax(thickness)*0.001
        return Area
    
    def surface_area(self,in_chord,out_chord):
        # Initialize empty arrays for positive and negative values
        positive_column1 = []
        positive_column2 = []
        negative_column1 = []
        negative_column2 = []

        # Read the .dat file
        with open(self.file_path, 'r') as file:
            for line in file:
                # Remove leading/trailing whitespaces and split the line by spaces
                data = line.strip().split()
                if len(data) >= 2:  # Ensure the line has at least two columns
                    try:
                        value1 = float(data[0])
                        value2 = float(data[1])
                        if in_chord <= value1 <= out_chord:
                            if value2 >= 0:
                                positive_column1.append(value1)
                                positive_column2.append(value2)
                            else:
                                negative_column1.append(value1)
                                negative_column2.append(value2)
                    except ValueError:
                        continue
        # Print the positive and negative arrays


        # Compute the surface using numpy
        postive_surface = -np.trapz(positive_column2, positive_column1)
        negative_surface = -np.trapz(negative_column2, negative_column1)

        # Plot the surface
        # plt.plot(positive_column1, positive_column2, label='Positive')
        # plt.plot(negative_column1, negative_column2, label='Negative')
        # plt.xlabel('Column 1 (X-axis)')
        # plt.ylabel('Column 2 (Y-axis)')
        # plt.title('Surface Plot')
        # plt.grid(True)
        # ax = plt.gca()
        # ax.set_aspect('equal', adjustable='box')
        # plt.draw()
        # plt.show()

        Area = (negative_surface+postive_surface)
        return Area

    def f(self,x): #In here x is the inner taper ratio
        Cri = 2 * self.S / ((self.taper_outer*x+x)*self.b_outer+(x+1)*self.b_inner)
        y=-self.V_tot + 2*self.Area_inner * (x**2*self.b_inner/2+0.5*(1-x)*x*self.b_inner+1/6*(1-x)**2*self.b_inner)*Cri**2+ 2*self.Area_outer * (self.taper_outer**2*self.b_outer/2+0.5*(1-self.taper_outer)*self.taper_outer*self.b_outer+1/6*(1-self.taper_outer)**2*self.b_outer)*x**2*Cri**2
        return y

    def gradient(self,f, x, step):
        return (f(x + step) - f(x - step)) / (step * 2)


    def newtonRaphson(self,f, x0, e, N, h, relax):
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
            i += 1
            #print('x1---', x1)

        if flag == 1:
            #print('\nRequired root is: %0.8f' , x1)
            return x0, i, x1
        else:
            #print('\nNot Convergent.')
            return 1000, i, "No solution found"
    
    def solve_equation(self):
        x1 = self.newtonRaphson(self.f,0.4,0.001,100, 0.01, 0.5)[2]
        if x1=="No solution found":
            self.flying_wing=True
            self.Cri=2 * self.S / ((self.taper_outer+1)*(self.b_outer+self.b_inner))

            return x1
        else:
            self.Cri =  2 * self.S / ((self.taper_outer*x1+x1)*self.b_outer+(x1+1)*self.b_inner)
            self.taper_inner = x1
    
    def makeplane(self):
        self.solve_equation()
        if self.flying_wing==False:
            self.plane=Plane(self.Cri,[self.taper_inner,self.taper_outer],[np.rad2deg(self.sweep_inner),np.rad2deg(self.sweep_outer)],[self.b_inner,self.b_outer+self.b_inner],self.V_bat)
        elif self.flying_wing==True:
            self.plane=Plane(self.Cri,[self.taper_outer],[np.rad2deg(self.sweep_outer)],[self.b_inner+self.b_outer])
        return self.plane
    
    def battery_placement(self):
        V_wing=2*self.Area_outer * (self.taper_outer**2*self.b_outer/2+0.5*(1-self.taper_outer)*self.taper_outer*self.b_outer+1/6*(1-self.taper_outer)**2*self.b_outer)*self.taper_inner**2*self.Cri**2
        V_pl_extra=self.V_bat_prop-V_wing
        a=(self.plane.c[-1]-self.plane.c[1])*2/self.plane.b[-1]
        if V_pl_extra<=0:
            y=np.roots([a**2,2*a,self.plane.c[-1]**2-self.V_bat_prop/self.Area_inner])
            y_wing=np.linspace(y,self.plane.b[-1],100)
            y_wing_x=y_wing-y
            chord_wing=a*y_wing+self.plane.c[-1]
            x_wing = (0.25 * self.plane.c[0] + np.tan(self.plane.sweep[1]) * self.plane.b[1]/2) + np.tan(self.plane.sweep[1]) * y_wing_x + 0.1 * chord_wing
            self.x_cg_prop_batt=np.trapz(x_wing*chord_wing**2,y_wing)


            V_batt_body=self.V_bat_pl

            y_tot=np.linspace(0,self.plane.b[1]/2,100)
            chord=np.linspace(self.plane.c[0],self.plane.c[1],100)
            thickness=self.Area_inner*chord/0.4
            offset=np.linspace(0,self.plane.offset[1],100)
            asd=self.plane.offset[1]+0.15*self.plane.c[1]-chord*0.15-offset
            V_body_triangle=2*np.trapz(thickness*(self.plane.offset[1]+0.15*self.plane.c[1]-chord*0.15-offset),y_tot)
            if V_body_triangle>V_batt_body:
                self.option=1
                x_cg_prop_batt_body=0.5*self.plane.offset[1]
                self.x_cg_batt=(x_cg_prop_batt_body*V_batt_body+self.x_cg_prop_batt*self.V_bat_prop)/(self.V_bat)
                self.cg_list=[[self.x_cg_batt[0],x_cg_prop_batt_body],[self.V_bat_prop,V_batt_body]]
                return self.x_cg_batt[0]
            else:
                self.option=2
                V_body_rectangle=V_batt_body-V_body_triangle
                A_x=2*np.trapz(thickness,y_tot)
                x_batt_body=V_body_rectangle/A_x
                x_cg_rect=x_batt_body/2+(self.plane.offset[1]+0.15*self.plane.c[1])
                self.x_cg_batt=(x_cg_rect*V_body_rectangle+self.x_cg_prop_batt*self.V_bat_prop+0.5*self.plane.offset[1]*V_body_triangle)/(self.V_bat)
                self.cg_list=[[self.x_cg_batt[0],x_cg_rect,0.5*self.plane.offset[1]],[self.V_bat_prop,V_body_rectangle,V_body_triangle]]
                return self.x_cg_batt[0]
        else:
            chord_wing_sections = np.linspace(self.plane.c[1], self.plane.c[2], 100)
            chord_y_wing = np.linspace(0, self.plane.b[2]/2 - self.plane.b[1]/2, 100)

            x_wing = (0.25 * self.plane.c[0] + np.tan(self.plane.sweep[1]) * self.plane.b[1]/2) + np.tan(self.plane.sweep[1]) * chord_y_wing + 0.1 * chord_wing_sections

            self.x_rect_batt=0
            wing_integrand = x_wing * chord_wing_sections**2
            x_cg_prop_batt_wing = np.trapz(wing_integrand, chord_y_wing)/np.trapz(chord_wing_sections**2, chord_y_wing)

            V_batt_body=V_pl_extra+self.V_bat_pl

            y_tot=np.linspace(0,self.plane.b[1]/2,100)
            chord=np.linspace(self.plane.c[0],self.plane.c[1],100)
            thickness=self.Area_inner*chord/0.4
            offset=np.linspace(0,self.plane.offset[1],100)
            V_body_triangle=2*np.trapz(thickness*(self.plane.offset[1]+0.15*self.plane.c[1]-chord*0.15-offset),y_tot)
            if V_body_triangle>V_batt_body:
                self.option=3
                x_cg_prop_batt_body=0.66666*self.plane.offset[1]
                self.x_cg_batt=(x_cg_prop_batt_body*V_batt_body+x_cg_prop_batt_wing*V_wing)/(self.V_bat)
                self.cg_list=[[x_cg_prop_batt_body,x_cg_prop_batt_wing],[V_batt_body,V_wing]]
                return self.x_cg_batt
            else:
                self.option=4
                V_body_rectangle=V_batt_body-V_body_triangle
                A_x=2*np.trapz(thickness,y_tot)
                x_batt_body=V_body_rectangle/A_x
                self.x_rect_batt=x_batt_body
                x_cg_rect=x_batt_body/2+(self.plane.offset[1]+0.15*self.plane.c[1])
                # plt.scatter(0,0.66666*self.plane.offset[1])
                # plt.scatter(0,x_cg_rect)
                x_cg_triangle=0.66666*self.plane.offset[1]
                self.x_cg_batt=(x_cg_rect*V_body_rectangle*0.9+x_cg_prop_batt_wing*V_wing+x_cg_triangle*V_body_triangle)/(self.V_bat)
                self.cg_list=[[x_cg_rect,x_cg_prop_batt_wing,x_cg_triangle],[V_body_rectangle*0.9,V_wing,V_body_triangle]]
                return self.x_cg_batt
            





            

        
        
