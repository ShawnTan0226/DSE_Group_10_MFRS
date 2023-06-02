import numpy as np

class LandingGear:
    def __init__(self, xcg, span, sweep_QC, MTOW, chord_tip, chord_root, MAC):
        self.xcg = xcg
        self.span = span
        self.sweep = sweep_QC
        self.MTOW = MTOW
        self.chord_tip = chord_tip
        self.chord_root = chord_root
        self.MAC = MAC

    def positioning(self):
        # Assumed initial MLG position
        self.pos_x_MLG = self.xcg + 0.15 * self.MAC
        # NLG distance flowing from this
        load_percentage_NLG = 0.1   # Set percentage of MTOW NLG supports
        self.pos_x_NLG = (self.xcg - self.pos_x_MLG*(1-load_percentage_NLG))/load_percentage_NLG

        # Find wingtip distance
        self.dis_x_cg_wingtip = 0.25 * self.chord_root + np.tan(self.sweep) * self.span/2 + 0.75 * self.chord_tip - self.xcg
        # Find lm and ln
        self.dis_x_cg_NLG = self.xcg - self.pos_x_NLG
        self.dis_x_cg_MLG = self.pos_x_MLG - self.xcg


    def loads(self):
        # Set number of struts
        n_s = 2
        # Loads on the landing gears
        self.load_NLG = (self.MTOW * self.dis_x_cg_MLG) / (self.dis_x_cg_NLG + self.dis_x_cg_MLG)
        self.load_MLG = (self.MTOW * self.dis_x_cg_NLG) / n_s * (self.dis_x_cg_NLG + self.dis_x_cg_MLG)

    def clearance(self, pusher = False):
        # Longitudinal tip-over criterion
        self.ang_lon = np.radians(15)   # degrees

        # Assess the clearance in case of TE pusher engines, most outer engines then critical
        if pusher:
            #self.pos_y_engine =
            pos_x_engine = 333333333
            prop_diameter = 33333333
            # Option if wingtip is constraining
            height_MLG_tip = np.tan(self.ang_lon) * (self.dis_x_cg_wingtip - self.dis_x_cg_MLG)
            # Option if pusher engine is constraining
            height_MLG_pusher = np.tan(self.ang_lon) * (pos_x_engine - self.pos_x_MLG) + prop_diameter
            self.height_MLG = max(height_MLG_pusher,height_MLG_tip)
        else:
            self.height_MLG = np.tan(self.ang_lon) * (self.dis_x_cg_wingtip - self.dis_x_cg_MLG)

        # Lateral tip-over criterion
        self.ang_lat = np.radians(55)  # degrees
        # Perpendicular distance
        w = self.height_MLG / np.tan(self.ang_lat)
        # Angle from NLG to half track width
        beta = np.arcsin(w / (self.xcg - self.pos_x_NLG))
        self.track_width_MLG = np.tan(beta) * (self.pos_x_MLG - self.pos_x_NLG)



drone = LandingGear(4.06, 22.05, 38, 19902, 1.5506, 5.800, 3)
drone.positioning()
drone.clearance()
print(drone.height_MLG, drone.xcg, drone.pos_x_MLG, drone.pos_x_NLG)

