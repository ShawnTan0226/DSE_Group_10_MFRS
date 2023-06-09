import numpy as np
class Stab:
    def __init__(self,plane,Coefficients,MTOW):
        self.plane=plane
        self.Coeff=Coefficients
        self.V0=plane.V

        self.muc=MTOW/(Coefficients.rho*plane.S*plane.MAC)
        self.mub=MTOW/(Coefficients.rho*plane.S*plane.b_tot)

    def get_asymm_eigen(self):
        
        self.A = 16*(self.mub**3)*(self.Coeff.KX2*self.Coeff.KZ2-self.Coeff.KXZ**2)  #squared KXZ

        self.B = -4*(self.mub**2)*(2*self.Coeff.C_Y_beta*(self.Coeff.KX2*self.Coeff.KZ2-self.Coeff.KXZ**2)+self.Coeff.C_n_r*self.Coeff.KX2+self.Coeff.C_l_p*self.Coeff.KZ2+(
            self.Coeff.C_l_r+self.Coeff.C_n_p)*self.Coeff.KXZ)  # squared KXZ

        self.C = 2*self.mub*((self.Coeff.C_Y_beta*self.Coeff.C_n_r-self.Coeff.C_Y_r*self.Coeff.C_n_beta)*self.Coeff.KX2+(
            self.Coeff.C_Y_beta*self.Coeff.C_l_p-self.Coeff.C_l_beta*self.Coeff.C_Y_p)*self.Coeff.KZ2+((
            self.Coeff.C_Y_beta*self.Coeff.C_n_p-self.Coeff.C_n_beta*self.Coeff.C_Y_p)+(
            self.Coeff.C_Y_beta*self.Coeff.C_l_r-self.Coeff.C_l_beta*self.Coeff.C_Y_r))*self.Coeff.KXZ+4*self.mub*self.Coeff.C_n_beta*self.Coeff.KX2+4*self.mub*self.Coeff.C_l_beta*self.Coeff.KXZ+0.5*(self.Coeff.C_l_p*self.Coeff.C_n_r-self.Coeff.C_n_p*self.Coeff.C_l_r))

        self.D = -4*self.mub*self.Coeff.CL*(self.Coeff.C_l_beta*self.Coeff.KZ2+self.Coeff.C_n_beta*self.Coeff.KXZ)+2*self.mub*(
            self.Coeff.C_l_beta*self.Coeff.C_n_p-self.Coeff.C_n_beta*self.Coeff.C_l_p)+0.5*self.Coeff.C_Y_beta*(
            self.Coeff.C_l_r*self.Coeff.C_n_p-self.Coeff.C_n_r*self.Coeff.C_l_p)+0.5*self.Coeff.C_Y_p*(
            self.Coeff.C_l_beta*self.Coeff.C_n_r-self.Coeff.C_n_beta*self.Coeff.C_l_r)+0.5*self.Coeff.C_Y_r*(
            self.Coeff.C_l_p*self.Coeff.C_n_beta-self.Coeff.C_n_p*self.Coeff.C_l_beta) # Added minus sign to the front

        self.E = self.Coeff.CL*(self.Coeff.C_l_beta*self.Coeff.C_n_r-self.Coeff.C_n_b*self.Coeff.C_l_r)

        self.eigenvalues_asymm=np.roots([self.A, self.B, self.C, self.D, self.E])

    def get_symm_eigen(self):
        self.A = 4 * (self.muc ** 2) * self.Coeff.KY2 * (self.Coeff.C_Z_alphadot - 2 * self.muc)

        self.B = self.Coeff.C_m_alphadot * 2 * self.muc * (self.Coeff.C_Z_q + 2 * self.muc) - self.Coeff.C_m_q * 2 * self.muc * (
                    self.Coeff.C_Z_alphadot - 2 * self.muc) - 2 * self.muc * self.Coeff.KY2 * (
                             self.Coeff.C_X_u * (self.Coeff.C_Z_alphadot - 2 * self.muc) - 2 * self.muc * self.Coeff.C_Z_alpha)

        self.C = self.Coeff.C_m_alpha * 2 * self.muc * (self.Coeff.C_Z_q + 2 * self.muc) - self.Coeff.C_m_alphadot * (
                    2 * self.muc * self.Coeff.C_X_0 + self.Coeff.C_X_u * (self.Coeff.C_Z_q + 2 * self.muc)) + self.Coeff.C_m_q * (
                             self.Coeff.C_X_u * (
                                 self.Coeff.C_Z_alphadot - 2 * self.muc) - 2 * self.muc * self.Coeff.C_Z_alpha) + 2 * self.muc * (
                     self.Coeff.KY2) * (self.Coeff.C_X_alpha * self.Coeff.C_Z_u - self.Coeff.C_Z_alpha * self.Coeff.C_X_u)

        self.D = self.Coeff.C_m_u * (self.Coeff.C_X_alpha * (self.Coeff.C_Z_q + 2 * self.muc) - self.Coeff.C_Z_0 * (
                    self.Coeff.C_Z_alphadot - 2 * self.muc)) - self.Coeff.C_m_alpha * (2 * self.muc * self.Coeff.C_X_0 + self.Coeff.C_X_u * (
                    self.Coeff.C_Z_q + 2 * self.muc)) + self.Coeff.C_m_alphadot * (
                             self.Coeff.C_X_0 * self.Coeff.C_X_u - self.Coeff.C_Z_0 * self.Coeff.C_Z_u) + self.Coeff.C_m_q * (
                             self.Coeff.C_X_u * self.Coeff.C_Z_alpha - self.Coeff.C_Z_u * self.Coeff.C_X_alpha)

        self.E = -self.Coeff.C_m_u * (self.Coeff.C_X_0 * self.Coeff.C_X_alpha + self.Coeff.C_Z_0 * self.Coeff.C_Z_alpha) + self.Coeff.C_m_alpha * (
                    self.Coeff.C_X_0 * self.Coeff.C_X_u + self.Coeff.C_Z_0 * self.Coeff.C_Z_u)
    
        self.eigenvalues_symm=np.roots([self.A, self.B, self.C, self.D, self.E])