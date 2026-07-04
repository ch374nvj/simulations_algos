"""
This library returns a class with an attribute dict `vmod` which contains the following intertial 
properties of the called object.
m_kg        = vmod['m_kg']
Jxz_b_kgm2  = vmod['Jxz_b_kgm2']
Jxx_b_kgm2  = vmod['Jxx_b_kgm2']
Jyy_b_kgm2  = vmod['Jyy_b_kgm2']
Jzz_b_kgm2  = vmod['Jzz_b_kgm2']
Terminal_velocity_mps = vmod['Vterm_mps']
CD_approx = vmod['CD_approx']
CL_approx = vmod['CL_approx']
CY_approx = vmod['CY_approx']
Aref_m2 = vmod['Aref_m2]
"""

import math

class BowlingBall():
    def __init__(self) -> None:
        self.r_sphere_m = 0.08
        self.m_sphere_kg = 5
        self.J_sphere_kgm2 = 0.4*self.m_sphere_kg*self.r_sphere_m
        self.Aref_m2 = math.pi * self.r_sphere_m**2
        
        self.CD_approx = 0.01
        self.CL_approx = 0.0
        self.CY_approx = 0.0


        self.Vterm_mps = math.sqrt((2*self.m_sphere_kg*9.81)/(1.20*self.CD_approx*self.Aref_m2))

        # For legacy
        self.vmod = {
            "m_kg" : 1,
            "Jxz_b_kgm2" : 0,
            "Jxx_b_kgm2" : self.J_sphere_kgm2,
            "Jyy_b_kgm2" : self.J_sphere_kgm2,
            "Jzz_b_kgm2" : self.J_sphere_kgm2,
            "Vterm_mps"  : self.Vterm_mps,
            "CD_approx"  : self.CD_approx,
            "CL_approx"  : 0.0,
            "CY_approx"  : 0.0,
            "Aref_m2"    : 1   
        }

class Lead_50Calib():
    def __init__(self) -> None:
    
        self.radius_m = 0.00635
        self.mass_kg = 0.012
        self.J_sphere_kgm2 = 0.4*self.mass_kg*self.radius_m
        self.Aref_m2 = math.pi * self.radius_m**2
        self.CD_approx = 0.47
        self.CL_approx = 0.0
        self.CY_approx = 0.0

        # Vterm = sqrt((2*m*G)/rho*C_D*Aref)
        self.Vterm_mps = math.sqrt((2*self.mass_kg*9.81)/(1.20*self.CD_approx*self.Aref_m2))

        # Lagacy
        self.vmod = {
            "m_kg" : self.mass_kg,
            "Jxz_b_kgm2" : 0,
            "Jxx_b_kgm2" : self.J_sphere_kgm2,
            "Jyy_b_kgm2" : self.J_sphere_kgm2,
            "Jzz_b_kgm2" : self.J_sphere_kgm2,
            "Vterm_mps"  : self.Vterm_mps,
            "CD_approx"  : 0.47,
            "CL_approx"  : 0.0,
            "CY_approx"  : 0.0,
            "Aref_m2"    : self.Aref_m2   
        }
