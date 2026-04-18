"""
This library returns a dict `vmod` which contains the following intertial 
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

def BowlingBall():
    r_sphere_m = 0.08
    m_sphere_kg = 5
    J_sphere_kgm2 = 0.4*m_sphere_kg*r_sphere_m

    vmod = {
        "m_kg" : 1,
        "Jxz_b_kgm2" : 0,
        "Jxx_b_kgm2" : J_sphere_kgm2,
        "Jyy_b_kgm2" : J_sphere_kgm2,
        "Jzz_b_kgm2" : J_sphere_kgm2,
        "Vterm_mps"  : 10,
        "CD_approx"  : 0.01,
        "CL_approx"  : 0.0,
        "CY_approx"  : 0.0,
        "Aref_m2"    : 1   
    }

    return vmod

def Lead_50Calib():
    radius_m = 0.00635
    mass_kg = 0.012
    J_sphere_kgm2 = 0.4*mass_kg*radius_m
    Aref_m2 = math.pi * radius_m**2

    vmod = {
        "m_kg" : 1,
        "Jxz_b_kgm2" : 0,
        "Jxx_b_kgm2" : J_sphere_kgm2,
        "Jyy_b_kgm2" : J_sphere_kgm2,
        "Jzz_b_kgm2" : J_sphere_kgm2,
        "Vterm_mps"  : 10,
        "CD_approx"  : 0.47,
        "CL_approx"  : 0.0,
        "CY_approx"  : 0.0,
        "Aref_m2"    : Aref_m2   
    }

    return vmod
