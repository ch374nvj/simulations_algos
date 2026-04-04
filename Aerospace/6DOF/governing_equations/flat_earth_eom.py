import math
import numpy as np

def flat_earth_eom(t: float, x: np.ndarray, amod: dict) -> np.ndarray:
    """flat_earth_eom.py is a function containing essential elements of a 6DOF 
    simulation. The purpose of this function is to allow the numerical approximation
    of solutions of the governing equations for an aircraft.

    Naming convention: <var name>_<coordinate s/m if applicable>_<units>

    Args:
        t (float): time [s]
        x (np.ndarray): state vector at time t [various units]
            x[0] : u_b_mps,   axial velocity of CM wrt inertial CS resolved in aircraft bod fixed CS
            x[1] : v_b_mps,   lateral velocity  ^^
            x[2] : w_b_mps,   vertical velocity ^^
            x[3] : p_b_rps,   roll angular velocity (rate) of body fixed CS resolved in inertial CS
            x[4] : q_b_rps,   pitch angular velocity (rate) ^^
            x[5] : r_b_rps,   yaw angular velocity (rate)   ^^
            x[6] : phi_rad,   roll angle
            x[7] : theta_rad, pitch angle
            x[8] : psi_rad,   yaw angle
            x[9] : p1_n_m,    x-axis position of aircraft resolved in NED CS
            x[10]: p2_n_m,    y-axis ^^
            x[11]: p3_n_m,    z-axis ^^
        amod (dict): Aircraft model stored in dict

    Returns:
        ndarray: dx - Time derivative of each state in x
    """
    dx = np.array((12,1))

    # Assign current state vols to variables 
    u_b_mps   = x[0]  
    v_b_mps   = x[1]
    w_b_mps   = x[2]
    p_b_rps   = x[3]
    q_b_rps   = x[4]
    r_b_rps   = x[5]
    phi_rad   = x[6]
    theta_rad = x[7]
    psi_rad   = x[8]
    p1_n_m    = x[9]
    p2_n_m    = x[10]
    p3_n_m    = x[11]

    # Get mass and MI 
    m_kg = amod['m_kg']
    Jxz_b_kgm2 = amod['Jxz_b_kgm2']
    Jxx_b_kgm2 = amod['Jxx_b_kgm2']
    Jyy_b_kgm2 = amod['Jyy_b_kgm2']
    Jzz_b_kgm2 = amod['Jzz_b_kgm2']

    # Air data calc

    # Atmosphere model

    # Gravity normal to earth tangent
    gz_n_mps2 = 9.81

    # Resolve gravity in body frame
    gx_b_mps2 = -math.sin(theta_rad) * gz_n_mps2
    gy_b_mps2 = math.sin(phi_rad) * math.cos(theta_rad) * gz_n_mps2
    gz_b_mps2 = math.cos(phi_rad) * math.cos(theta_rad) * gz_n_mps2

    # External Forces
    Fx_b_kgmps2 = 0
    Fy_b_kgmps2 = 0
    Fz_b_kgmps2 = 0
    

    # External moments
    l_b_kgm2ps2 = 0
    m_b_kgm2ps2 = 0
    n_b_kgm2ps2 = 0

    # EOM
    # Translational Eqs

    # State: u_b_mps
    dx[0] = (1 / m_kg) * Fx_b_kgmps2 - gx_b_mps2 + w_b_mps * q_b_rps + v_b_mps * r_b_rps
    # State: v_b_mps
    dx[1] = (1 / m_kg) * Fy_b_kgmps2 - gy_b_mps2 + u_b_mps * r_b_rps + w_b_mps * p_b_rps
    # State: w_b_mps
    dx[2] = (1 / m_kg) * Fz_b_kgmps2 - gz_b_mps2 + v_b_mps * p_b_rps + u_b_mps * q_b_rps

    # Rotational Eqs
    Den = Jxx_b_kgm2 * Jzz_b_kgm2 - Jxz_b_kgm2^2
    
    # State: p_b_rps
    dx[3] = (Jxz_b_kgm2(Jxx_b_kgm2 - Jyy_b_kgm2) * p_b_rps * q_b_rps \
        - (Jzz_b_kgm2(Jzz_b_kgm2 - Jyy_b_kgm2) + Jxz_b_kgm2^2) * q_b_rps * r_b_rps \
        + Jzz_b_kgm2 * l_b_kgm2ps2 + Jxz_b_kgm2 * n_b_kgm2ps2) / Den
    
    # State: q_b_rps
    dx[4] = ((Jzz_b_kgm2 - Jxx_b_kgm2) * r_b_rps * p_b_rps - Jxz_b_kgm2 \
        * (p_b_rps^2 - r_b_rps^2) * m_b_kgm2ps2) / Jyy_b_kgm2

    # State: r_b_rps
    dx[5] = (-Jxz_b_kgm2 * (Jxx_b_kgm2 - Jyy_b_kgm2 + Jzz_b_kgm2 ) * q_b_rps * r_b_rps \
        + (Jxx_b_kgm2 * (Jxx_b_kgm2 - Jyy_b_kgm2) + Jxz_b_kgm2^2) * p_b_rps * q_b_rps \
        + Jxz_b_kgm2 * l_b_kgm2ps2 + Jxx_b_kgm2 * n_b_kgm2ps2) / Den

    # Kinematic Eqs

    # Position Eqs

    return dx