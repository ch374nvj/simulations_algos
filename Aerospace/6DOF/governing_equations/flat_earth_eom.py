import math
import numpy as np
from tools.Interpolators import interp

def flat_earth_eom(t: float, x: np.ndarray, amod: dict, vmod: dict) -> np.ndarray:
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
        amod (dict): atmospheric model stored in dict
        vmod (dict): Vehicle (aircraft) model stored in dict

    Returns:
        ndarray: dx - Time derivative of each state in x
    """
    dx = np.zeros(12)

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
    m_kg = vmod['m_kg']
    Jxz_b_kgm2 = vmod['Jxz_b_kgm2']
    Jxx_b_kgm2 = vmod['Jxx_b_kgm2']
    Jyy_b_kgm2 = vmod['Jyy_b_kgm2']
    Jzz_b_kgm2 = vmod['Jzz_b_kgm2']

    # Euler angles trignometry pre-calc
    s_phi   = math.sin(phi_rad)
    c_phi   = math.cos(phi_rad)
    s_theta = math.sin(theta_rad)
    c_theta = math.cos(theta_rad)
    t_theta = math.tan(theta_rad)
    s_psi   = math.sin(psi_rad)
    c_psi   = math.cos(psi_rad)
    
    # current altitude
    h_m = -p3_n_m
    
    # Atmosphere model
    rho_interp_kgpm3 = interp(amod['alt_m'], amod['rho_kgpm3'], h_m)
    # rho_interp_kgpm3 = 1.2
    c_interp_mps = interp(amod['alt_m'], amod['c_mps'], h_m)

    # Air data calc
    TAS_mps = math.sqrt(u_b_mps**2 + v_b_mps**2 + w_b_mps**2)
    qbar_kgpm2 = 0.5*rho_interp_kgpm3*TAS_mps**2

    if w_b_mps == 0 and u_b_mps == 0:
        w_over_u = 0
    else:
        w_over_u = w_b_mps/u_b_mps
    alpha_rad = math.atan(w_over_u)

    if v_b_mps == 0 and TAS_mps == 0:
        v_over_VT = 0
    else:
        v_over_VT = v_b_mps/TAS_mps
    beta_rad = math.asin(v_over_VT)

    s_alpha = math.sin(alpha_rad)
    c_alpha = math.cos(alpha_rad)
    s_beta  = math.sin(beta_rad)
    c_beta  = math.cos(beta_rad)

    # Gravity normal to earth tangent
    # gz_n_mps2 = 9.81
    gz_n_mps2 = interp(amod['alt_m'], amod['g_mps2'], h_m)

    # Resolve gravity in body frame
    gx_b_mps2 = -s_theta * gz_n_mps2
    gy_b_mps2 = s_phi * c_theta * gz_n_mps2
    gz_b_mps2 = c_phi * c_theta * gz_n_mps2

    # Aerodynamic forces
    drag_kgmps2 = vmod['CD_approx']*qbar_kgpm2*vmod['Aref_m2']
    side_kgmps2 = vmod['CY_approx']*qbar_kgpm2*vmod['Aref_m2']
    lift_kgmps2 = vmod['CL_approx']*qbar_kgpm2*vmod['Aref_m2']

    # External Forces
    # C_w/b
    C_w_b = np.array([[ c_alpha*c_beta,s_beta, s_alpha*c_beta],
                      [-c_alpha*s_beta,c_beta, s_alpha*s_beta],
                      [-s_alpha       ,0     , c_alpha       ]])
    
    # F^D_A
    F_D_A = np.array([drag_kgmps2, side_kgmps2, lift_kgmps2]).reshape(-1,1)
    # F^b_A
    F_b_A = -C_w_b.T @ F_D_A

    [Fx_b_kgmps2,
    Fy_b_kgmps2,
    Fz_b_kgmps2] = F_b_A.flatten().tolist()

    # External moments
    l_b_kgm2ps2 = 0
    m_b_kgm2ps2 = 0
    n_b_kgm2ps2 = 0

    # EOM
    # Translational Eqs

    # State: u_b_mps
    dx[0] = (1 / m_kg) * Fx_b_kgmps2 + gx_b_mps2 - w_b_mps * q_b_rps + v_b_mps * r_b_rps
    # State: v_b_mps
    dx[1] = (1 / m_kg) * Fy_b_kgmps2 + gy_b_mps2 - u_b_mps * r_b_rps + w_b_mps * p_b_rps
    # State: w_b_mps
    dx[2] = (1 / m_kg) * Fz_b_kgmps2 + gz_b_mps2 - v_b_mps * p_b_rps + u_b_mps * q_b_rps

    # Rotational Eqs
    Den = Jxx_b_kgm2 * Jzz_b_kgm2 - Jxz_b_kgm2**2.0
    
    # State: p_b_rps
    dx[3] = (Jxz_b_kgm2*(Jxx_b_kgm2 - Jyy_b_kgm2) * p_b_rps * q_b_rps \
        - (Jzz_b_kgm2*(Jzz_b_kgm2 - Jyy_b_kgm2) + Jxz_b_kgm2**2.0) * q_b_rps * r_b_rps \
        + Jzz_b_kgm2 * l_b_kgm2ps2 + Jxz_b_kgm2 * n_b_kgm2ps2) / Den
    
    # State: q_b_rps
    dx[4] = ((Jzz_b_kgm2 - Jxx_b_kgm2) * r_b_rps * p_b_rps - Jxz_b_kgm2 \
        * (p_b_rps**2.0 - r_b_rps**2.0) * m_b_kgm2ps2) / Jyy_b_kgm2

    # State: r_b_rps
    dx[5] = (-Jxz_b_kgm2 * (Jxx_b_kgm2 - Jyy_b_kgm2 + Jzz_b_kgm2 ) * q_b_rps * r_b_rps \
        + (Jxx_b_kgm2 * (Jxx_b_kgm2 - Jyy_b_kgm2) + Jxz_b_kgm2**2.0) * p_b_rps * q_b_rps \
        + Jxz_b_kgm2 * l_b_kgm2ps2 + Jxx_b_kgm2 * n_b_kgm2ps2) / Den

    # Kinematic Eqs
    # W: Euler angle rate mapping matrix
    W = np.array([[1, s_phi*t_theta, c_phi*t_theta ],
                  [0,      c_phi   , -s_phi        ],
                  [0, s_phi/c_theta, c_phi/c_theta ]])

    pqr = np.array([p_b_rps, q_b_rps, r_b_rps]).reshape(-1,1)

    # eul_dot: phi_dot, theta_dot, psi_dot
    eul_dot = W @ pqr

    dx[6:9] = eul_dot.flatten().tolist()

    # Position Eqs
    # C_n/b
    C_n_b = np.array([[c_theta*c_psi, (-c_phi*s_phi + s_phi*s_theta*c_psi), ( s_phi*s_psi + c_phi*s_theta*c_psi)],
                      [c_phi*s_psi  , ( c_phi*c_psi + s_phi*s_theta*s_psi), (-s_phi*c_psi + c_phi*s_theta*s_psi)],
                      [-s_theta     ,   s_phi*c_theta,                        c_phi*c_theta                     ]])

    uvw_dot = np.array(dx[:3]).reshape(-1,1)

    pos_dot = C_n_b @ uvw_dot 

    dx[9:] = pos_dot.flatten().tolist()

    return dx