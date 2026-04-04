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

    return dx