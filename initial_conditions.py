##########################################
# code for returning the initial orbits
#
# Gabe
##########################################

import numpy as np

from constants import *

state_vector_dtype = np.dtype([('t',np.int32),('x',np.float64),('y',np.float64),('z',np.float64),('vx',np.float64),('vy',np.float64),('vz',np.float64)])
position_dtype = np.dtype([('t',np.int32),('x',np.float64),('y',np.float64),('z',np.float64)])

def from_target(target: np.ndarray(dtype=position_dtype)) -> np.ndarray(dtype=state_vector_dtype):
    """
    In progress; no-op rn
    :param target: position in time to target
    :return: the initial state vector for the orbit
    """
    # given an input position in time, returns a state vector that will reach the target given the following restrictions:
    # [_, -200 [AU], _, _, _, 0 [km/s], 0 [km/s]]
    # first value to set is H, the angular momentum vector. The position vector only determines a 2D plane of values
    # which can then be refined by the restriction to have [u, 0, 0] at [-200 AU, y, z]
    # [Hx, Hy, Hz] = [x0,y0,z0] X [vx0, vy0, vz0] = [K, y0, z0] X [vx0, 0, 0] = [0, z0*vx0, -y0*vx0]
    # [Hx, Hy, Hz] = [x,y,z] X [u1, v1, w1] = [y*w1-z*v1, z*u1-x*w1, x*v1-y*u1];
    # y*w1-z*v1 = 0, z*u1-x*w1 = z0*vx0, y*u1-x*v1 = y0*vx0
    # y/z = v1/w1,                     , y*(u1-w1*x/z) = y0*vx0
    #*y/z = v1/w1 = y0/z0
    # transform into a perifocal plane, p1 x p2. periapse is [-1, 0]. The transformation from XYZ to P1P2P3 is X-axis,
    #   then around the resulting Z axis
    #   so introduce a variable d = hyp(y,z), and the next work is going to be in the XDH system
    # p = vXh/mu - r^ //note: scaled by the eccentricity)
    #   = (r*||v||^2 - v(r*v) - mu*r/||r||)/mu // use the bac-cab rule to get rid of H-vector
    #   = ([-200, d, 0]*u - [u,0,0]*200u - mu * [200/hyp(200,d), d/hyp(200,d), 0])/mu
    #   = [-200u - 200u^2 - 200*mu/hyp(200,d), du-d/hyp(200,d), 0]/mu
    #   = [-200(u^2 + u + mu/hyp(200,d))/mu, d(u-1/hyp(200,d))/mu, 0]
    # now we have p and h, derive q
    # q = hXp = [0,0,1]x[-200(u^2 + u + mu/hyp(200,d))/mu, d(u-1/hyp(200,d))/mu, 0]
    #   = [-d(u-1/hyp(200,d))/mu, -200(u^2 + u + mu/hyp(200,d))/mu, 0]
    # The goal is apparent at this point; find some equation relating u and a scaling factor for d to find arcs to [x1,y1,z1]
    # let s be that scaling factor
    # h(scalar) = su
    # e = ||p|| = sqrt([derive yourself if you actually care, there's no need to spell it out])
    # r = hh/mu * 1/(1+e*cos(theta)
    # 1+e*cos(theta) = hh/(r*mu) => theta = acos(e*ssuu/(hyp(200,s)*mu)-1))
    pass
