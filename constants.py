#!/usr/bin/python
## holds the constant parameters for the Sun
from typing import NamedTuple

from astropy.constants import Constant
from boinor.bodies import Body


class Sun:
    gravity = 139348062043.343  # km3/s2
    au = 149597870.691  # km
    day = 86400  # s
    year = 365.25  # days
    end_sim_time = 200*day*year


Altaira = Body(
    parent=None,
    k=Constant(
        abbrev="GM_altaira",  # abbrev
        name="Altaira centric gravitational constant",  # name
        value=139348062043.343e9,  # value (given by problem) (python is nice to us :))
        unit="m3/s2",  # unit
        uncertainty=0,  # uncertainty
        reference="",  # reference
        system="si",  # system
    ),
    name="Altaira",
    R=Constant(
        "R_altaira",
        "Altaira equatorial radius",
        6.95700e8,  # Given
        "m",
        0,
        "",
        system="si",
    ),
    mass=Constant(
        "M_altaira",
        "Solar mass",
        139348062043.343e9 * 6.67430e-11,  # GM (Given) * Gravitational constant
        "kg",
        "0",
        "",
        system="si",
    ),
)

class C7(NamedTuple):
    t: int
    x: float
    y: float
    z: float
    u: float
    v: float
    w: float

type UnnamedTuple = tuple[int, float, float, float, float, float, float]

class Classical(NamedTuple):
    a: float
    ecc: float
    inc: float
    raan: float
    argp: float
    nu: float


def unname_tuple(state: C7) -> UnnamedTuple:
    return state.t, state.x, state.y, state.z, state.u, state.v, state.w

def name_tuple(state: UnnamedTuple) -> C7:
    return C7(
        t=state[0],
        x=state[1],
        y=state[2],
        z=state[3],
        u=state[4],
        v=state[5],
        w=state[6],
    )

if __name__ == "__main__":
    print("This module is not meant to be run.")

class SailParameters:
    r0 = 149597870.691  # km
    flux = 5.4026e-6  # N/m2
    area = 15000  # m2
    mass = 500  # kg
