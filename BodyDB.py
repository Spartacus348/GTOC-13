#!/usr/bin/env python
import os.path
from typing import NamedTuple
import math
import pprint
import random
import pickle
from sys import getsizeof
import csv
import boinor.core.angles
import boinor.twobody
import boinor.bodies
import boinor.frames
import boinor.twobody.sampling
import boinor.util
from astropy import units as u
from astropy.constants import Constant
import astropy.time
import astropy.coordinates
import multiprocessing

import timeit

file = os.path.join("ext-resources", "gtoc13_planets.csv")

step = 864
workers = 4


class Sun:
    gravity = 139348062043.343  # km3/s2
    au = 149597870.691  # km
    day = 86400  # s
    year = 365.25  # days


class SailParameters:
    r0 = 149597870.691  # km
    flux = 5.4026e-6  # N/m2
    area = 15000  # m2
    mass = 500  # kg


class C7(NamedTuple):
    t: float = math.nan
    x: float = math.nan
    y: float = math.nan
    z: float = math.nan
    u: float = math.nan
    v: float = math.nan
    w: float = math.nan


class __C6(NamedTuple):
    x: float = math.nan
    y: float = math.nan
    z: float = math.nan
    u: float = math.nan
    v: float = math.nan
    w: float = math.nan


initialState = C7(x=-200, v=0, w=0)

if not isinstance(workers, int):
    raise ValueError("Number of workers should be an integer")

# Validate worker/step distrobution
if ((Sun.day * Sun.year * 200 / step) % workers) != 0:
    raise ValueError(
        f"Sun.day*Sun.year/step)%workers should be 0, was {(Sun.day*Sun.year/step)%workers}"
    )

worker_ranges = [(0, 0, 0)] * workers
for i in range(workers):
    worker_ranges[i] = (
        int(i * Sun.day * Sun.year * 200 / workers),
        int(((i + 1) * Sun.day * Sun.year * 200 / workers)),
        step,
    )
worker_ranges[-1] = (worker_ranges[-1][0], int(worker_ranges[-1][1] + step), step)
# Python range is end is non-inclusive, so we bump the end to be included
db: dict[
    float,
    tuple[
        __C6,
        __C6,
        __C6,
        __C6,
        __C6,
        __C6,
        __C6,
        __C6,
        __C6,
        __C6,
        __C6,
        __C6,
        __C6,
        __C6,
    ],
] = dict()
# time is the key
# A tuple with __C6 of the planets follows

Altaira = boinor.bodies.Body(
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

orbits: list[boinor.twobody.Orbit] = list()
with open(file, newline="") as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        try:
            int(row[0])
        except:
            continue  # The first row is headers

        # Boinor uses radians and GTOC uses degrees
        nu = boinor.core.angles.D_to_nu(
            boinor.core.angles.M_to_D(float(row[9]) * (math.pi / 360))
        ) * (180 / math.pi)
        a = float(row[4])
        e = float(row[5])
        inc = float(row[6])
        r = (a * (1 - (e**2))) / (1 + (e * math.cos(nu)))
        raan = float(row[7])
        argp = float(row[8])
        orbits.append(
            boinor.twobody.Orbit.from_classical(
                attractor=Altaira,
                a=a << u.km,
                ecc=e << u.one,
                inc=inc << u.deg,
                raan=raan << u.deg,
                argp=argp << u.deg,
                nu=nu << u.deg,
                epoch=astropy.time.Time("0", format="unix"),
                # plane=boinor.frames.Planes.BODY_FIXED, # Custom bodies not implemented. IDK if this is a problem.
            )
        )
# ephem = orbits[0].to_ephem(
#     strategy=boinor.twobody.sampling.EpochsArray(
#         epochs=boinor.util.time_range(
#             astropy.time.Time("0", format="unix"),
#             end=astropy.time.Time(f"{200 * 365.25 * 86400}", format="unix"),
#         )
#     )
# )


def procjob(
    conditions: tuple[int, int, int], orbits: list[boinor.twobody.Orbit]
) -> None:
    db: dict[
        int,
        tuple[
            __C6,
            __C6,
            __C6,
            __C6,
            __C6,
            __C6,
            __C6,
            __C6,
            __C6,
            __C6,
            __C6,
            __C6,
            __C6,
            __C6,
        ],
    ] = dict()
    bodycount = len(orbits)
    for t in range(*conditions):  # SPLAT!
        tmplist = [(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)] * bodycount
        for i in range(bodycount):
            r, v = orbits[i].propagate(t << u.s).rv()
            tmplist[i] = (*r, *v)
        db[t] = tuple(tmplist)


print(timeit.timeit("procjob((0, 10, 1), orbits)", globals=globals(), number=100))

# Demo data to test object size. Obj is 2621528 bytes.
#
# for t in range(0, int(200 * 365.25 * 86400) + 86400, 864):
#     if t % (86400 * 365.25) == 0:
#         print((t / (86400 * 365.25)) / 2)
#     db[t] = (
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#         __C6(
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#             random.uniform(-200, 200),
#         ),
#     )


# def get_C7_from_orbit(t: float, orbitparms: None) -> C7:
#     pass


# print(getsizeof(db))

# with open("out.pickle", "wb") as file:
#     pickle.dump(db, file)
