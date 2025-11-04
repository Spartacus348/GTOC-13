#!/usr/bin/env python
from astropy import units as u
from astropy.constants import Constant
from typing import NamedTuple
import astropy.time
import boinor.bodies
import boinor.core.angles
import boinor.twobody
import csv
import math
import multiprocessing
import os.path
import pprint

file = os.path.join("ext-resources", "gtoc13_planets.csv")


class Sun:
    gravity = 139348062043.343  # km3/s2
    au = 149597870.691  # km
    day = 86400  # s
    year = 365.25  # days


start = 0
# end = 200 * Sun.year * Sun.day
# step = 8640
end = 1000
step = 1
workers = 1


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


class __Classical(NamedTuple):
    nu: float
    a: float
    e: float
    inc: float
    r: float
    raan: float
    argp: float


initialState = C7(x=-200, v=0, w=0)

if not isinstance(workers, int):
    raise ValueError("Number of workers should be an integer")

# Validate worker/step distrobution
if ((Sun.day * Sun.year * 200 / step) % workers) != 0:
    raise ValueError(
        f"Sun.day*Sun.year/step)%workers should be 0, was {(Sun.day*Sun.year/step)%workers}"
    )


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

# Altaira = boinor.bodies.Body(
#     parent=None,
#     k=Constant(
#         abbrev="GM_altaira",  # abbrev
#         name="Altaira centric gravitational constant",  # name
#         value=139348062043.343e9,  # value (given by problem) (python is nice to us :))
#         unit="m3/s2",  # unit
#         uncertainty=0,  # uncertainty
#         reference="",  # reference
#         system="si",  # system
#     ),
#     name="Altaira",
#     R=Constant(
#         "R_altaira",
#         "Altaira equatorial radius",
#         6.95700e8,  # Given
#         "m",
#         0,
#         "",
#         system="si",
#     ),
#     mass=Constant(
#         "M_altaira",
#         "Solar mass",
#         139348062043.343e9 * 6.67430e-11,  # GM (Given) * Gravitational constant
#         "kg",
#         "0",
#         "",
#         system="si",
#     ),
# )

f: list[__Classical, ...] = list()
with open(file, newline="") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        f.append(
            __Classical(
                nu=boinor.core.angles.D_to_nu(
                    boinor.core.angles.M_to_D(
                        float(row["Mean Anomaly at t=0 (deg)"]) * (math.pi / 360)
                    )
                )
                * (180 / math.pi),
                a=float(row["Longitude of the Ascending Node (deg)"]),
                e=float(row["Eccentricity ()"]),
                inc=float(row["Inclination (deg)"]),
                r=(
                    float(row["Longitude of the Ascending Node (deg)"])
                    * (1 - (float(row["Eccentricity ()"]) ** 2))
                )
                / (
                    1
                    + (
                        float(row["Eccentricity ()"])
                        * math.cos(
                            boinor.core.angles.D_to_nu(
                                boinor.core.angles.M_to_D(
                                    float(row["Mean Anomaly at t=0 (deg)"])
                                    * (math.pi / 360)
                                )
                            )
                            * (180 / math.pi)
                        )
                    )
                ),
                raan=float(row["Longitude of the Ascending Node (deg)"]),
                argp=float(row["Argument of Periapsis (deg)"]),
            )
        )
body_classics: tuple[__Classical, ...] = tuple(f)

#         # Boinor uses radians and GTOC uses degrees
#         nu = boinor.core.angles.D_to_nu(
#             boinor.core.angles.M_to_D(float(row[9]) * (math.pi / 360))
#         ) * (180 / math.pi)
#         a = float(row[4])
#         e = float(row[5])
#         inc = float(row[6])
#         r = (a * (1 - (e**2))) / (1 + (e * math.cos(nu)))
#         raan = float(row[7])
#         argp = float(row[8])
#         orbits.append(
#             boinor.twobody.Orbit.from_classical(
#                 attractor=Altaira,
#                 a=a << u.km,
#                 ecc=e << u.one,
#                 inc=inc << u.deg,
#                 raan=raan << u.deg,
#                 argp=argp << u.deg,
#                 nu=nu << u.deg,
#                 epoch=astropy.time.Time("0", format="unix"),
#                 # plane=boinor.frames.Planes.BODY_FIXED, # Custom bodies not implemented. IDK if this is a problem.
#             )
#         )
# ephem = orbits[0].to_ephem(
#     strategy=boinor.twobody.sampling.EpochsArray(
#         epochs=boinor.util.time_range(
#             astropy.time.Time("0", format="unix"),
#             end=astropy.time.Time(f"{200 * 365.25 * 86400}", format="unix"),
#         )
#     )
# )


worker_jobs = [(0, 0, 0, body_classics)] * workers
for i in range(workers):
    worker_jobs[i] = (
        int(i * start / workers),
        int(((i + 1) * end / workers)),
        step,
        body_classics,
    )
worker_jobs[-1] = (
    worker_jobs[-1][0],
    int(worker_jobs[-1][1] + step),
    step,
    body_classics,
)


def procjob(inputs: tuple[int, int, int, tuple[__Classical, ...]]) -> dict[
    int,
    tuple[__C6, ...],
]:
    db: dict[
        int,
        tuple[__C6, ...],
    ] = dict()
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

    orbits: tuple[boinor.twobody.Orbit] = list()
    for x in inputs[3]:
        orbits.append(
            boinor.twobody.Orbit.from_classical(
                attractor=Altaira,
                a=x.a << u.km,
                ecc=x.e << u.one,
                inc=x.inc << u.deg,
                raan=x.raan << u.deg,
                argp=x.argp << u.deg,
                nu=x.nu << u.deg,
                epoch=astropy.time.Time("0", format="unix"),
                # plane=boinor.frames.Planes.BODY_FIXED, # Custom bodies not implemented. IDK if this is a problem.
            )
        )
    bodycount = len(orbits)
    pprint.pp(orbits)
    for t in range(inputs[0], inputs[1], inputs[2]):
        tmplist = [__C6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)] * bodycount
        for i in range(bodycount):
            r, v = orbits[i].propagate(t << u.s).rv()
            tmplist[i] = __C6(*r, *v)
        db[t] = tuple(tmplist)
    return db


if __name__ == "__main__":
    with multiprocessing.Pool(workers) as p:
        a = p.map(procjob, worker_jobs)

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
