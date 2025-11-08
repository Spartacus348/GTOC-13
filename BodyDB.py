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
import math
import multiprocessing
import os.path
import pickle

file = os.path.join("ext-resources", "gtoc13_planets.csv")


class Sun:
    gravity = 139348062043.343  # km3/s2
    au = 149597870.691  # km
    day = 86400  # s
    year = 365.25  # days


start = 0
end = 200 * Sun.year * Sun.day
step = Sun.day
workers = 20


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
    a: float
    ecc: float
    inc: float
    raan: float
    argp: float
    nu: float


initialState = C7(x=-200, v=0, w=0)

if not isinstance(workers, int):
    raise ValueError("Number of workers should be an integer")

if ((end - start) % step) != 0:
    raise ValueError(
        f"Step size should fully encompass interval from start to end. Remainder: {((start - end) % step)}"
    )

step = int(step)


body_classics = list()
with open(file, newline="", encoding="iso-8859-1") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        body_classics.append(
            __Classical(
                nu=boinor.core.angles.D_to_nu(
                    boinor.core.angles.M_to_D(
                        float(row["Mean Anomaly at t=0 (deg)"]) * (math.pi / 360)
                    )
                )
                * (180 / math.pi),
                a=float(row["Semi-Major Axis (km)"]),
                ecc=float(row["Eccentricity ()"]),
                inc=float(row["Inclination (deg)"]),
                raan=float(row["Longitude of the Ascending Node (deg)"]),
                argp=float(row["Argument of Periapsis (deg)"]),
            )
        )
body_classics: tuple[__Classical, ...] = tuple(body_classics)

worker_jobs = [(0, 0, 0, body_classics)] * workers
for i in range(workers):
    worker_jobs[i] = (
        int(start + step * math.floor((end - start) / (workers * step)) * i),
        int(start + step * math.floor((end - start) / (workers * step)) * (i + 1)),
        step,
        body_classics,
    )
worker_jobs[-1] = (
    worker_jobs[-1][0],
    int(end + step),
    step,
    body_classics,
)


def procjob(inputs: tuple[int, int, int, tuple[__Classical, ...]]) -> dict[
    int,
    tuple[__C6, ...],
]:
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
    for x in inputs[3]:
        orbits.append(
            boinor.twobody.Orbit.from_classical(
                attractor=Altaira,
                a=x.a << u.km,
                ecc=x.ecc << u.one,
                inc=x.inc << u.deg,
                raan=x.raan << u.deg,
                argp=x.argp << u.deg,
                nu=x.nu << u.deg,
                epoch=astropy.time.Time("0", format="unix"),
                # plane=boinor.frames.Planes.BODY_FIXED, # Custom bodies not implemented. IDK if this is a problem.
            )
        )
    db: dict[
        int,
        tuple[__C6, ...],
    ] = dict()
    for t in range(inputs[0], inputs[1], inputs[2]):
        tmplist = [__C6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)] * len(orbits)
        for i in range(len(orbits)):
            r, v = orbits[i].propagate(t << u.s).rv()
            tmplist[i] = __C6(
                r[0].value, r[1].value, r[2].value, v[0].value, v[1].value, v[2].value
            )
        db[t] = tuple(tmplist)
    return db


if __name__ == "__main__":
    with multiprocessing.Pool(workers) as p:
        a = p.map(procjob, worker_jobs)
    ret: dict[int, tuple[__C6, ...]] = dict()
    for x in a:
        ret.update(x)
    with open("out.pickle", "wb") as file:
        pickle.dump(ret, file)
