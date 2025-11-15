#!/usr/bin/env python
import csv
import math
import multiprocessing
import os.path
import pickle

import astropy.time
import boinor.bodies
import boinor.core.angles
import boinor.twobody
import numpy as np
from astropy import units as u
from astropy.constants import Constant

from constants import Classical, Sun, UnnamedTuple

print("Done importing!")


def procjob(inputs) -> None:
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
    for t in range(inputs[0], inputs[1], inputs[2]):
        tmplist = [
            (
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            )
        ] * len(orbits)
        for i in range(len(orbits)):
            r, v = orbits[i].propagate(t << u.s).rv()
            tmplist[i] = (
                float(r[0].value),
                float(r[1].value),
                float(r[2].value),
                float(v[0].value),
                float(v[1].value),
                float(v[2].value),
            )
        inputs[4].put((t, tuple(tmplist)))
    return None


if __name__ == "__main__":
    start = 0
    end = 200 * Sun.year * Sun.day
    step = Sun.day
    workers = os.cpu_count() - 1 if os.cpu_count() is not None else 4

    if not isinstance(workers, int):
        raise ValueError("Number of workers should be an integer")

    if ((end - start) % step) != 0:
        raise ValueError(
            f"Step size should fully encompass interval from start to end. Remainder: {((start - end) % step)}"
        )

    step = int(step)

    file = os.path.join("ext-resources", "gtoc13_planets.csv")
    body_classics = list()
    with open(file, newline="", encoding="iso-8859-1") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ecc = float(row["Eccentricity ()"])
            mea = np.deg2rad(float(row["Mean Anomaly at t=0 (deg)"]))
            if ecc > 1:
                nu = boinor.core.angles.F_to_nu(
                    boinor.core.angles.M_to_F(mea, ecc), ecc
                )
            elif ecc < 1:
                nu = boinor.core.angles.E_to_nu(
                    boinor.core.angles.M_to_E(mea, ecc), ecc
                )
            else:
                nu = boinor.core.angles.D_to_nu(boinor.core.angles.M_to_D(mea))

            body_classics.append(
                Classical(
                    nu=np.rad2deg(nu),
                    a=float(row["Semi-Major Axis (km)"]),
                    ecc=float(row["Eccentricity ()"]),
                    inc=float(row["Inclination (deg)"]),
                    raan=float(row["Longitude of the Ascending Node (deg)"]),
                    argp=float(row["Argument of Periapsis (deg)"]),
                )
            )
    body_classics: tuple[Classical, ...] = tuple(body_classics)

    queue = multiprocessing.SimpleQueue()
    worker_jobs = [(0, 0, 0, body_classics, queue)] * workers
    for i in range(workers):
        worker_jobs[i] = (
            int(start + step * math.floor((end - start) / (workers * step)) * i),
            int(start + step * math.floor((end - start) / (workers * step)) * (i + 1)),
            step,
            body_classics,
            queue,
        )
    worker_jobs[-1] = (worker_jobs[-1][0], int(end + step), step, body_classics, queue)
    multiprocessing.freeze_support()
    proclist: list[multiprocessing.Process] = list()
    for i in worker_jobs:
        proclist.append(multiprocessing.Process(target=procjob, args=(i,)))
    for i in proclist:
        i.start()
    ret: dict[int, tuple[UnnamedTuple, ...]] = dict()
    for i in range(int((end - start) / step) + 1):
        item = queue.get()
        ret[item[0]] = item[1]
    with open("BodyDB.pickle", "wb") as file:
        pickle.dump(ret, file)
