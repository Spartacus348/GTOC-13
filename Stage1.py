#!/usr/bin/env python

import csv
import multiprocessing
import os
import os.path
import pickle
from collections.abc import Callable
from sys import stderr

import numpy as np
from astropy import units as u
from boinor.twobody import Orbit

import initial_conditions
import runner
from constants import C7, Altaira, Sun

with open("BodyDB.pickle", "rb") as file:
    db: dict[int, tuple[tuple[float, float, float, float, float, float], ...]] = (
        pickle.load(file)
    )


radiusList: list[float] = list()
with open(
    os.path.join("ext-resources", "gtoc13_planets.csv"),
    newline="",
    encoding="iso-8859-1",
) as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        radiusList.append(float(row["Radius (km)"]))


def effect(msg: runner.Message):
    if msg.txt not in ("", "ERR"):
        past = msg.past[0][0]
        print(msg, file=stderr)
        print(
            f"{float(past[0])},{float(past[1])},{float(past[2])},{float(past[3])},{float(past[4])},{float(past[5])},{float(past[6])},{int(msg.past[0][1])}"
        )
    return get_job(job)


def job(message: runner.Message):
    current = message.next.pop()[0]
    r = [current[1], current[2], current[3]] << u.km
    v = [current[4], current[5], current[6]] << u.km / u.s
    t = current[0] << u.s
    txt = ""
    orb = Orbit.from_vectors(Altaira, r, v, epoch=t)
    intersect = False
    for t in range(0, int(70 * Sun.year * Sun.day), Sun.day):
        for p in range(len(db[t])):
            probe = np.array(r)
            q = np.array((db[t][p][0:3]))
            R = radiusList[p] + 100

            inside = np.sum((probe - q) ** 2) <= R**2
            if inside:
                intersect = p
                break
        if intersect:
            break
        try:
            orb = orb.propagate(Sun.day << u.s)
        except AssertionError:
            txt = "ERR"
            break
        r, v = orb.rv()

    txt = "INTERSECT" if intersect else txt
    message.past.append((current, intersect if intersect else None))
    return runner.Message(
        past=message.past,
        txt=txt,
        next=runner.empty_nexts(),
        func=effect,
    )


def get_job(job: Callable[[runner.Message], runner.Message]) -> runner.Message:
    rng = np.random.default_rng()
    t_range = [10 * Sun.year, 70 * Sun.year]
    t_target = int(np.floor(rng.uniform(low=t_range[0], high=t_range[1])) * Sun.day)
    t_start = rng.uniform(
        low=0, high=t_target - t_range[0]
    )  # t_start is no nearer than 10 years before target
    planet_id = np.floor(rng.uniform(low=0, high=10))
    entry = db[t_target][int(planet_id)]
    state_start = initial_conditions.from_target(
        target=np.array([float(t_target), entry[0], entry[1], entry[2]]),
        initial_t=t_start,
    )

    return runner.Message(
        past=list(),
        txt="",
        next=[
            (
                C7(
                    t=state_start[0],
                    x=state_start[1],
                    y=state_start[2],
                    z=state_start[3],
                    u=state_start[4],
                    v=0,
                    w=0,
                ),
                None,
            ),
        ],
        func=job,
    )


if __name__ == "__main__":
    multiprocessing.freeze_support()
    in_queue = multiprocessing.SimpleQueue()
    out_queue = multiprocessing.Queue()
    workers = os.cpu_count() - 1 if os.cpu_count() is not None else 4
    for i in range(workers):
        in_queue.put(get_job(job))
    print("t,x,y,z,u,v,w,planetID")
    runner.runner(in_queue, out_queue)
