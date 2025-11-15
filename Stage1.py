#!/usr/bin/env python

import csv
import math
import multiprocessing
import os
import os.path
import pickle
from sys import stderr
from typing import NamedTuple

import numpy as np
from astropy import units as u
from boinor.twobody import Orbit

from constants import Altaira, Sun

with open("BodyDB.pickle", "rb") as file:
    db: dict[int, tuple[tuple[float, float, float, float, float, float]]] = pickle.load(
        file
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


class C7(NamedTuple):
    t: float = math.nan
    x: float = math.nan
    y: float = math.nan
    z: float = math.nan
    u: float = math.nan
    v: float = math.nan
    w: float = math.nan


class Message(NamedTuple):
    past: list[C7]
    txt: str
    next: list[C7]


def job(message: Message):
    current = message.next.pop()
    r = [current[1], current[2], current[3]] << u.km
    v = [current[4], current[5], current[6]] << u.km / u.s
    t = current[0] << u.s

    orb = Orbit.from_vectors(Altaira, r, v, epoch=t)
    intersect = False
    for t in range(0, int(70 * Sun.year * Sun.day), Sun.day):
        for p in range(len(db[t])):
            probe = np.array(r)
            q = np.array((db[t][p][0:3]))
            R = radiusList[p] + 500

            inside = np.sum((probe - q) ** 2) <= R**2
            if inside:
                intersect = True
                break
        if intersect:
            break
        orb = orb.propagate(Sun.day << u.s)
        r, v = orb.rv()
    #### ACTUAL CONE OF CONCERN CALC GOES ABOVE
    # Size matters. Try to use python primatives.
    # The message gets pickled behind the scenes
    # and appended to the out_queue.
    # This doubling message reached >12k loops in
    # less than 10 min with 60 workers using 15GB ram.

    message.past.append(current)
    txt = "INTERSECT" if intersect else ""
    return Message(message.past, txt, message.next)


def worker(
    out_queue,
):
    rng = np.random.default_rng()
    while True:
        task = Message(
            past=list(),
            txt="",
            next=list(
                (
                    C7(
                        t=0,
                        x=-200 * Sun.au,
                        y=rng.uniform(low=-200 * Sun.au, high=200 * Sun.au),
                        z=rng.uniform(low=-200 * Sun.au, high=200 * Sun.au),
                        u=rng.uniform(low=1, high=200),
                        v=0,
                        w=0,
                    ),
                ),
            ),
        )
        try:
            out_queue.put(job(task), True, 10)
        except ValueError:
            print("Queue full or message too big!")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    workers = os.cpu_count()
    out_queue = multiprocessing.Queue()
    proclist: list[multiprocessing.Process] = list()
    for i in range(workers):
        proclist.append(multiprocessing.Process(target=worker, args=(out_queue,)))
        proclist[-1].start()
    counter = 0
    print("t,x,y,z,u,v,w")
    while True:
        counter += 1
        print(counter, file=stderr)
        msg = out_queue.get()
        if msg.txt != "":
            past = msg.past[0]
            print(msg, file=stderr)
            print(
                f"{float(past[0])},{float(past[1])},{float(past[2])},{float(past[3])},{float(past[4])},{float(past[5])},{float(past[6])}"
            )
