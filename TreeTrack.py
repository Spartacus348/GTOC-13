#!/usr/bin/env python

import multiprocessing
import math
from typing import NamedTuple
from astropy import units as u
from boinor.twobody import Orbit
from astropy.constants import Constant
import boinor.bodies
import pprint
import time
import os


class C7(NamedTuple):
    t: int
    x: float
    y: float
    z: float
    u: float
    v: float
    w: float

type UnnamedTuple = tuple[int, float, float, float, float, float, float]

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

class Message(NamedTuple):
    past: list[C7]
    txt: str
    next: list[C7]


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


def job(message: Message):
    current = message.next.pop()
    r = [current[1], current[2], current[3]] << u.km
    v = [current[4], current[5], current[6]] << u.km / u.s
    t = current[0] << u.s

    ### This will change from Altaira if we're
    ### in a different sphere of influence.
    orb = Orbit.from_vectors(Altaira, r, v, epoch=t)

    #### ACTUAL CONE OF CONCERN CALC GOES BELOW
    ar, av = orb.propagate(60 * 60 * 24 << u.s).rv()
    message.next.append(
        (
            int(t.value + 60 * 60 * 24),
            float(ar[0].value),
            float(ar[1].value),
            float(ar[2].value),
            float(av[0].value + 0.1),
            float(av[1].value),
            float(av[2].value),
        )
    )
    message.next.append(
        (
            int(t.value + 60 * 60 * 24),
            float(ar[0].value),
            float(ar[1].value),
            float(ar[2].value),
            float(av[0].value),
            float(av[1].value),
            float(av[2].value),
        )
    )
    #### ACTUAL CONE OF CONCERN CALC GOES ABOVE
    # Size matters. Try to use python primatives.
    # The message gets pickled behind the scenes
    # and appended to the out_queue.
    # This doubling message reached >12k loops in
    # less than 10 min with 60 workers using 15GB ram.

    message.past.append(current)
    return message


def worker(
    in_queue,
    out_queue,
):
    while True:
        task: Message = in_queue.get()
        if task.txt == "END":
            "Ending worker"
            return None
        try:
            out_queue.put(job(task), True, 10)
        except ValueError:
            print("Queue full or message too big!")
        in_queue.task_done()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    workers = os.cpu_count()
    in_queue = multiprocessing.JoinableQueue()
    out_queue = multiprocessing.Queue()
    ### Initial seed. This should be changed.
    # Right now it's a planet.
    in_queue.put(
        Message(
            past=[],
            txt="",
            next=[
                C7(
                    0,
                    125887833.47010976,
                    -23962861.25159607,
                    -5739952.7434539795,
                    7.737457193582947,
                    32.0566489717605,
                    -0.9681529768116752,
                )
            ],
        )
    )
    proclist: list[multiprocessing.Process] = list()
    for i in range(workers):
        proclist.append(
            multiprocessing.Process(target=worker, args=(in_queue, out_queue))
        )
        proclist[-1].start()
    inc = 0
    while True:
        in_queue.join()
        print(f"Loop {inc}")
        inc += 1
        time.sleep(0.01)  # Sleep to accomodate IPC weirdness.
        in_queue.join()
        if out_queue.empty():
            print("Out queue empty. Quiting.")
            break
        while not out_queue.empty():
            ### Message handling should accomodate whatever we're
            ### searching for and the result we want back.
            # In this dummy code we don't record results and just
            # pass the new jobs from the message into the in_queue.
            msg = out_queue.get()
            for i in msg.next:
                in_queue.put(Message(msg.past, "", list((i,))))
            if msg.txt == "END":
                pprint.pp(msg)
                break
    # This politely stops the workers.
    for i in range(workers):
        in_queue.put(Message(past=list(), txt="END", next=list()))
    for i in proclist:
        i.join()
