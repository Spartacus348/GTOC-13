#!/usr/bin/env python

import multiprocessing
import os
import pprint
import time
from collections.abc import Callable
from multiprocessing.queues import JoinableQueue, Queue
from sys import stderr
from typing import NamedTuple, Self

from astropy import units as u
from boinor.twobody import Orbit

from constants import C7, Altaira, UnnamedTuple


class Message(NamedTuple):
    past: list[UnnamedTuple]
    txt: str
    next: list[UnnamedTuple]
    func: Callable[[Self], Self]


def __noop(Message) -> Message:
    raise Exception("THIS FUNCTION SHOULD NOT BE CALLED.")
    return Message(past=list(), txt="END", next=list(), func=__noop)


def end_runner() -> Message:
    return Message(past=list(), txt="END", next=list(), func=__noop)


def fake_job(message: Message):
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
    in_queue: JoinableQueue[Message],
    out_queue: Queue[Message],
):
    while True:
        task: Message = in_queue.get()
        if task.txt == "END":
            "Ending worker"
            return None
        try:
            out_queue.put(task.func(task), True, 10)
        except ValueError:
            print("Queue full or message too big!", file=stderr)
        in_queue.task_done()


def runner(in_queue: JoinableQueue[Message], out_queue: Queue[Message]) -> None:
    workers = os.cpu_count() - 1 if os.cpu_count() is not None else 4
    proclist: list[multiprocessing.Process] = list()
    for i in range(workers):
        proclist.append(
            multiprocessing.Process(target=worker, args=(in_queue, out_queue))
        )
        proclist[-1].start()
    inc = 0
    while True:
        in_queue.join()
        print(f"Loop {inc}", file=stderr)
        inc += 1
        time.sleep(0.01)  # Sleep to accomodate IPC weirdness.
        in_queue.join()
        if out_queue.empty():
            print("Out queue empty. Quiting.", file=stderr)
            break
        while not out_queue.empty():
            ### Message handling should accomodate whatever we're
            ### searching for and the result we want back.
            # In this dummy code we don't record results and just
            # pass the new jobs from the message into the in_queue.
            msg = out_queue.get()
            for i in msg.next:
                in_queue.put(
                    Message(past=msg.past, txt=msg.txt, next=list((i,)), func=msg.func)
                )
            if msg.txt == "END":
                pprint.pp(msg)
                break
    # This politely stops the workers.
    for i in range(workers):
        in_queue.put(end_runner())
    for i in proclist:
        i.join()


if __name__ == "__main__":
    multiprocessing.freeze_support()
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
            func=fake_job,
        )
    )
    runner(in_queue=in_queue, out_queue=out_queue)
