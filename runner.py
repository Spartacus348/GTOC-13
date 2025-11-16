#!/usr/bin/env python

import multiprocessing
import os
from collections.abc import Callable
from multiprocessing.queues import Queue, SimpleQueue
from sys import stderr
from typing import Any, NamedTuple, Self

from constants import C7


class Message(NamedTuple):
    past: list[tuple[C7, Any]]
    txt: str
    next: list[tuple[C7, Any]]
    func: Callable[[Self], Self]


def noop(Message) -> Message:
    """
    NO-OP function to satisfy Message definition when Message.Next is empty.
    """
    raise Exception("THIS FUNCTION SHOULD NOT BE CALLED.")
    return Message(past=list(), txt="END", next=list(), func=noop)


def empty_message() -> Message:
    """
    Return a message with empty past and next, no txt, and no-op func.
    """
    return Message(past=list(), txt="", next=list(), func=noop)


def empty_nexts(num: int = 1) -> list[tuple[C7, Any]]:
    return [(C7(0, 0, 0, 0, 0, 0, 0), None)] * num


def end_runner() -> Message:
    return Message(past=list(), txt="END", next=list(), func=noop)


def worker(
    in_queue: SimpleQueue[Message],
    out_queue: Queue[Message],
) -> None:
    while True:
        task: Message = in_queue.get()
        if task.txt == "END":
            break
        try:
            out_queue.put(task.func(task), True, 10)
        except ValueError:
            print("Queue full or message too big!", file=stderr)
    return None


# TODO: Add iterator func.
def runner(in_queue: SimpleQueue[Message], out_queue: Queue[Message]) -> None:
    workers = os.cpu_count() - 1 if os.cpu_count() is not None else 4
    proclist: list[multiprocessing.Process] = list()
    for i in range(workers):
        proclist.append(
            multiprocessing.Process(target=worker, args=(in_queue, out_queue))
        )
        proclist[-1].start()
    jobs = 0
    while True:
        print(f"Completed jobs {jobs}", file=stderr)
        msg = out_queue.get()
        for i in msg.next:
            in_queue.put(Message(past=msg.past, txt=msg.txt, next=[i], func=msg.func))
        if msg.txt == "END":
            break
        jobs += 1
    # This politely stops the workers.
    for i in range(workers):
        in_queue.put(end_runner())
    for i in proclist:
        i.join()
