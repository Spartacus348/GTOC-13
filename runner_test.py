#!/usr/bin/env python
import multiprocessing
import pprint

import runner


def job(msg: runner.Message) -> runner.Message:
    pprint.pp(msg)
    if msg.txt == "Job Run":
        print("Second round")
        return runner.end_runner()
    else:
        print("First round")
        return runner.Message(
            past=list(), txt="Job Run", next=list(((0, 0, 0, 0, 0, 0, 0),)), func=job
        )


if __name__ == "__main__":
    multiprocessing.freeze_support()
    in_queue = multiprocessing.JoinableQueue()
    out_queue = multiprocessing.Queue()
    in_queue.put(
        runner.Message(
            past=list(), txt="Job Run", next=list(((0, 0, 0, 0, 0, 0, 0),)), func=job
        )
    )
    runner.runner(in_queue, out_queue)
