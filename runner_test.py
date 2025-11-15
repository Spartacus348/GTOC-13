#!/usr/bin/env python
import multiprocessing
import pprint

import TreeTrack


def job(msg: TreeTrack.Message) -> TreeTrack.Message:
    pprint.pp(msg)
    if msg.txt == "Job Run":
        print("Second round")
        return TreeTrack.end_runner()
    else:
        print("First round")
        return TreeTrack.Message(
            past=list(), txt="Job Run", next=list(((0, 0, 0, 0, 0, 0, 0),)), func=job
        )


if __name__ == "__main__":
    multiprocessing.freeze_support()
    in_queue = multiprocessing.JoinableQueue()
    out_queue = multiprocessing.Queue()
    in_queue.put(
        TreeTrack.Message(
            past=list(), txt="Job Run", next=list(((0, 0, 0, 0, 0, 0, 0),)), func=job
        )
    )
    TreeTrack.runner(in_queue, out_queue)
