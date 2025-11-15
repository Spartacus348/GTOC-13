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
            past=list(), txt="Job Run", next=runner.empty_nexts(), func=job
        )


# def fake_job(message: Message):
#     current = message.next.pop()[0]
#     r = [current[1], current[2], current[3]] << u.km
#     v = [current[4], current[5], current[6]] << u.km / u.s
#     t = current[0] << u.s

#     ### This will change from Altaira if we're
#     ### in a different sphere of influence.
#     orb = Orbit.from_vectors(Altaira, r, v, epoch=t)

#     #### ACTUAL CONE OF CONCERN CALC GOES BELOW
#     ar, av = orb.propagate(60 * 60 * 24 << u.s).rv()
#     message.next.append(
#         (
#             int(t.value + 60 * 60 * 24),
#             float(ar[0].value),
#             float(ar[1].value),
#             float(ar[2].value),
#             float(av[0].value + 0.1),
#             float(av[1].value),
#             float(av[2].value),
#         )
#     )
#     message.next.append(
#         (
#             int(t.value + 60 * 60 * 24),
#             float(ar[0].value),
#             float(ar[1].value),
#             float(ar[2].value),
#             float(av[0].value),
#             float(av[1].value),
#             float(av[2].value),
#         )
#     )
#     #### ACTUAL CONE OF CONCERN CALC GOES ABOVE
#     # Size matters. Try to use python primatives.
#     # The message gets pickled behind the scenes
#     # and appended to the out_queue.
#     # This doubling message reached >12k loops in
#     # less than 10 min with 60 workers using 15GB ram.

#     message.past.append(current)
#     return message


if __name__ == "__main__":
    multiprocessing.freeze_support()
    in_queue = multiprocessing.JoinableQueue()
    out_queue = multiprocessing.Queue()
    in_queue.put(
        runner.Message(past=list(), txt="Job Run", next=runner.empty_nexts(), func=job)
    )
    runner.runner(in_queue, out_queue)
