"""
Holds the job types
Function signatures:

f(task: Message) -> Message(past: history+self, next: all future jobs)

Current types:
    Search_and_collide:
        iterates an orbit until it runs within a certain radius of any planet
    Search_and_launch:
        tries to fire orbits at various planets, if the resulting trajectory is outside the cone of possibility it tries
         again. Stops after some constant number of sent trajectories.
"""

from runner import Message

N_SENDS: int = 9 * 10


def search_and_collide(task: Message) -> Message:
    pass


def search_and_send(task: Message) -> Message:
    pass
