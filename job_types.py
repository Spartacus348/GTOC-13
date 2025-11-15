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
import boinor.twobody.states
import numpy as np
from astropy import units as u

import constants
from bodyDB import BODY_DB

N_SENDS: int = 9 * 10

def planet_reachable(orbit, planet):
    pass


def search_and_collide(task: Message) -> Message:
    """
    Propagates an orbit forward until it runs within a certain radius of any planet. Then it generates a new message
        calling search_and_send with itself as the most recent historical message
    :param task: Message containing its history and the orbit
    :return: Message containing the intersection state at the next planet
    """
    # iterate through the states at the t_step
    state = task.next.pop()
    orbit = c7_to_orbit(constants.name_tuple(state))
    home_planet_id = None
    for time, planets in BODY_DB.items():
        if time < orbit.epoch:
            continue
        if time == orbit.epoch:
            for id, planet in enumerate(planets):
                planet_reachable(orbit, planet)
                if planet_reachable(orbit, planet):
                    home_planet_id = id
                    continue
        this_orbit = orbit.propagate(time)
        for id, planet in enumerate(planets):
            if id == home_planet_id:
                continue
            if planet_reachable(orbit, planet):
                return Message(
                    past = task.past + [state_to_tuple(time, state)],
                    txt  = str(time),
                    next = [state_to_tuple(time, planet)],
                    func = search_and_send
                )


def search_and_send(task: Message) -> Message:
    pass
