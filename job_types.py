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
from copy import copy

import numpy as np

import runner
from constants import Altaira, Sun, UnnamedTuple
from runner import Message
import boinor.twobody.states
from astropy import units as u

import constants
from bodyDB import BODY_DB


PROP_MAX_YR = 70 #years
MIN_SUN_SKIM_AU = 0.05 #au
CART_RANGE = 100 #km
N_SENDS: int = 9 * 10 # sent orbits

def planet_reachable(orbit, planet) -> bool:
    sat_r = np.array(orbit.r.to_value(u.km))
    planet_r = np.array(planet[0:2])

    return np.linalg.norm(sat_r - planet_r) > CART_RANGE


def state_to_tuple(time: int, state) -> constants.UnnamedTuple:
    return time, state[0], state[1], state[2], state[3], state[4], state[5]


def c7_to_orbit(state: constants.C7) -> boinor.twobody.Orbit:
    return boinor.twobody.Orbit.from_vectors(
        attractor=Altaira,
        r= [state.x, state.y, state.z] << u.km,
        v= [state.u, state.v, state.w] << u.km/u.s,
        epoch = state.t,
    )


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
    this_orbit = copy(orbit)
    for time, planets in BODY_DB.items():
        # skip timesteps in the past
        if time < orbit.epoch:
            continue
        # use this timestep to find where we're starting from
        if time == orbit.epoch:
            for p_id, planet in enumerate(planets):
                planet_reachable(orbit, planet)
                if planet_reachable(orbit, planet):
                    home_planet_id = p_id
                    continue
        # end if we run out of time
        if (time-orbit.epoch)/(Sun.day * Sun.year) > PROP_MAX_YR:
            return runner.empty_message()

        this_orbit = orbit.propagate(time - this_orbit.epoch)
        # end if we fly within 0.05 AU
        if np.linalg.norm(np.array(this_orbit.r))/Sun.au < MIN_SUN_SKIM_AU:
            return runner.empty_message()

        for p_id, planet in enumerate(planets):
            if p_id == home_planet_id:
                continue
            if planet_reachable(this_orbit, planet):
                return Message(
                    past = task.past + [state_to_tuple(time, state)],
                    txt  = str(time),
                    next = [state_to_tuple(time, planet)],
                    func = search_and_launch
                )

    # end of sim
    time = Sun.end_sim_time
    x, y, z = this_orbit.r
    u, v, w = this_orbit.v
    return Message(
        past = task.past + [state_to_tuple(-1, state)],
        txt = str(time),
        next = [(time, x, y, z, u, v, w)],
        func = end_of_run
    )


def search_and_launch(task: Message) -> Message:
    pass

def end_of_run(task: Message) -> Message:
    pass