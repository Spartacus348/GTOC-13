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

import multiprocessing
import os
import random
from copy import copy

import boinor.iod.izzo
import boinor.twobody.states
import numpy as np
from astropy import units as u

import constants
import runner
from bodyDB import BODY_DB
from constants import C7, Altaira, Sun
from runner import Message
from Stage1 import get_job

PROP_MAX_YR = 70 #years
MIN_SUN_SKIM_AU = 0.05 #au
CART_RANGE = 100 #km
N_SENDS: int = 9 * 10 # sent orbits
FLYBY_TOL = 1e-7 # 0.1 mm/s in km/s

rand = np.random.default_rng()

def planet_reachable(orbit, planet) -> bool:
    sat_r = np.array(orbit.r.to_value(u.km))
    planet_r = np.array(planet[0:2])

    return np.linalg.norm(sat_r - planet_r) > CART_RANGE

def state_to_c7(time, state) -> C7:
    return C7(
        time, state[0], state[1], state[2], state[3], state[4], state[5]
    )

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
    state = task.next.pop()[0]
    orbit = c7_to_orbit(state)
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
                    past = task.past + [(state_to_c7(time, state), p_id)],
                    txt  = str(p_id),
                    next = [(state_to_c7(time, planet), p_id)],
                    func = search_and_launch
                )

    # end of sim
    time = Sun.end_sim_time
    x, y, z = this_orbit.r
    u, v, w = this_orbit.v
    return Message(
        past = task.past + [(state_to_c7(Sun.end_sim_time, state),None)],
        txt = str(time),
        next = [(state_to_c7(time, (x, y, z, u, v, w)),None)],
        func = end_of_run
    )

def search_and_launch(task: Message) -> Message:
    """
    Produces N_SENDS tasks targeting other planets.
    Targets planets and timesteps, finds the trajectory to meet them, then tests if the trajectory is reachable
        through gravity assist. If it is, it adds the trajectory. If not, it drops it.
    :param task: Message containing the approach vector and affecting planet
    :return: Message containing the exit vector and source planet
    """
    #
    state, home_id = task.next.pop()
    orbit = c7_to_orbit(state)

    home_planet_state = state_to_c7(orbit.epoch, BODY_DB[orbit.epoch][home_id])
    home_planet_v = np.array([home_planet_state.u, home_planet_state.v, home_planet_state.w])
    approach_rel_vel = np.array(orbit.v.to_value(u.km/u.s)) - home_planet_v
    approach_speed = np.linalg.norm(approach_rel_vel)

    planet_options = list(range(len(BODY_DB[orbit.epoch])))
    planet_options.remove(home_id)
    time_options = [key for key in BODY_DB.keys() if orbit.epoch < key < orbit.epoch + 70 * Sun.day * Sun.year]
    future_tasks = []
    while len(future_tasks) < N_SENDS:
        planet = random.choice(planet_options)
        time = random.choice(time_options)
        v0, _ = boinor.iod.izzo.lambert(
            k= Sun.gravity,
            r0 = orbit.r,
            r = planet[0:3] << u.km,
            tof= (time-orbit.epoch) << u.s,
        )
        exit_rel_vel = np.array(v0.to_value(u.km/u.s)) - home_planet_v
        # determine if this exit rel vel is achievable given the approach vel and radius range
        # first filter: determine if the two velocity vectors are the same magnitude
        if np.abs(np.linalg.norm(exit_rel_vel) - approach_speed) > FLYBY_TOL:
            continue

    pass


def end_of_run(task: Message) -> Message:
    print(json.dumps(task.past))
    return empty_message()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    in_queue = multiprocessing.JoinableQueue()
    out_queue = multiprocessing.Queue()
    workers = os.cpu_count() - 1 if os.cpu_count() is not None else 4
    for i in range(workers):
        in_queue.put(get_job(search_and_collide))
    print("t,x,y,z,u,v,w,planetID")
    runner.runner(in_queue, out_queue)
