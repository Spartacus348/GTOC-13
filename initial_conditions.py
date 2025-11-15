##########################################
# code for returning the initial orbits
#
# Gabe
##########################################
import csv
import os

import numpy as np
from astropy import units as u

from constants import Altaira, Sun


def check_constraints(guess, d1, target, x0):
    # checks various bounds and returns a guess at the edge of that boundary
    u0_inner = guess[0]
    d0_inner = guess[1]
    mu = Sun.gravity

    r0 = np.hypot(d0_inner, x0)
    sma = 1 / ((2 / r0) - (u0_inner**2 / mu))

    # first check: eccentricity has to be > 1
    if sma > 0:
        v_parabola = np.sqrt(2 * mu / r0)
        u0_inner = v_parabola * 1.01

        # recalculate value
        sma = 1 / ((2 / r0) - (u0_inner**2 / mu))

    # second check: periapse has to be below the target altitude
    p = (u0_inner * d0_inner) ** 2 / mu
    ecc = np.sqrt(1 - p / sma)
    periapse = -sma * (ecc - 1)
    r1 = np.hypot(d1, target[3])

    return np.array([u0_inner, d0_inner])


def from_target(
    target: np.ndarray,
    initial_t: float = 0,
    initial_y: float = None,
    initial_z: float = None,
    initial_u: float = None,
) -> np.ndarray:
    """
    Returns the initial position and velocity for an orbit that reaches
        the target position at the target time, leaving at initial time
    :param target: position in time to target
    :param initial_t: initial time, optional, defaults to 0
    :param initial_y: initial velocity guess, optional
    :param initial_z: initial velocity guess, optional
    :param initial_u: initial velocity guess, optional
    :return: the initial state vector for the orbit
    """
    # given an input position in time, returns a state vector that will reach the target given the following restrictions:
    # r0 = [_, -200 [AU], _, _, _, 0 [km/s], 0 [km/s]]
    # uses python-constraint to solve two constraining equations
    # first convert the problem into a 2d solution, distance along yz => d
    roll_radians = np.arctan2(target[3], target[2])  # opposite is z, adjacent is y
    d1 = np.hypot(target[3], target[2])
    x0_km = -200 * Sun.au

    def trajectory_constraint(u0_inner, d0_inner, mu, x0, d1_inner, x1) -> float:
        return (
            ((d0_inner * u0_inner) ** 2) / mu
            - np.hypot(x1, d1_inner)
            - (
                (x0 * x1 - d1_inner * d0_inner) / np.hypot(x0, d0_inner)
                + d0_inner * d1_inner * u0_inner * u0_inner / mu
            )
        )

    def time_constraint(u0_inner, d0_inner, mu, x0, r1, side_of_periapse, dt) -> float:
        r0 = np.hypot(d0_inner, x0)
        sma = 1 / ((2 / r0) - (u0_inner**2 / mu))
        semilatus_rectum = ((d0_inner * u0_inner) ** 2) / mu
        ecc = np.sqrt(1 - semilatus_rectum / sma)
        t0 = np.arccos((semilatus_rectum / r0 - 1) / ecc)
        if pos0 := (semilatus_rectum / r1 - 1) / ecc > 1:
            t1 = np.arccos(pos0)
        else:
            t1 = 0

        h0 = -2 * np.arctanh(np.sqrt((ecc - 1) / (ecc + 1)) * np.tan(t0 / 2))
        h1 = 2 * np.arctanh(np.sqrt((ecc - 1) / (ecc + 1)) * np.tan(t1 / 2))

        if side_of_periapse:
            h1 *= -1

        return (
            np.sqrt(np.abs(sma) ** 3 / mu)
            * ((ecc * np.sinh(h1) - h1) - (ecc * np.sinh(h0) - h0))
            - dt
        )

    def traj_derivatives(
        u0_inner, d0_inner, mu, x0, d1_inner, x1
    ) -> tuple[float, float]:
        dtraj_du = 2 * u0_inner * d0_inner * (d0_inner - d1_inner) / mu

        r0 = np.hypot(d0_inner, x0)
        dtraj_dd = (
            (2 * u0_inner * u0_inner * d0_inner / mu)
            + (x0 * x1 * d0_inner / r0**3)
            - (d1_inner * (1 - (d0_inner / r0) ** 2) / r0)
            + u0_inner * u0_inner * d1_inner / mu
        )
        return dtraj_du, dtraj_dd

    def time_derivatives(
        u0_inner, d0_inner, mu, x0, r1, side_of_periapsis, _
    ) -> tuple[float, float]:
        de_du = u0_inner

        r0 = np.hypot(d0_inner, x0)
        de_dd = mu * d0_inner / r0**3

        energy = u0_inner * u0_inner / 2 - mu / r0
        sma = 1 / ((2 / r0) - (u0_inner**2 / mu))
        da_du = -2 * sma * sma * u0_inner / mu
        da_dd = -2 * sma * sma * d0_inner / (r0**3)

        p = ((d0_inner * u0_inner) ** 2) / mu
        dp_du = 2 * p / u0_inner
        dp_dd = 2 * p / d0_inner

        ecc = np.sqrt(1 - p / sma)
        decc_du = 2 * (p * de_du + energy * dp_du) / (mu * mu * ecc)
        decc_dd = 2 * (p * de_dd + energy * dp_dd) / (mu * mu * ecc)

        m_dot = np.sqrt(np.abs(sma) ** 3 / mu)
        t0 = np.arccos((p / r0 - 1) / ecc)
        if pos0 := (p / r1 - 1) / ecc > 1:
            t1 = np.arccos(pos0)
        else:
            t1 = 0

        if np.isnan(t1):
            t1 = 0

        h0 = -2 * np.arctanh(np.sqrt((ecc - 1) / (ecc + 1)) * np.tan(t0 / 2))
        h1 = 2 * np.arctanh(np.sqrt((ecc - 1) / (ecc + 1)) * np.tan(t1 / 2))

        if side_of_periapsis:
            h1 *= -1

        m1 = ecc * np.sinh(h1) - h1
        m0 = ecc * np.sinh(h0) - h0
        tau = m_dot * (m1 - m0)

        dmd_du = 3 * tau * da_du / (2 * sma * (m1 - m0))
        dmd_dd = 3 * tau * da_dd / (2 * sma * (m1 - m0))

        dh0_du = -(da_du * (ecc * np.cosh(h0) - 1) + sma * np.cosh(h0) * decc_du) / (
            sma * ecc * np.sinh(h0)
        )
        dh0_dd = -(da_dd * (ecc * np.cosh(h0) - 1) + sma * np.cosh(h0) * decc_dd) / (
            sma * ecc * np.sinh(h0)
        )

        dm0_du = decc_du * np.sinh(h0) + dh0_du * (ecc * np.cosh(h0) - 1)
        dm0_dd = decc_dd * np.sinh(h0) + dh0_dd * (ecc * np.cosh(h0) - 1)

        if h1 == 0:
            dh1_du = 0
            dh1_dd = 0
        else:
            dh1_du = -(
                da_du * (ecc * np.cosh(h1) - 1) + sma * np.cosh(h1) * decc_du
            ) / (sma * ecc * np.sinh(h1))
            dh1_dd = -(
                da_dd * (ecc * np.cosh(h1) - 1) + sma * np.cosh(h1) * decc_dd
            ) / (sma * ecc * np.sinh(h1))

        dm1_du = decc_du * np.sinh(h1) + dh1_du * (ecc * np.cosh(h1) - 1)
        dm1_dd = decc_dd * np.sinh(h1) + dh1_dd * (ecc * np.cosh(h1) - 1)

        dtau_du = dmd_du * (m1 - m0) + m_dot * (dm1_du - dm0_du)
        dtau_dd = dmd_dd * (m1 - m0) + m_dot * (dm1_dd - dm0_dd)

        return dtau_du, dtau_dd

    # initial guess is close enough in to pass within the target orbit,
    # with enough speed to definitely be hyperbolic
    vel = (
        np.sqrt(40 + np.abs(2 * Sun.gravity / x0_km))
        if initial_u is None
        else initial_u
    )
    dist = (
        0.9 * d1
        if initial_y is None or initial_z is None
        else np.hypot(initial_y, initial_z)
    )
    guess = np.array([vel, dist])
    before_periapse = (
        target[2] < 0
    )  # if y-component is < 0 it's on the rising side of periapsis

    rtol = np.array([1e-10, 1e-10])
    error = np.array(
        [
            trajectory_constraint(
                guess[0], guess[1], Sun.gravity, x0_km, d1, target[1]
            ),
            time_constraint(
                guess[0],
                guess[1],
                Sun.gravity,
                x0_km,
                d1,
                before_periapse,
                target[0] - initial_t,
            ),
        ]
    )

    dx = newton_solve(
        d1, error, guess, target, time_derivatives, traj_derivatives, x0_km, initial_t
    )
    while np.all(np.abs(dx / guess) > rtol):
        guess += dx
        guess = check_constraints(guess, d1, target, x0_km)

        error = np.array(
            [
                trajectory_constraint(
                    guess[0], guess[1], Sun.gravity, x0_km, d1, target[1]
                ),
                time_constraint(
                    guess[0],
                    guess[1],
                    Sun.gravity,
                    x0_km,
                    d1,
                    before_periapse,
                    target[0] - initial_t,
                ),
            ]
        )

        dx = newton_solve(
            d1,
            error,
            guess,
            target,
            time_derivatives,
            traj_derivatives,
            x0_km,
            initial_t,
        )

    u0, d0 = guess[0], guess[1]
    y0, z0 = d0 * np.cos(roll_radians), d0 * np.sin(roll_radians)

    return np.array([0, x0_km, y0, z0, u0, 0, 0])


def newton_solve(
    d1, error, initial_guess, target, time_derivatives, traj_derivatives, x0_km, t0
) -> np.ndarray:
    dtraj_du, dtraj_dd = traj_derivatives(
        initial_guess[0], initial_guess[1], Sun.gravity, x0_km, d1, target[1]
    )

    dtime_du, dtime_dd = time_derivatives(
        initial_guess[0],
        initial_guess[1],
        Sun.gravity,
        x0_km,
        d1,
        target[1],
        target[0] - t0,
    )

    j_matrix = np.array([[dtraj_du, dtraj_dd], [dtime_du, dtime_dd]])
    dx = np.linalg.solve(j_matrix, error)
    return dx


def test_from_target():
    from boinor.twobody import Orbit

    output = {}
    with open(
        os.path.join("ext-resources", "gtoc13_planets.csv"),
        newline="",
        encoding="iso-8859-1",
    ) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            output[row["Name"]] = row.copy()

    row = output["Yavin"]
    target = Orbit.from_classical(
        attractor=Altaira,
        a=float(row["Semi-Major Axis (km)"]) * u.km,
        ecc=float(row["Eccentricity ()"]) * u.one,
        inc=float(row["Inclination (deg)"]) * u.deg,
        raan=float(row["Longitude of the Ascending Node (deg)"]) * u.deg,
        argp=float(row["Argument of Periapsis (deg)"]) * u.deg,
        nu=float(row["Mean Anomaly at t=0 (deg)"]) * u.deg,
    )

    target_time = target.epoch + 30 * Sun.day * Sun.year * u.s

    propagated_orbit = target.propagate(target_time)
    print(propagated_orbit)
    pos = propagated_orbit.r.to_value(u.km)
    print(pos)

    target_input = np.array([20 * Sun.day * Sun.year, pos[0], pos[1], pos[2]])
    output = from_target(target_input)
    print(output)
