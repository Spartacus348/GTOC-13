# Overall Approaches

| Approach | Explanation | Pros | Cons |
| ---- | ----- | ----- | ----- |
| Initial parameter search | Randomly iterate initial parameter space, iterating down on good solutions | Guaranteed solution; can reduce search space to initial conditions that interact with planets; comparatively easy to calculate entire trajectory from initial conditions | hard to use solar sail; solution space chaotic at the scale of 200 years (small changes in start not guaranteed to result in small changes in results) | 

# Free Parameters
- y_0, z_0, vx_0, t_0
- solar sail (5% local solar gravity * cos(angle_off_up))

# Toolbox

## Software Libraries
- Boinor: poliastro replacement, orbit objects and plotting: https://github.com/boinor/boinor
- numpy: vectorized functions

## Algorithms
- kepler's equation (orbit -> position as a function of time) 
- Lambert solver (two positions in time -> velocity/orbit) 
- max_displacement_angle = f(approach_velocity, body(gm, radius))
- radius_of_influence(body(mass), distance_from_sun)