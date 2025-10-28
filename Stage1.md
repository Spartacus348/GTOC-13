# Definitions:
X, Y, Z, Vx, Vy, Vz, t (Cartesian vector set)

# Pipeline:
1. Generate fine body positions
2. Generate coarse hitboxes
3. Initial Search

# Planet location
## Fine
t, Planet1(X, Y, Z), ... Planet9(X, Y, Z)

t: int (seconds)

dt = 1 day

tMax = 200y = 73050d

## Coarse
t, Planet1((X1Y1Z1), (X2Y2Z2), (X3Y3Z3)), ... Planet9((X1Y1Z1), (X2Y2Z2), (X3Y3Z3))

t: int (seconds)

dt = 1 day

tMax = 200y = 73050d

# Initial Condition Search
## Constants
X = -200

Vy = 0

Vz = 0

## Inputs
Y, Z, Vx, ti

### Precision
100m, 100m, 0.1mm/s

## Outputs
X, Y, Z, Vx, Vy, Vz, ti

## Pre-search reject
- Solar periapsis > gas giant
- Maximal
- Vx < 0
- t > 70yrs (25567.5d)
## Search reject
- no encounter in 70 years
- Calc perihelion, drop at t where coord <= 0.5 au of the sun

## Search Return
- Encounter UNLESS dual encounter w/ yandi+another, then return second

## Search Pseudocode
```
For t from ti to (lesser of 70yrs and perihelion reject) in steps dt:
  find coord at t
  retval = false
  for planetHitbox in database at dt:
    if coord in planetHitbox:
      retval = coord
      if match = yandi, continue
      else break
  return retval
```
