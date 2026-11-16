# Barnes-Hut N-body Simulation (C++ TBB / Java Fork-Join)

2D [Barnes-Hut](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation) gravity simulation,
implemented twice. Assignment 3 for HY-342 (Parallel Programming), University of Crete.

- `main.cpp`: C++ with Intel oneTBB (`parallel_for`, `parallel_invoke`)
- `BarnesHutSimulation.java`: Java with the Fork/Join framework and parallel streams

Each time step:

1. Build a quadtree over the universe `[-R, R] x [-R, R]`.
2. Compute the total mass and center of mass of every node, with subtrees processed in parallel.
3. Compute the net force on every body in parallel. A node far enough away is treated as a
   single mass, using opening angle θ = 0.5.
4. Update velocities and positions in parallel.

## Files

| File | Description |
|---|---|
| `main.cpp` | C++ / TBB implementation |
| `BarnesHutSimulation.java` | Java / Fork-Join implementation |
| `inputs/` | Sample inputs (5 to 1002 bodies) |
| `Report_4-5.pdf` / `Report_4-5.odt` | Report with timings (in Greek) |

## Build & run

```bash
# C++ (requires oneTBB: `brew install tbb` or `apt install libtbb-dev`)
g++ -std=c++17 -O2 main.cpp -o barnes_hut -ltbb
./barnes_hut inputs/input1.txt 1262 4 25000

# Java
javac BarnesHutSimulation.java
java BarnesHutSimulation inputs/input1.txt 1262 4 25000
```

Arguments: `input_file iterations [threads] [dt]`

- `threads`: number of worker threads (default: number of CPU cores)
- `dt`: time step in seconds (default: 1)

The example above simulates one year of the inner solar system (1262 × 25000 s). Afterwards,
Earth is back near its starting position. The result is written to `output.txt`.

## Input / output format

```
N                       # number of bodies
R                       # radius of the universe (m)
x y vx vy mass name     # one line per body, SI units
```

`output.txt` has the same format, with the final positions and velocities.

## Note

The timings in the report were measured with an earlier version of the code, before the force
calculation was fixed.
