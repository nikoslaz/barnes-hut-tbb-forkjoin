# Barnes-Hut N-body Simulation (C++ TBB / Java Fork-Join)

Parallel 2D [Barnes-Hut](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation) gravity simulation,
implemented twice for CS342 (Parallel Programming) at the University of Crete:

- `main.cpp` – C++ with Intel oneTBB (`parallel_for`, `parallel_invoke`)
- `BarnesHutSimulation.java` – Java with the Fork/Join framework and parallel streams

Each step builds a quadtree over `[-R, R] x [-R, R]`, computes the center of mass of every node
in parallel, computes the net force on every body in parallel (opening angle θ = 0.5), and then
updates velocities and positions.

## Build & run

```bash
# C++ (requires oneTBB, e.g. `brew install tbb` or `apt install libtbb-dev`)
g++ -std=c++17 -O2 main.cpp -o barnes_hut -ltbb
./barnes_hut inputs/input1.txt 1000 4 25000

# Java
javac BarnesHutSimulation.java
java BarnesHutSimulation inputs/input1.txt 1000 4 25000
```

Arguments: `input_file iterations [threads] [dt]`

- `threads` – worker threads (default: number of CPU cores)
- `dt` – time step in seconds (default: 1)

The result is written to `output.txt`.

## Input format

```
N            # number of bodies
R            # radius of the universe
x y vx vy mass name   # one line per body (SI units)
```

Sample inputs are in [`inputs/`](inputs/). The report with timings is `Report_4-5.pdf` (in Greek).
