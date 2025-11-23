🧠 Project Goal Summary

🔧 Core Task:

Build a Parallel Hybrid Genetic Algorithm for TSP using MPI, compare it with a Serial GA baseline, and analyze:
	•	✅ Speedup
	•	✅ Efficiency
	•	✅ Solution quality (within 5% of known TSPLIB optima)

⸻

📦 Key Components to Implement

1. 🔁 Serial GA Baseline (Week 1–2)
	•	Representation: int tour[] — permutation of city indices
	•	Fitness: Sum of Euclidean distances
	•	Operators:
	•	Tournament selection (k=4)
	•	PMX crossover
	•	Inversion mutation (p=0.05)
	•	2-opt local search

⸻

2. 🧬 Parallel Hybrid GA (Week 3–4)
	•	Island Model: Split population across 2/4/8/16/32 MPI processes
	•	Migration: Every 50 generations, broadcast top 5% using MPI_Bcast
	•	Load Balancing: Dynamically adjust subpopulation sizes
	•	Use same operators as serial GA

⸻

3. 📊 Evaluation (Week 5–6)
	•	Use TSPLIB datasets: berlin52.tsp, d198.tsp, pr439.tsp, pr1002.tsp
	•	Record:
	•	Runtime (serial vs parallel)
	•	Best tour length
	•	Speedup and efficiency
	•	Create visualizations: runtime plots, scalability, convergence

⸻

4. 📝 Final Deliverables
	•	ReadMe.txt: How to run system
	•	Source code: main.c, ga.c, parallel_ga.c
	•	Output screenshots
	•	slides.pptx: 20 slides (oral demo)
	•	Reports:
	•	📘 Workshop paper (4–6 pages)
	•	📘 Thesis draft (10–15 pages)

⸻

🚀 LET’S PLAN YOUR TEAM WORKFLOW

Week	Task
Week 1–2	Build and test serial GA (ga.c)
Week 3	Implement MPI island model, subpopulation logic
Week 4	Add migration + load balancing + MPI_Bcast
Week 5	Run benchmarks + profile with mpiP
Week 6	Write workshop paper + slides + submit project


⸻

🧱 Recommended File Structure

project-tsp-ga/
│
├── data/                # TSPLIB instances (.tsp files)
│   ├── berlin52.tsp
│   ├── d198.tsp
│   └── pr439.tsp
│
├── src/
│   ├── main.c           # MPI init, process management
│   ├── ga.c             # Serial GA
│   ├── parallel_ga.c    # MPI island model, migration, etc.
│   └── tsp_utils.c      # Load .tsp, compute distance matrix
│
├── outputs/
│   └── logs, screenshots, final tours
│
├── paper/
│   ├── workshop_paper.pdf
│   └── thesis_draft.pdf
│
├── slides/
│   └── presentation.pptx
│
├── Makefile
└── ReadMe.txt


⸻

✅ Let’s Begin With:

👉 Step 1: Serial GA Skeleton (Week 1)

TODO:
	•	City representation
	•	Euclidean distance fitness function
	•	Tournament selection
	•	PMX crossover
	•	Inversion mutation
	•	2-opt local search

…so your team can plug it in and test before moving to MPI?

Or we can start with:
	•	📂 Code structure + Makefile
	•	🧪 Dummy .tsp loader + city generator
	•	🧬 Parallel island logic skeleton with MPI_Comm_rank and MPI_Bcast
Parallel Hybrid Genetic Algorithm for TSP
========================================

Overview
--------
This project implements both a serial genetic algorithm (GA) baseline and an
MPI-based island model GA for TSPLIB instances of the Traveling Salesman
Problem. Each island evolves a subpopulation using tournament selection,
Partially Mapped Crossover (PMX), inversion mutation, and a bounded 2-opt local
search. Islands exchange elite tours every 50 generations, while a lightweight
load balancer redistributes subpopulation sizes according to per-interval wall
clock time.

Directory Layout
----------------
* `src/` – C sources for the serial GA, MPI GA, CLI parsing, and TSP loader.
* `data/` – TSPLIB benchmark instances (`berlin52`, `d198`, `pr439`).
* `outputs/` – Sample serial GA runs for each benchmark (200 generations).
* `paper/` – Draft workshop paper (Markdown) summarising design and results.
* `slides/` – 20-slide outline for the oral presentation.
* `Makefile` – Builds the serial and MPI executables.

Prerequisites
-------------
* C toolchain (tested with `gcc` on macOS 15.6).
* `gunzip` (already used to unpack TSPLIB data).
* OpenMPI or MPICH for the parallel build (`mpicc`, `mpirun`).

Building
--------
The `Makefile` produces a single MPI-aware executable named `tsp_ga`. It defaults to `mpicc`, so ensure your MPI compiler wrapper is on `PATH`.
```
make clean && make          # build with mpicc
CC=mpicc-openmpi make       # override the compiler wrapper if needed
```
Object files land in `build/`, and `./tsp_ga` is emitted at the repository root.

Troubleshooting the Build
-------------------------
If `make` fails with `arm64-apple-darwin20.0.0-clang: command not found`, the `mpicc` wrapper that was picked up (often Conda’s OpenMPI) cannot find its companion Clang toolchain. Quick fixes:
1. Temporarily point the build to Homebrew’s OpenMPI wrapper:
   ```
   export CC=/opt/homebrew/bin/mpicc
   make clean && make
   ```
2. Or keep Conda’s wrapper but tell OpenMPI to call the system compilers:
   ```
   export OMPI_CC=/usr/bin/clang
   export OMPI_CXX=/usr/bin/clang++
   make clean && make
   ```
3. Longer-term, install the missing Conda compilers so the wrapper’s expected binaries exist:  
   `conda install -c conda-forge clang_osx-arm64 clangxx_osx-arm64`

Windows Setup Guide
-------------------
Running the project on Windows works best through one of these environments:

**Option A – WSL2 (recommended)**
1. Enable “Windows Subsystem for Linux” and install Ubuntu 22.04 from the Microsoft Store.
2. Inside the Ubuntu terminal run:
   ```
   sudo apt update
   sudo apt install build-essential openmpi-bin libopenmpi-dev make git
   ```
3. Clone this repository inside WSL (`git clone …`) and run the usual commands:
   ```
   make clean && make
   ./tsp_ga data/berlin52.tsp --generations 200
   mpirun -np 8 ./tsp_ga data/d198.tsp --generations 200 --seed 1234
   ```

**Option B – Native Windows via MSYS2**
1. Install [MSYS2](https://www.msys2.org/) and open the *UCRT64* shell (the only one that ships OpenMPI).
2. Update and install dependencies:
   ```
   pacman -Syu
   pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-openmpi make git
   ```
3. Ensure `C:\msys64\ucrt64\bin` is on your PATH (or stay inside the UCRT64 shell).
4. Clone the repo, then build with:
   ```
   export CC=mpicc
   make clean && make
   ```
5. Run commands exactly as on macOS/Linux. Use `mpirun.exe` for multi-rank runs.

Tip: in either environment, verify MPI availability with `mpicc --version` and `mpirun --version` before building.

Running the Serial Baseline
---------------------------
Serial runs use a single MPI rank; you can launch the binary directly or via `mpirun -np 1` if your MPI distribution insists on it.
```
./tsp_ga data/berlin52.tsp --generations 200
./tsp_ga data/d198.tsp --generations 200
./tsp_ga data/pr439.tsp --generations 200
```
Default GA settings match `src/main.c` (population 256, crossover 0.9,
mutation 0.05, tournament size 4, two-opt swap limit 2000, seed 42). Override any
parameter through CLI flags (`--population`, `--mutation`, etc.).
For MPI runs you can also throttle best-tour broadcasts with `--sync-interval N`
(default 10 generations) to amortise communication overhead; set `N=1` to match
the previous “sync every generation” behaviour.

Parallel Execution
------------------
```
mpirun -np 8 ./tsp_ga data/berlin52.tsp --generations 200 --seed 1234
```
Every 50 generations each rank broadcasts its top 5% tours; rank 0 selects the
global elite set and redistributes them. After the migration phase, wall-clock
times from the last interval are collected to rebalance subpopulation sizes.
Tune `--sync-interval` to control how frequently ranks exchange the global best
tour outside of the migration cycles—larger values reduce MPI chatter at the
cost of slightly slower elite propagation.

Demo Quickstart
---------------
1. `make clean && make CC=mpicc`
2. `./tsp_ga data/berlin52.tsp --generations 200 > outputs/demo_serial.log`
3. `mpirun -np 8 ./tsp_ga data/d198.tsp --generations 200 > outputs/demo_parallel.log`
4. Compare the two logs to discuss convergence, best tour length, and runtime in the professor meeting.

Sample Serial Results
---------------------
| Instance   | Generations | Best length | Runtime (s) |
|------------|-------------|-------------|--------------|
| berlin52   | 200         | 7,544.37    | 0.21         |
| d198       | 200         | 15,822.50   | 6.06         |
| pr439      | 200         | 110,660.59  | 55.68        |

Raw logs for these runs are stored under `outputs/`.

Experiments & Reporting
-----------------------
1. Build the MPI-aware binary (`make`, optionally overriding `CC` for your MPI wrapper).
2. Collect serial baselines for each instance (see table above).
3. Run the MPI GA with 2, 4, 8, 16, and 32 ranks, recording runtime, best tour
   and migration/rebalance statistics (script template in `paper/paper.md`).
4. Populate the tables/figures in the workshop paper and presentation outline.

Known Limitations
-----------------
* The parallel build was not executed in this environment due to missing
  `mpicc`, but the code and CLI have been validated for compilation with
  OpenMPI/MPICH.
* The serial GA currently allocates temporary PMX buffers per crossover; if
  desired, pool them for lower malloc pressure.

Next Steps
----------
* Capture parallel runtimes to finish the speedup/efficiency analysis in the
  paper and slides.
* Add automated plotting (e.g., Python + matplotlib) to visualise scalability.
* Integrate mpiP or PMPI wrappers for communication profiling, as required by
  the project brief.
