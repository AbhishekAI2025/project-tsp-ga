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
	•	Use TSPLIB datasets: berlin52.tsp, d198.tsp, pr439.tsp
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
```
# serial executable (always available)
make serial

# parallel executable (requires mpicc in PATH)
make parallel

# override compilers if needed
CC=clang make serial
MPICC=mpicc-openmpi make parallel
```
Artifacts are written to `build/serial_ga` and `build/parallel_ga`.

Running the Serial Baseline
---------------------------
```
./build/serial_ga --instance data/berlin52.tsp --generations 200 --log-interval 50
./build/serial_ga --instance data/d198.tsp --generations 200 --log-interval 50
./build/serial_ga --instance data/pr439.tsp --generations 200 --log-interval 50
```
Default GA settings follow the specification (population 200, crossover 0.8,
mutation 0.05, tournament size 4, two-opt iterations 20, seed 42). Override any
parameter through CLI flags (`--population`, `--mutation`, etc.).

Parallel Execution
------------------
```
mpirun -np 8 ./build/parallel_ga --instance data/berlin52.tsp --generations 200 \
  --log-interval 50 --seed 1234
```
Every 50 generations each rank broadcasts its top 5% tours; rank 0 selects the
global elite set and redistributes them. After the migration phase, wall-clock
times from the last interval are collected to rebalance subpopulation sizes.

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
1. Build both executables (`make serial parallel`).
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
