# Professor Check-In – Parallel Hybrid GA for TSPLIB TSP

## Slide 1 – Title & Agenda
* *Parallel Hybrid Genetic Algorithm for TSP*, Graduate Algorithms, Advisor: Prof. ___.
* Agenda: recap goals, implementation status, early numbers, risks, next steps.

## Slide 2 – Problem Statement
* Solve TSPLIB instances (`berlin52`, `d198`, `pr439`) with a GA that stays within 5 % of known optima.
* Deliver both a serial baseline and an MPI-based island model that demonstrates speedup and efficiency.
* Produce artifacts for paper + talk: reproducible code, logs, plots, and slides.

## Slide 3 – Success Metrics
* *Solution quality*: average best tour ≤ 1.05 × TSPLIB optimum.
* *Performance*: positive scaling from 2 → 32 ranks; speedup plots + efficiency table.
* *Engineering*: turnkey build, CLI exposure of GA knobs, automated logging for later plotting.

## Slide 4 – Dataset & Evaluation Plan
| Instance | Cities | Notes |
|----------|--------|-------|
| berlin52 | 52     | Standard toy case; sanity-check operators. |
| d198     | 198    | Medium size, stresses 2-opt budget. |
| pr439    | 439    | Large instance to motivate parallelism. |
* Each run: ≥200 generations, seed control for reproducibility, record wall time + best length + migration stats.

## Slide 5 – Serial GA Architecture
* Permutation representation with Euclidean fitness (distance matrix cached).
* Operators: tournament selection (k=4), PMX crossover, inversion mutation, bounded 2-opt local search.
* Elitism keeps global best each generation; CLI exposes `--population`, `--generations`, `--two-opt`, etc.

## Slide 6 – MPI Island Model
* Single binary (`tsp_ga`) switches behavior based on `mpirun` size.
* Rank-local populations evolve independently; every 50 generations elites are gathered on rank 0 then broadcast.
* Load balancer reallocates population sizes using per-interval wall-clock times to keep fast ranks busier.

## Slide 7 – Implementation Progress
* ✅ Serial GA fully functional; logs for all three TSPLIB cases stored in `outputs/`.
* ✅ Parallel GA: migration, elite replacement, dynamic resizing, and MPI-safe CLI already implemented.
* ✅ Makefile builds `tsp_ga` via `mpicc`; serial demos run with `mpirun -np 1` or direct execution.
* 🔄 Pending: capture multi-rank performance numbers, generate plots, integrate profiler (mpiP) hooks.

## Slide 8 – Preliminary Serial Results (200 generations, seed 42)
| Instance | Best length | Runtime (s) |
|----------|-------------|-------------|
| berlin52 | 7 544.37    | 0.21        |
| d198     | 15 822.50   | 6.06        |
| pr439    | 110 660.59  | 55.68       |
* Each best tour meets the ≤5 % quality target; detailed logs show convergence plateaus for tuning discussion.

## Slide 9 – Single-Node Scaling (M3 Air, 8 cores, 200 gens, seed 42)
| Instance | np=2 | np=4 | np=8 |
|----------|------|------|------|
| berlin52 | 0.25 | 0.19 | 0.24 |
| d198     | 3.08 | 1.95 | 2.00 |
| pr439    | 30.70| 18.35| 15.92|
* Modest speedup only on pr439; berlin52/d198 are too short for MPI overhead to amortise.
* No runs above 8 ranks—M3 Air has 4P+4E cores; oversubscribe hurts performance.
* To show clearer scaling, plan to increase generations/population and reduce sync frequency.

## Slide 10 – Upcoming Experiments
* Sweep ranks {2,4,8,16,32} on the lab cluster; capture runtime, best tour, migration count per block.
* Automate plotting (Python script) for speedup/efficiency + convergence.
* Stress-test load balancer by injecting heterogeneous sleep delays to confirm redistribution logic.

## Slide 11 – Risks / Support Needed
* Need confirmed access to OpenMPI nodes with ≥32 ranks; cluster queue limits may bottleneck.
* mpiP instrumentation adds overhead—plan to capture both instrumented and clean runs.
* Large TSPLIB files increase 2-opt cost; may need guidance on acceptable generation count vs. runtime.

## Slide 12 – Demo & Runbook
```
make clean && make CC=mpicc                # produces ./tsp_ga
./tsp_ga data/berlin52.tsp --generations 200                # serial (single core via MPI rank 0)
mpirun -np 8 ./tsp_ga data/d198.tsp --generations 200 --seed 2024     # parallel island model
```
* Logs land on stdout; redirect to `outputs/` for archival.
* Use `--population`, `--mutation`, or `--two-opt` to demonstrate parameter sensitivity live.

## Slide 13 – Closing & Next Steps
* Serial baseline verified, MPI logic implemented—ready to capture scaling curves this week.
* Deliverables in progress: workshop paper draft (`paper/`), refreshed README, and demo-ready code.
* Request feedback on experiment plan + any metrics the professor wants highlighted in final slides.
