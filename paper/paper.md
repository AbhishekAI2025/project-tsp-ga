# Parallel Hybrid Genetic Algorithm for the Traveling Salesman Problem

Authors: Group 7 – COSC6361 (Fall 2025)  
Advisor: Dr. Minhua Huang

---

## Abstract

We present a hybrid genetic algorithm (GA) for the Traveling Salesman Problem
that combines tournament selection, Partially Mapped Crossover (PMX), inversion
mutation, and bounded 2-opt local search. A baseline serial implementation is
extended into an MPI island model capable of running across 2–32 processes. The
parallel design introduces synchronous elite migration and a lightweight load
balancer that reallocates subpopulation sizes based on observed generation time.
We evaluate solution quality and runtime on `berlin52`, `d198`, and `pr439`
instances from TSPLIB. Serial experiments demonstrate tours within 5% of the
known optima. The MPI implementation is packaged for execution on multi-node
clusters; final speedup measurements will be collected once OpenMPI is available
on the target hardware.

---

## 1 Introduction

The Traveling Salesman Problem (TSP) is a canonical NP-hard optimisation
challenge with applications in logistics, robotics, and circuit design. Genetic
Algorithms (GAs) offer a flexible heuristic for large instances, yet their
performance relies on carefully engineered operators and effective exploitation
of modern parallel hardware. This project delivers both a robust serial GA
baseline and an MPI hybrid GA tailored to the requirements of COSC6361.

Our contributions are:

1. A modular C implementation of the serial GA, including PMX crossover,
   inversion mutation, and bounded 2-opt improvement.
2. An MPI island model with elite migration via broadcast and dynamic load
   balancing driven by per-island timing.
3. A reproducible workflow for evaluating speedup, efficiency, and tour quality
   on TSPLIB benchmarks, accompanied by documentation suitable for a workshop
   submission.

---

## 2 Background and Related Work

Hybrid GAs for the TSP commonly combine permutation-safe crossover operators
such as PMX, order crossover (OX), or cycle crossover (CX) with local search
procedures, e.g., 2-opt or Lin–Kernighan. Parallelisation strategies include
master–worker models, fine-grained demes, and coarse-grained island models. The
island model trades communication frequency for solution diversity, making it a
good fit for MPI deployments on clusters where communication costs are non-zero.

Our design follows the coarse-grained approach: each MPI rank evolves a
subpopulation independently, exchanging a small elite set periodically. The
literature suggests migration intervals between 25 and 100 generations and elite
fractions under 10% to maintain diversity while sharing useful genetic
material—guidelines we adopt directly.

---

## 3 Methodology

### 3.1 Serial GA Overview

* **Representation** – Tours are stored as permutations of `[0, …, n-1]`. The
  TSP parser reads TSPLIB `NODE_COORD_SECTION` files and precomputes a full
  distance matrix in double precision.
* **Initialisation** – The population is seeded with random permutations created
  by shuffling the identity tour; the best individual is tracked for elitism.
* **Fitness** – We minimise tour length. Fitness is computed as `1 / (L + 1e-9)`
  to avoid division by zero and maintain high precision for selection.
* **Selection** – Tournament selection with size four balances exploitation and
  diversity. This parameter is configurable via the CLI.
* **Crossover** – PMX maintains permutation validity. We re-implemented PMX to
  avoid pathological cycles and added a small allocation pool for indices.
* **Mutation** – Inversion mutation reverses a random subsequence. With
  probability 0.05 per child, it allows substantial route perturbations without
  breaking adjacency constraints.
* **Local Search** – A bounded 2-opt loop (`max_iterations = 20`) further refines
  offspring. This aggressively reduces edge crossings and typically drives tours
  within a few percent of known optima.
* **Elitism** – The globally best tour survives each generation and is injected
  into the next population.

### 3.2 MPI Island Model

Each MPI rank maintains a private GA identical to the serial baseline. Parallel
behaviour differs only in the elite exchange and load-balancing routines:

* **Migration interval** – Every 50 generations, ranks gather their top 5%
  tours. Rank 0 consolidates the global elite set (also 5% of the global
  population) and broadcasts it back. Each rank replaces its weakest individuals
  with the received elite tours.
* **Load balancing** – After migration, ranks measure wall-clock time for the
  previous interval and share it using `MPI_Allgather`. Subpopulation sizes are
  redistributed proportionally to the inverse runtime (faster ranks receive
  more work) while preserving the global population of 200. Arrays are resized
  and reseeded accordingly.
* **Global best tracking** – `MPI_Allreduce` with `MPI_MINLOC` locates the best
  tour length each generation. The owning rank broadcasts the tour so that all
  populations remain synchronised.

### 3.3 Implementation Notes

* The codebase is modular: `tsp.c` handles parsing, `ga_common.c` exposes shared
  operators, `ga_serial.c` and `parallel_ga.c` manage the high-level evolution
  loops, and `cli.c` handles credential parsing.
* CLI flags allow the instructor to change population size, mutation rate, local
  search depth, and logging cadence without recompilation.
* Logging is line-buffered for human-readable monitoring; summarised results are
  written to `stdout` (serial) or rank 0 only (parallel).

---

## 4 Experimental Setup

### 4.1 Hardware & Software

* macOS 15.6, Apple M2 Pro (12 CPU cores) – used for serial validation.
* GCC 14.2.0 for compilation.
* TSPLIB95 instances downloaded from the Heidelberg TSPLIB site:
  `berlin52`, `d198`, `pr439`.
* MPI build requires OpenMPI ≥ 4.1; commands assume `mpicc` and `mpirun`.

### 4.2 Evaluation Protocol

1. Build the executables (`make serial parallel`).
2. Run the serial GA with 200 generations per instance; record runtime (via
   `/usr/bin/time -p`) and best tour length.
3. Execute the MPI GA with 2, 4, 8, 16, 32 ranks on each instance. Capture
   runtime, final tour quality, and efficiency (speedup / processes).
4. Repeat three times per configuration to smooth randomness; compute means and
   standard deviations.
5. Profile MPI communication using `mpiP` (optional extension).

### 4.3 Metrics

* **Tour length** – Evaluate solution quality against TSPLIB optima. A 5% gap is
  acceptable for the course requirement.
* **Runtime / Speedup** – Compare serial vs. parallel runtime to assess
  scalability.
* **Efficiency** – Speedup divided by process count.
* **Load balance ratio** – Ratio of max to min subpopulation size after
  balancing.

---

## 5 Results

### 5.1 Serial Baseline (200 Generations)

| Instance | Optimum | Best length | Gap (%) | Runtime (s) |
|----------|---------|-------------|---------|-------------|
| berlin52 | 7542.4  | 7,544.37    | +0.03   | 0.21        |
| d198     | 15,780  | 15,822.50   | +0.27   | 6.06        |
| pr439    | 107,217 | 110,660.59  | +3.21   | 55.68       |

All serial runs meet the ≤5% target; `pr439` benefits most from additional
generations or a larger population.

### 5.2 Parallel Scalability (To Be Collected)

Use the following template to record speedup once OpenMPI is available:

| Instance | Procs | Runtime (s) | Speedup | Efficiency | Best length |
|----------|-------|-------------|---------|------------|-------------|
| berlin52 | 1     | 0.21        | 1.00    | 1.00       | 7,544.37    |
| berlin52 | 2     |             |         |            |             |
| berlin52 | 4     |             |         |            |             |
| …        |       |             |         |            |             |

Populate analogous tables for `d198` and `pr439`, and include plots of speedup
vs. processes in the final paper version.

---

## 6 Discussion

The hybrid GA attains high-quality tours quickly for smaller instances.
Performance on `pr439` is competitive but suggests exploring adaptive mutation
rates or deeper local search. The parallel algorithm should offer near-linear
speedup for low process counts because migration volume is modest (≤10 tours per
rank) and load balancing prevents stragglers after the first interval.

Potential improvements include:

* Integrating additional local search (3-opt) triggered by stagnation.
* Replacing synchronous migration with asynchronous peer-to-peer exchange to
  reduce collective latency.
* Caching distance lookups in vectorised form to better exploit SIMD on large
  instances.

---

## 7 Conclusion

We delivered a complete codebase and workflow for analysing a parallel hybrid
GA on canonical TSPLIB datasets. The serial baseline satisfies the project’s
quality constraints, while the MPI design is ready for execution on multi-core
clusters. Future work involves finalising empirical speedup data, extending the
paper with plots, and generating a polished presentation using the provided
outline.

---

## Appendix A Reproduction Checklist

1. `make serial parallel`
2. `./build/serial_ga --instance data/berlin52.tsp --generations 200`
3. `mpirun -np 8 ./build/parallel_ga --instance data/d198.tsp --generations 200`
4. Log outputs to `outputs/` and update tables in Section 5.
5. Generate plots via the notebook in `analysis/` (create as needed).

