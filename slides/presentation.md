# Parallel Hybrid GA for TSPLIB TSP Instances

## Slide 1 – Title
* Project title, team members, course, advisor.

## Slide 2 – Motivation
* Importance of TSP; why hybrid GA; need for scalability.

## Slide 3 – Objectives
* Serial baseline, MPI island model, metrics (speedup, efficiency, quality).

## Slide 4 – Dataset Overview
* TSPLIB instances (`berlin52`, `d198`, `pr439`) with size and optimal cost.

## Slide 5 – Genetic Algorithm Recap
* Representation, fitness, selection, crossover, mutation basics.

## Slide 6 – Operator Details
* Tournament selection (k=4), PMX implementation tweaks.

## Slide 7 – Local Search
* Inversion mutation + 2-opt loop, rationale for 20 iterations.

## Slide 8 – Serial Architecture
* Data structures, elitism, logging, configuration via CLI.

## Slide 9 – MPI Island Model
* Process layout, independent evolution, migration cadence.

## Slide 10 – Migration Mechanics
* Elite fraction (5%), broadcast cycle, replacement strategy.

## Slide 11 – Load Balancing
* Timing collection, proportional reallocation, resizing buffers.

## Slide 12 – Implementation Notes
* Code organisation, reusable modules (`tsp`, `ga_common`, `cli`).

## Slide 13 – Experimental Setup
* Hardware, compiler versions, run parameters, measurement tools.

## Slide 14 – Serial Results
* Table of best tour lengths and runtime (use numbers from paper).

## Slide 15 – Parallel Plan
* Expected commands (`mpirun`), ranks tested (2–32), metrics to report.

## Slide 16 – Scalability Expectations
* Discussion of communication volume, expected speedup trend.

## Slide 17 – Challenges
* PMX corner cases, ensuring permutation validity, managing memory.

## Slide 18 – Future Enhancements
* Adaptive operators, asynchronous migration, GPU exploration.

## Slide 19 – Conclusions
* Summary of accomplishments, readiness for experiment completion.

## Slide 20 – Q&A
* Contact info, repository location, acknowledgements.
