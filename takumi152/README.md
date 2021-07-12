# ICFPC 2021 solver and tools (takumi152)

## Solvers

Two types of solver was written by me during the contest.

### Simulated annealing solver (`icfpc2021.cpp`)

This solver uses simulated annealing to optimize vertex placement.  
The solver consists of two steps:

1. Solve the problem in real number domain (i.e. ignores a constraint that vertex coordinates must be integer)
2. Round coordinates of all vertex to integer, then solve the problem in integer domain.

After certain time have passed (which is defined as `time_limit` in source code), the solver terminates and outputs the best solution.

### Brute force solver (`icfpc2021_brute_force.cpp`)

This solver will search all possible vertex placements (using backtrack search), and returns the best solution that was found.

## Jupyter Notebook tools

For convenience, I have created some tools as Jupyter Notebooks.

### `get_problem.ipynb`

Get problems from server using HTTP API.

### `submit-direct.ipynb`

Submit solutions to server using HTTP API.

### `try_all.ipynb`

A simple code for solving all problems using the simulated annealing solver.

### `try_and_internal_submit.ipynb`

Repeatedly solve each problems using the simulated annealing solver, keep the best results and send the result to our team's internal leaderboard server.

#### Side notes

* The input should be provided to each solver through standard input.
  * Both solver uses text format converted from problem json files. You can convert problem json files to text using `../server/problem2txt.py`.
