# ICFP Programming Contest 2021

## foota's solver

## Build

```
$ make
```

## Usage

```
burain_wall_solver [time=300.0] < input_file
```

Output is vertices of the `pose`.

* example

```
$ ./burain_wall_solver 180 < problem_001.in
```

`input_file` is generated from a problem json file by `problems/convet_bonus_format.py`.

## Solutions

* A figure is translated to the center of the hole.
* Vertices of the figure outside the hole move to nearest edges of the hole.
* If the vertex of the figure is in the hole, the edge is in the hole, or the distance of the edge is within the limit, update as follows:
  * Stretching or compressing the distance between the sides to be within the limit.
  * Searching the best point by the evalulation function (energy).
  * With some probability, other than the best point will be selected.
  * With some probability, the vertex of figure will move to the vertex of hole.
