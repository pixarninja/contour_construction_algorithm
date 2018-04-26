# Contour Construction Algorithm

## Usage:

```
make
./construct_contours [flag] [number of datapoints]
```
You may enter one of the following flags, followed by the number of datapoints
to generate for the given shape:

1. -c for one circle
1. -2 for two circles
1. -3 for three circles
1. -d for a donut
1. -e for an ellipse
1. -p for a cardioid
1. -m for your own dataset (named \"my_data.dat\" and placed in the ./datapoints directory)
1. -s for a square
1. -t for a triangle

Please note: the GNUplot plotting utility must be installed for the path to be plotted.

## Brief:

This algorithm calculates the shortest path between a set of datapoints
using only curvature. I have named this method of distance computation
"Tao Distance", because it deals with calculating a constant defined
below:
```
Tao = dot_product(T1, T2)/(T1.length * T2.length)
```
Vectors T1 and T2 are unit tangent vectors calculated by the program. But since
the length of a unit vector is 1, we can write Tao as:
```
Tao = dot_product(T1, T2)
```
We then calculate curvature by using the length of the displacement
vector (T2 - T1) divided by the angle between T2 and T1. Given the previously
stated equivalence, Tao can be used to simplify the calculation of the angle between the two unit
tangent vectors. The final Tao-Distance Equation is defined below:
```
Tao-Distance = V.length + curvature + theta
```
Thus the program calculates this distance for each possibility,
chooses the point with the smallest Tao-Distance, and builds the
path without having to compare the found path with any others.

This program is designed to work for "convex datasets", meaning the
optimal solution will be a convex shape. This algorithm is less
accurate for shapes that have concave features, or for shapes that
contain points "inside" them.
The Time Complexity of this algorithm,
given that the dataset follows the above restrictions,
is O(n<sup>2</sup>).

## TODO:

My goal with this algorithm is to be able to reconstruct a planar, convex
shape given only its vertices. I have achieved this goal, as proven
by the various test cases (which are set as flags) for the program.
The next step is to write a MEL or C++ script as a plugin for
Autodesk Maya, implementing this algorithm to construct (or re-
construct) 3-dimensional meshes.

## Applications:

Some applications of this algorithm would be
the construction of 3-dimensional meshes from photographic references,
polarization of U/V components of objects (both polygons and NURBS),
and dynamic modeling of polygons at rendertime.

There are other applications not involving 3D Graphics;
because this algorithm calculates the shortest path between a set of
vertices that are known to be convex, it is a "solution"
to the Traveling Salesman Problem under very specific conditions.
My work for completing an algorithm that solves the TSP for
any set of panar Cartesian coordinates can be found in my project repo [pixarninja/shortest_path_algorithm](https://github.com/pixarninja/shortest_path_algorithm).
My work on the TSP uses concepts discovered through working on
the Contour Construction Algorithm, however a much different approach. My goal
with the Shortest Path Alogrithm is to be able to write
programs to complete all of the above applications, for any
3D pointcloud of data.
