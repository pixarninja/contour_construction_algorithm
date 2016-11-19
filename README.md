Contour Construction Algorithm

BRIEF:

This algorithm calculates the shortest path between a set of data
using only curvature. I have named this method of distance computation
"Tao Distance", because it deals with calculating a constant defined
below:

		Tao = dot_product(T1, T2)/(T1.length * T2.length)

T1 and T2 are unit tangent vectors calculated by the program. But since
the length of a unit vector is 1, we can write Tao as:

		Tao = dot_product(T1, T2)

We then calculate curvature by using the length of the displacement
vector (T2 - T1) divided by the angle between T2 and T1. Tao
is used to simplify the calculation of the angle between the two unit
tangent vectors. The final Tao-Distance Equation is defined below:

		Tao-Distance = V.length + curvature + theta

Thus the program calculates this distance for each possibility,
chooses the point with the smallest Tao Distance, and builds the
path without having to compare the found path with any others.

This program is designed to work for "round datasets", meaning the
optimal solution will be a convex shape. This algorithm is less
accurate for shapes that have concave features, or for shapes that
contain points "inside" them.

TODO:

My goal with this algorithm is to be able to reconstruct a planar,
convex
shape given only its vertices. I have achieved this goal, as proven
by the various test cases (which are set as flags) for the program.
The next step is to write a MEL or C++ script as a plugin for
Autodesk Maya, implementing this algorithm to construct (or re-
construct) 3D polygons. Some applications of this would be
polarization of U/V components on objects (both polygons and NURBS),
or dynamic modeling of polygons at rendertime.
