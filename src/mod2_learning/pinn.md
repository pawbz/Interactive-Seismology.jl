---
title: "To be updated"
order: 1
chapter: 1
section: 1
layout: "md.jlmd"
tags: ["module2"]
---


Many real-time optimization problems involve minimizing the cost or time taken to reach a goal. In this section,
we test the PINN method on the task of fnding the shortest-time path in several environments where an analytic
solution exists. One signifcant diference from the previous example (swinging up a pendulum) is that the time
taken is a variable and the subject to be minimized.
We present two famous examples of fnding the shortest-time path. (1) Fermat’s principle, or the principle of
least time, is that the path of the light ray between two given points is determined as the path that can be traveled
in the shortest time. (2) Another one is the shortest-time descent path between two given points under gravity,
known as the brachistochrone curve. It has been mathematically proven that it is part of a cycloid.
Te refraction of light can be explained by Fermat’s principle. When the refraction index in a medium varies
according to the location, the speed of the ray changes. In this case, the path of the ray bends to fnd a shortertime path in the medium. An example of the refraction of light within a medium where the refraction index
varies is shown in Fig. 3a.
Figure 3. Applications of PINN for determining the shortest-time path under specifc environments. (a)
Finding the shortest-time path of a light ray within the medium where the refraction index varies along y. (b)
Finding the shortest-time descent path between two given points under constant gravity. In the two fgures, the
analytic solutions are given in yellow dashed lines and the converged solutions by PINN are shown in blue. For
baselines, the results using RL with 105
 iterations are shown in magenta lines.
5
Vol.:(0123456789)
Scientifc Reports | (2024) 14:202 | https://doi.org/10.1038/s41598-023-49977-3
www.nature.com/scientificreports/
In the 2D xy space of Fig. 3a, the refraction index (n) has a sinusoidal profle with respect to y, as represented
by the contour color. Te ray slows down around y = 0.5 where the n is high. Terefore, the path of light that
passes in the shortest time will be determined to become more vertical near the slowest region, y = 0.5. Tus,
the ray trajectory will bend according to the changes in the refraction index, consistent with the analytic solution dictated by the law of refraction, as shown in yellow dashed line in Fig. 3a. When using RL to search for the
shortest-time arrival trajectory, even afer 105
 iterations of training, a trajectory far from the analytic solution
emerges, as indicated by the magenta line. While this trajectory does become vertical around y = 0.5 where n
is high, it is still far from being the shortest possible time.
To fnd the shortest-time path using PINN, the input of the neural network of Fig. 1 is set to the normalized
time (tN ), and the outputs are the x and y coordinates. The governing equation is given by
F =
 1
T
dx
dtN
2
+
 1
T
dy
dtN
2
−  c
n
2
, where c = 1 is the speed of light in vacuum and T is the time taken from
the starting point to the endpoint. Here, T is unknown and a trainable variable during the neural network training. Constraints are set as the coordinates of the starting and ending points, respectively (0, 0), (1, 1), and the
goal is set as Lgoal = T, aiming to minimize the time taken. More specifc settings can be found in the "Methods"
section.
Te solid lines in Fig. 3a depict the solution paths at each checkpoint (iteration I = 500, 1000, 1500, 2000),
with the blue solid line representing the converged fnal solution (I = 2232). Unlike traditional methods that
explore the shortest-time curve from a fxed start point to an endpoint, the PINN approach shows the characteristic of searching for the shortest-time curve while simultaneously approaching the start and endpoint. Te
converged solution shows a path that almost coincides with the analytic solution.
Figure 3b illustrates the search for the shortest-time descent curve (brachistochrone curve) connecting two
given points under gravity (g = 9.8 m/s2). Here, the yellow dashed line is the analytic solution for the shortesttime curve, which is a cycloid. In the RL solution with 105
 iterations shown in magenta, the path more vertically
free-falls to increase velocity and then changes direction horizontally at the lower altitude to move at a high
speed. While this seems intuitively plausible, it deviates from the actual shortest curve, the cycloid. On the other
hand, for PINN, the governing equation is set as the mechanical energy conservation law,
F = gy0 −

gy + 1
2
 1
T
dx
dtN
2
+
 1
T
dy
dtN
2
, and the constraints were set using the start ((x0, y0) = (0, 1))
and endpoint coordinates ((x1, y1) = (1, 0)). Te goal, similar to Fig. 3a, was set as Lgoal = T to fnd the shortesttime curve. In Fig. 3b, at the initial phase (I = 500), a non-physical path that fails to align the start and endpoints
is derived. But it gradually converges to a path close to the analytic solution shown in yellow. Tis approach using
PINN to search for the shortest-time path can be utilized in chip design to fnd the minimum-loss circuit.
