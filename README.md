# dubins3D
A Dubins path is the shortest length path for an object with a bounded curvature (minimum turning radius).  Our ICRA 2024 paper shows how to generate all the Curve-Straight-Curve paths between two starting and ending positions and orientations in 3D.  We show a one-to-one mapping between a RRPRR robot and a 3D CSC Dubins path.  The solution involves forming and solving a system of multivariate polynomials in a way that is similar to Gauss elimination but for nonlinear equations – techniques deeply rooted in resultant-based elimination methods pioneered by Newton and Euler, and later popularized mostly by Sylvester. The goal is to eliminate all but one variable resulting in a polynomial whose roots can be substituted back into the system to solve for the remaining variables.

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/aWfmgsal0JU/0.jpg)](https://www.youtube.com/watch?v=aWfmgsal0JU)

[ICRA, May 14 2024](https://2024.ieee-icra.org/) "An Analytic Solution to the 3D CSC Dubins Path Problem" by Victor M. Baez, Nikhil Navkar, Aaron T. Becker.


This work was initiated at the [Dagstuhl Seminar 23091, "Algorithmic Foundations of Programmable Matter" in 2023](https://www.dagstuhl.de/23091). We thank [Cynthia Sung](https://sung.seas.upenn.edu/research/) for presenting this open problem and for discussions on practical applications.  
Support by the National Priority Research Program (NPRP) award (NPRP13S-0116-200084) from the Qatar National Research Fund (a member of The Qatar Foundation), 
the [Alexander von Humboldt Foundation](https://www.humboldt-foundation.de/en/), and the National Science Foundation under
[IIS-1553063](http://nsf.gov/awardsearch/showAward?AWD_ID=1553063),
[1849303](https://nsf.gov/awardsearch/showAward?AWD_ID=1849303), and
[2130793](https://nsf.gov/awardsearch/showAward?AWD_ID=2130793).

>Paper abstract: We present an analytic solution to the 3D Dubins path problem for paths composed of an initial circular arc, a straight component, and a final circular arc. These are commonly called CSC paths.  By modeling the start and goal configurations of the path as the base frame and final frame of an RRPRR manipulator, we treat this as an inverse kinematics problem. The kinematic features of the 3D Dubins path are built into the constraints of our manipulator model. 
Furthermore, we show that the number of solutions is not constant, with up to seven valid CSC path solutions even in non-singular regions.

For programming a 2D Dubins path, see

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/tjUVE3jto98/0.jpg)](https://www.youtube.com/watch?v=tjUVE3jto98)

For programming a 2D Reeds Shepp path, see

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/fAqh_cy7ePI/0.jpg)](https://www.youtube.com/watch?v=fAqh_cy7ePI)

