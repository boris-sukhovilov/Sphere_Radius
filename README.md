# Sphere_Radius

adjointMatrix.m - Calculate adjoint matrix

calc_Radius.m - Calculate radius by all methods

CenterOfGravityCoordFromPairDistance.m - Calculate the coordinates of the vectors of n points in the coordinate system
                                         with the center at the center of gravity of the coordinates of the vectors of points
                                         and the axes directed along the eigenvectors, corresponding to the 3 largest eigenvalues
                                         matrix of scalar products of the coordinates of the vectors of points

Earth_radius.m - Determination Earth radius by flight distances

example_for_suh.m - (c) Vasiliev U.S example for Sukhovilov,
                    that it is undesirable to locate points near the poles
                    small measurement errors cause a large radius error?

final_rank.m - Calculation final rank of the measured matrix of half-squared distances

FitUsingLengths.m - Calculating the radius of a sphere using David Eberly method

generate_optim_tetrahedron_points.m - Generation of tetrahedron points with pairwise equal opposite edges 
                                      on a sphere of radius R

generate_sphere_points_on_isosceles_pyramid.m - Generate sphere points on isosceles pyramid

generateMatrixDistance.m - generate Matrix Distance

generateRandomPointsInSolidAngle.m - Generate Random Points In Solid Angle

getTetrahedronEdges.m - Get Tetrahedron Edges from Matrix Distance

MonteKarloTest.m - Perform Monte Carlo tests

optim_sphere.m - Determining the optimal configuration of "n" points on a sphere

optim_sphere_fun.m - Minimized functional for determining the optimal configuration 
                     of "n" points on a sphere

optimal_n_points_on_sphere.m, optimal_n_points_on_sphere_fix.m - Calculation of optimal 
                                                                 coordinates of n points on a sphere of radius r

paint.m - Display the surface of a sphere with points on it

pseudo_inv.m - Pseudo-inversion of a matrix using SVD

radius_optimal_n_points_on_sphere.m - Calculation sphere radius using formula (43)

regular_polyhedron.m - Calculation sphere radius for platonic solids

sensitivity_matrix.m - Calculation of the sensitivity matrix for the solution R from the implicit equation: 
                       f(R) = det(R^2*ones(n)-S.*S/2) = 0

sigma_eigenvalues.m - Calculating the standard deviation of the eigenvalues of the matrix of half-squared pairwise distances

SphereRadius_Carnot.m - Calculating the radius of a sphere circumscribing a tetrahedron using the Carnot method

SphereRadius_Cayley_Menger.m - Calculating the radius of a sphere circumscribing a tetrahedron using the Cayley-Menger determinants

SphereRadius_Euler_Grelle.m - Calculating the radius of a sphere circumscribing a tetrahedron using the Euler, Grelle method

SphereRadius_Sukhovilov.m - Calculating the radius of a sphere using Boris Sukhovilov method

SphereRadiusFromDistance.m - Function to call in Monte Carlo tests

SphereRadiusTest.m - Main function for experiments

sumith_fit.m - Calculating the radius of a sphere using Sumith YD method

IsFloatEqualRelative.m - Function for comparing (normalized) numbers with relative error maxRelDiff

kahanSumMatrix.m, kahanSumIterative.m - Kahan summation

## Citations

If you are interested in the scientific or methodological basis of the project, you can read the sources below:

Boris Sukhovilov [Determining Sphere Radius through Pairwise Distances](https://arxiv.org/abs/2504.02334)

```bibtex
@article{sukhovilov2025determiningsphereradiuspairwise,
         title={Determining Sphere Radius through Pairwise Distances}, 
         author={Boris Sukhovilov},
         year={2025},
         eprint={2504.02334},
         archivePrefix={arXiv},
         primaryClass={cs.CG},
         url={https://arxiv.org/abs/2504.02334}, 
}
