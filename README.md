# LearningToIntegrate

The repository contains three folders which we briefly introduce.

## FM_Code

Our python code for learning the transport maps (via ACF, CFM and OTCFM) and the transformation of the quadrature points. As input for the learning process, modes of random samples of both Gaussian and Lévy fields must be provided. Afterwards, the nodes of a sparse Gauss-Hermite quadrature rule can be transformed to quadrature nodes for the target distribution.

Approximation degree A2 = 9 modes
Approximation degree A3 = 25 modes

## Simulation

Our R script for the generation of sparse grids as well as generating and interpolating of random fields. It provides functions that can generate realizations of generalized random fields and convolve them with Matern kernels for various approximation degrees (including non-truncated fields). The mode vector (needed for learning in using the FM code) can be exatracted or the fields can be interpolated to the FEM meshes (evaluated at the center of mass). The random fields can also be interpolated given a list of modes (e.g. transformed quadrature points). The code accepts a list containing the following columns: ID (1 column), entries of the mode vector (9 for A2 or 25 columns for A3), SG level + weights (1 column each).

## PDESolver

Matlab code that computes the solution of the Darcy flow equation with given realizations of random fields. As a input, a table containing the values of the random fields interpolated to the centers of mass of the mesh is needed. When using this code make sure to:
- add the curlcurl directory to the working space, e.g. with addpath('curlcurl')
- select the correct mesh in the EvaluateQuadNodes.m script
- check the data path in the import_realizations function of rf.MaternByCell
- check the data path in the export_results function in uq.EvaluatePredefinedSamples
- check if the number of samples (nsample) in the solve function in uq.EvaluatePredefinedSamples matches the number of realizations
