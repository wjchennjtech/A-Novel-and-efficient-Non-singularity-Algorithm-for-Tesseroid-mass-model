This script accompanies the paper “A Novel and Efficient Non-Singularity Algorithm for the Tesseroid Mass Model.”

Operating System: Windows 10/11, MATLAB2024

Input data format:
Column 1: Longitude
Column 2: Latitude
Column 3: Boundary (unit: meters) or Density (unit: kilograms per cubic meter)


Adjustable parameters:
nmax: Maximum degree of the spherical harmonic expansion
Ntheta: Number of Gaussian–Legendre quadrature points for colatitude integration
height: Computation height

Code Description：
This program consists of one main script and four subroutines:
GaussLegendre.m: Gaussian–Legendre quadrature function for computing the associated Legendre integrals with respect to colatitude.
plm.m: Associated Legendre function.
SHA_single_tess.m: Spherical harmonic expansion function for a single tesseroid
SHS_regular.m: Spherical harmonic synthesis function.
The main script is named Main_script.m, which calls all subroutines to perform the computations described in the paper.

