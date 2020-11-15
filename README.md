# QENS_utilities
A collection of python functions that may be useful as building blocks for QENS data analysis

The file utitlities_qens.py contains a collection of small functions that may help in writing scripts for the analysis of quasi-elastic neutron scattering (QENS) data.

These functions have been inspired by the analysis of QENS data from proteins in solution, the concept of which is explained in
M. Grimaldo et al., EPJ Web of Conferences 83, 02005 (2015) published by the European Physical Society;
https://dx.doi.org/10.1051/epjconf/20158302005

Central to this approach to QENS data analysis is the description of the spectrometer energy resolution function R by a sum of Gaussian functions such that this resolution function can be accounted for analytically - as opposed to a numerical convolution - when fitting a model S to the observed spectra. This model scattering function S can consist of a sum of an aribrary number of Lorentzian functions that account for diffusion processes, of an elastic scattering contribution modeled by a Dirac function, as well as of an affine background defined by a slope and a constant offset. The Lorentzian and Dirac funtions are centered at zero energy transfer by the definition of quasi-elastic scattering. The Lorentzian and Dirac functions are transformed into sums of Voigt functions by the analytical convolution with R.

QENS spectra are generally contained in 2-dimensional matrices M storing the scattering intensity for different momentum and energy transfers. Typically, there are on the order of 20 momentum transfers q and 1e3 energy transfers w. 

The collection of functions contains a wrapper named wrapper_fit_func_parser that can be used to interface with scipy.optimize.curve_fit. Using this wrapper, either vectors containing a slice along w at fixed q, or the entire matrix M rearranged into a single long vector can be passed to curve_fit. In the latter case, the function array2vector is used to rearrange M.

A typical call of the fit function can be written as follows:

popt, pcov = curve_fit( lambda x, *p: wrapper_fit_func_parser( x, q, n, r, len( f0 ), p, **protein_fit ), x, y, p0=f0, bounds=(l,u), sigma=dy )

Therein, the vectors x and y contain the energies and corresponding scattering intensities. q is the vector containing the momentum transfer(s), n = range(...,...) contains the indices of the momentum transfers within q to be taken into account for the fit, and f0 the initial guess.


