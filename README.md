# QENS_utilities
A collection of python functions that may be useful as building blocks for QENS data analysis

The file utilities_qens.py contains a collection of small functions that may help in writing scripts for the analysis of quasi-elastic neutron scattering (QENS) data. These small functions are designed to be able to employ curve_fit from scipy.optimize for the data fitting.

The functions have been inspired by the analysis of QENS data from proteins in solution, the concept of which is explained in
M. Grimaldo et al., EPJ Web of Conferences 83, 02005 (2015) published by the European Physical Society;
https://dx.doi.org/10.1051/epjconf/20158302005

Central to this approach to QENS data analysis is the description of the spectrometer energy resolution function R by a sum of Gaussian functions such that this resolution function can be accounted for analytically - as opposed to a numerical convolution - when fitting a model S to the observed spectra. This model scattering function S can consist of a sum of an arbitrary number of Lorentzian functions that account for diffusion processes, of an elastic scattering contribution modeled by a Dirac function, as well as of an affine background defined by a slope and a constant offset. The Lorentzian and Dirac functions are centered at zero energy transfer by the definition of quasi-elastic scattering. The Lorentzian and Dirac functions are transformed into sums of Voigt functions by the analytical convolution with R.

QENS spectra are generally contained in 2-dimensional matrices M storing the scattering intensity for different momentum and energy transfers. Typically, there are on the order of 20 momentum transfers q and 1e3 energy transfers w. 

The collection of functions contains a wrapper named wrapper_fit_func_parser that can be used to interface with scipy.optimize.curve_fit. Using this wrapper, either vectors containing a slice along w at fixed q, or the entire matrix M rearranged into a single long vector can be passed to curve_fit. In the latter case, the function array2vector is used to rearrange M. When M is passed subsequent to this rearrangement, a global fit of the model function for all energy and momentum transfers at once can be performed. This global fit is used in order to impose the dependence on the momentum transfer of the diffusion processes described by the model function.

A typical call of the fit function can be written as follows:

popt, pcov = curve_fit( lambda x, *p: wrapper_fit_func_parser( x, q, n, r, len( f0 ), p, **protein_fit ), x, y, p0=f0, bounds=(l,u), sigma=dy )

Therein, the vectors x, y, and dy contain the energies, corresponding scattering intensities, and their errors, respectively. q is the vector containing the absolute values of the momentum transfer(s), n = range(...,...) contains the indices of the momentum transfers within q to be taken into account for the fit, and f0 the vector for the initial guess for the fit parameters p.
r is the matrix containing the parameters of the Gaussians describing the spectrometer resolution function.

Using a dictionary of keyword arguments, denoted **protein_fit in the above example, the model for the fit can be chosen, and information required by that specific model can be passed to the model function, such as a fixed parameters.

Specific models are selected by their name assigned within the function model_sqw_parser. This name is passed via the keyword SQWmodel.
For any model chosen, additional keywords have to specify the associated number of Lorentzians as well as whether or not a Dirac, a sloped, or a flat background are to be used.

Depending on the specific model, additional keywords may apply that for instance pass fixed parameters.

For instance, 

protein_fit = dict( SolventIntensities=Isolv, SolventWidths=Wsolv, NumberLorentzians=4, UseDirac=1, UseSlopedBackground=1, 
				       UseFlatBackground=1, BGSlope=bgs, BGFlat=bgc, SolventDirac=dsolv, SQWModel='BATSGlobal2State' )
               
selects a model that performs a global fit for several momentum transfers at once that contains the following components, as explained in 
M. Grimaldo et al., Phys.Chem.Chem.Phys. 17, 4645  (2015), https://dx.doi.org/doi:10.1039/C4CP04944F, equation 6:

- One Lorentzian that accounts for the center-of-mass diffusion of tracer proteins in an aqueous solution sample and obeys gamma = D q^2, where gamma is the Lorentzian half width at half maximum;

- Two Lorentzians that account for the internal diffusive dynamics. The amplitudes and widths of these Lorentzians are coupled as explained in the above reference, based on F. Roosen-Runge et al., J. Chem. Phys. 144, 204109 (2016), https://doi.org/10.1063/1.4950889. This coupling is implemented in the function two_state;

- A fourth Lorentzian that accounts for the solvent water (D2O) contribution;

- A Dirac accounting for the scattering signal from the sample container;

- An affine background defined by slope and offset. 

The solvent Lorentzian intensities and widths as well as the elastic scattering intensity are fixed parameters that are passed to the model function as vectors via the keywords SolventIntensities, SolventWidths, and SolventDirac, respectively.

Importantly, the entries in the vectors for the initial guess p and bounds (l, u) as well as the lengths of these vectors depend on the chosen model. In the above example, the first 5 entries account for the global parameters that apply for all q, i.e. the center-of-mass diffusion coefficient D, coupled internal diffusion coefficients D1, D2, and associated residence times tau1, and tau2. The following entries account for the q-dependent parameters beta(q) and A0(q), where beta is a free amplitude scaling parameter and 0<=A0(q)<=1 the elastic incoherent structure factor. There are n entries for beta(q) followed by n entries for A0(q) in the parameter vectors.

The keyword arguments as passed via the function wrapper_fit_func_parser onwards to several functions that use the arguments.  

Plots of the fit results and their components can be generated by calling

fitres = model_sqw_subplots_parser( x, popt, q, n, r, **protein_fit ) ,

and the goodnes of the fit can be calculated by

gof = goodness_of_fit( x, y, popt, q, n, r, dy, **protein_fit ) .





