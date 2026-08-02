# QENS_utilities
## A collection of python functions that may be useful as building blocks for QENS data analysis

**Overview** 

The file *utilities_qens.py* contains a collection of small functions that may help in writing scripts for the analysis of quasi-elastic neutron scattering (QENS) data. These small functions are designed to be able to employ *curve_fit* from *scipy.optimize* for the data fitting.

The functions have been inspired by the analysis of QENS data from proteins in solution, the concept of which is explained in
M. Grimaldo et al., EPJ Web of Conferences 83, 02005 (2015) published by the European Physical Society;
https://dx.doi.org/10.1051/epjconf/20158302005 .
Further general information may for instance be found in M. Grimaldo et al., Quarterly Reviews of Biophysics 52, e7 (2019) published by Cambridge University Press;
https://dx.doi.org/10.1017/S0033583519000027 . Additional details are described in C. Beck et al., EPJ Web of Conferences 272, 01004 (2022); https://doi.org/10.1051/epjconf/202227201004 .
The article I. Mosca et al., The Journal of Physical Chemistry Letters 15, 12051 (2024), https://doi.org/10.1021/acs.jpclett.4c02533 , contains an example of application.

**Disclaimer**

Although the concepts used for the functions within _utilities_qens.py_ are well established (see above), the functions themselves constitute new implementations that have not yet been fully tested. They are made available for the sole purpose of an exchange of ideas.

**Documentation**

Central to this approach to QENS data analysis is the description of the spectrometer energy resolution function _R_ by a sum of Gaussian functions such that this resolution function can be accounted for analytically - as opposed to a numerical convolution - when fitting a model _S_ to the observed spectra. This model scattering function _S_ can consist of a sum of an arbitrary number of Lorentzian functions that account for diffusion processes, of an elastic scattering contribution modeled by a Dirac function, as well as of an affine background defined by a slope and a constant offset. The Lorentzian and Dirac functions are centered at zero energy transfer by the definition of quasi-elastic scattering. The Lorentzian and Dirac functions are transformed into sums of Voigt functions by the analytical convolution with _R_.

QENS spectra are generally contained in 2-dimensional matrices _M_ storing the scattering intensity for different momentum and energy transfers. Typically, there are on the order of 20 momentum transfers _q_ and 1e3 energy transfers _w_. 

The collection of functions contains a wrapper named *wrapper_fit_func_parser* that can be used to interface with *scipy.optimize.curve_fit*. Using this wrapper, either vectors containing a slice along _w_ at fixed _q_, or the entire matrix _M_ rearranged into a single long vector can be passed to *curve_fit*. In the latter case, the function _array2vector_ is used to rearrange _M_. When _M_ is passed subsequent to this rearrangement, a global fit of the model function for all energy and momentum transfers at once can be performed. This global fit is used in order to impose the dependence on the momentum transfer of the diffusion processes described by the model function. (More recent versions of *curve_fit* can take matrices as arguments directly, circumventing the need for _array2vector_.)

A typical call of the fit function can be written as follows:

_popt, pcov = curve_fit( lambda x, *p: wrapper_fit_func_parser( x, q, n, r, len( f0 ), p, **protein_fit ), x, y, p0=f0, bounds=(l,u), sigma=dy )_

Therein, the vectors _x_, _y_, and _dy_ contain the energies, corresponding scattering intensities, and their errors, respectively. _q_ is the vector containing the absolute values of the momentum transfer(s), _n = range(...,...)_ contains the indices of the momentum transfers within _q_ to be taken into account for the fit, and _f0_ the vector for the initial guess for the fit parameters _p_.
_r_ is the matrix containing the parameters of the Gaussians describing the spectrometer resolution function.

Using a dictionary of keyword arguments, denoted _**protein_fit_ in the above example, the model for the fit can be chosen, and information required by that specific model can be passed to the model function, such as a fixed parameter.

Specific models are selected by their name assigned within the function *model_sqw_parser*. This name is passed via the keyword *SQWmodel*.
For any model chosen, additional keywords have to specify the associated number of Lorentzians as well as whether or not a Dirac, a sloped, or a flat background are to be used.

Depending on the specific model, additional keywords may apply that for instance pass fixed parameters.

For example, 

_protein_fit = dict( SolventIntensities=Isolv, SolventWidths=Wsolv, NumberLorentzians=4, UseDirac=1, UseSlopedBackground=1, 
				       UseFlatBackground=1, BGSlope=bgs, BGFlat=bgc, SolventDirac=dsolv, SQWModel='BATSGlobal2State' )_
               
selects a model performing a global fit for several momentum transfers at once that contains the following components, as explained in 
M. Grimaldo et al., Phys.Chem.Chem.Phys. 17, 4645  (2015), https://dx.doi.org/doi:10.1039/C4CP04944F, equation 6:

- One Lorentzian that accounts for the center-of-mass diffusion of tracer proteins in an aqueous solution sample and obeys _gamma = D q<sup>2</sup>_, where _gamma_ is the Lorentzian half width at half maximum;

- Two Lorentzians that account for the internal diffusive dynamics. The amplitudes and widths of these Lorentzians are coupled as explained in the above reference, based on F. Roosen-Runge et al., J. Chem. Phys. 144, 204109 (2016), https://doi.org/10.1063/1.4950889. This coupling is implemented in the function two_state;

- A fourth Lorentzian that accounts for the solvent water (D<sub>2</sub>O) contribution;

- A Dirac accounting for the scattering signal from the sample container;

- An affine background defined by slope and offset. 

The solvent Lorentzian intensities and widths as well as the elastic scattering intensity are fixed parameters that are passed to the model function as vectors via the keywords _SolventIntensities_, _SolventWidths_, and _SolventDirac_, respectively.

Importantly, the entries in the vectors for the initial guess _p_ and bounds _(l, u)_ as well as the lengths of these vectors depend on the chosen model. In the above example, the first 5 entries account for the global parameters that apply for all _q_, i.e. the center-of-mass diffusion coefficient _D_, coupled internal diffusion coefficients _D<sub>1</sub>_, _D<sub>2</sub>_, and associated residence times _tau<sub>1</sub>_, and _tau<sub>2</sub>_. The following entries account for the _q_-dependent parameters _beta(q)_ and _A<sub>0</sub>(q)_, where _beta_ is a free amplitude scaling parameter and 0 <= _A<sub>0</sub>(q)_ <= 1 the elastic incoherent structure factor. There are _n_ entries for beta(q) followed by _n_ entries for _A<sub>0</sub>(q)_ in the parameter vectors.

In the models implemented in *model_sqw_parser*, _p_ represent the free parameters of the model scattering function _S_, and _fp_ the associated fit parameters. The equations _p[ ...] = ... fp[...] ..._ connect the fixed and fit parameters to these parameters _p_ of _S_.
For global fits for all _q_ at once, _p[k,...]_ becomes two-dimensional, the first dimension representing the index of _q_. In the models, the order of the components is fixed: First, all Lorentzians, subsequently the Dirac, and finally the sloped and flat background.

The keyword arguments are passed via the function *wrapper_fit_func_parser* onwards to several functions that use the arguments.  

Plots of the fit results and their components can be generated by calling

_fitres = model_sqw_subplots_parser( x, popt, q, n, r, **protein_fit )_ ,

and the goodnes of the fit can be calculated by

_gof = goodness_of_fit( x, y, popt, q, n, r, dy, **protein_fit )_ .

The function *convoluted_model_parser* builds the actual fit function from a sum of Voigt functions and the background. It is called by *model_sqw_parser* which contains the different models for the scattering function. Importantly, *model_sqw_parser* can be extended by adding new models following the existing scheme.

The spectrometer resolution function expressed by a sum of an arbitrary number of Gaussians can be fitted by the function *resolution_function*. This function is quite involved, as it has been attempted to accommodate several different spectrometers. A successful fit of a resolution function by a sum of Gaussians in general requires a considerable effort in adapting initial guesses and boundaries. It is noted that the Gaussians may have non-zero center positions along the energy axis.

The other functions in the collection, such as load routines for data, are rather self-explaining.

It is noted that the low-level interfacing to *curve_fit* can be simplified by employing the python package *lmfit*.

**Module to fit the Mittag-Leffler function**

The *utilities_qens.py* also allow to fit the Fourier transform of the Mittag-Leffler function (MLF), which can be seen as a generalized Lorentzian function. 
In this case, in the absence of a "generalized Voigt" function, a numerical convolution with the resolution function is required.
This numerical convolution is only applied to the MLF part of the model, whereas the other components continue to be convoluted analytically.
The numerical convolution uses the description of the resolution function by the sum of Gaussians, as does the analytical convolution, and performs the convolution
on a longer energy axis than the data axis to reduce numerical truncation effects.
The Mittag-Leffler fit can be included as illustrated by the lines below, where the first line allocates a vector with the length of the number of used _q_-values. 
This vector _Lw_ will be used to store the _q_-dependent widths of the Mittag-Leffler function.

_Lw = np.zeros( np.shape( qp ) )_

_protein_fit = dict( SolventIntensities=Isolv, SolventWidths=Wsolv, NumberLorentzians=2, UseDirac=1, UseSlopedBackground=1, UseMLF=1, GammaL=GLw,
		       UseFlatBackground=1, BGSlope=bgs, BGFlat=bgc, SolventDirac=dsolv, SQWModel='BATSGlobalMLF' )_

Within the model implemented inside *model_sqw_parser*, the MLF has to follow immediately after the Lorentzians and precede the Dirac and sloped and flat backgrounds, as illustrated by the example model _BATSGlobalMLF_.
The MLF can be fitted for all _q_ at once with a global parameter _alpha_, such that a global fractional diffusion coefficient _gamma<sup>alpha</sup> = D<sub>alpha</sub> q<sup>2</sup>_ as a solution of the fractional diffusion equation can be fitted.
For further reading see, e.g., I. Krasnov et al., J. R. Soc. Interface 13, 20160506 (2016), https://dx.doi.org/10.1098/rsif.2016.0506 , and references therein.

