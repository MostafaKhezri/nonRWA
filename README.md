# nonRWA level crossing

This is a mathematica package that calulates level crossings and effective couplings for a trasnmon coupled to a resoantor. The paper detailing the Physics of this can be find at [arXiv:1606.05721](http://arxiv.org/abs/1606.05721).

## System
The system consists of a [transmon](https://arxiv.org/abs/cond-mat/0703002) coupled to a readout resonator through a charge-charge interaction. The qubit state can be [read](http://arxiv.org/abs/1401.0257) by populating the readout resonator and measureing amplitude/phase shift of the leaking field of the resoantor.

Here I (and therefore this code) assume that the eigenenergies of the system do not change _much_ doe to non-RWA terms. This is a very good assumption for the typical parameters in the circuit QED measurements of transmons, and I have also checked it separately (see [Future improvement](#Sorting-nonRWA-Hamiltonian-eigenstates/eigenstates))

## Goal
Find the qubit levels, as well as the phton numbers at which level crossing between qubit-resonator levels occures.

## The code
The package consists of functions that do the following:
1. Diagonalize the qubit-resonator system in the RWA case, and correctly sort (label) the resulting eigenstate/eigenenergy pairs.
2. Find the energy resonance conditions, to get the photon numbers at which the level crossing occurs.
3. Calculate the effective coupling between the resonant levels, due to non-RWA terms.

### Disclaimer
- This package solves a specific set of Physics problems, and is written based on some Physical assumptions and in some parameter regime. The output of the functions in this package may be mathematically correct, but it may not be Physically meaningfull for arbitrary parameter regimes.
- This package is not written to produce error messages when the functions are not used properly.
- Although I did my best to have all the different functions as independant and modular as possible, due to the nature of the problem that is being solved, this was only doable up to a limit.
- Following the previous point, this means that when calling functions that use the output and/or parameters of some previous functions, things must be consistent.

## Usage
The repo includes a Mathematica notebook named [test.nb](./test.nb), which contains examples of all the functions of the package.

In general, to get a complete manula for each function, run (in notebook or console of Mathematica) `?FunctionName`

## Known issues
##### Failing at large n
When number of photon increases above some value (this value depends on parameters), the diagonalizer fails. This issue does not depen on how powerfull your machine is. This happens due to the way that I am sorting (labeling) the eigenvalues, and there is a physical reason for this. I may know a way to fix this, see the [future idea section](#Different-digonalizing/sorting).

## Future improvements

##### Different digonalizing/sorting
The system can be diagonalized in block matrices, and then all these matrices can be combined toghether to form the total system. This is not expected to make things faster, but it most probably prevents the issue of failing to sort the eigenenergies at large photon numbers.

##### Sorting nonRWA Hamiltonian eigenstates/eigenstates
This code diagonalizes and sorts the RWA Hamiltonian of the system. Although for the typical parameters in the circuit QED measurements of transmons this is a very good approximation, but one can also calculate the non-RWA eigenenergies and eigenstates. However, one can no longer use the same [method](https://arxiv.org/abs/1606.04204) to sort these non-RWA eigenenergies. There exists a way to sort these new eigenenergies, and I have a (not-cleaned-up) code that does this. I will add that as a new function to this package in the future.
