k_short.py
==========

The k_short.py calculates intersystem crossing rates (k_ISC) in the short-time
approximation following Eq (11) of JCP 134, 154105 (2011).

It expects a few Turbomole 6.3 output files in two directories for the initial
singlet and the final triplet states.

Please put the following files in each of the directories:

control
coord
vib_normal_modes
vibspectrum

The values of SOC (in here the sum of squares of absolute magnitudes of SOC
matrix elements between S and all three components of T state) and of the
adiabatic energy gap E(S)-E(T) must be provided as well.

## Usage ##

Usage: k_short.py [options] singlet_dir triplet_dir

Options:
  -h, --help       show this help
  -d, --debug      print all calculated matrices
  -u, --unweighted calculate constants in dimension-less mass-unweighted
                   coordinates (by default the mass-weighted coordinates are used)
  -a, --aeg        adiabatic energy gap E(S)-E(T) in [au]
  -s, --soc        spin orbit coupling given as the sum of squares of matrix
                   elements for all three components of T state in [cm^-2]
