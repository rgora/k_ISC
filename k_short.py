#!/usr/bin/env python
"""
This script calculates k_ISC in the short-time approximation from eq (11) of
JCP 134, 154105 (2011). It expects a few Turbomole v. 6.3 output files in two
directories for the initial singlet and final triplet states. Please put the
following files in each of the directories:

control
coord
vib_normal_modes
vibspectrum

The values of SOC and of the adiabatic energy gap E(S)-E(T) must be provided as
well.

Usage: k_short.py [options] singlet_dir triplet_dir

Options:
  -h, --help       show this help
  -d, --debug      print all calculated matrices
  -u, --unweighted calculate constants in dimension-less mass-unweighted
                   coordinates (by default the mass-weighted coordinates are used)
  -a, --aeg        adiabatic energy gap E(S)-E(T) in [au]
  -s, --soc        spin orbit coupling given as the sum of squares of matrix
                   elements for all three components of T state in [cm^-2]
"""

#     Copyright (C) 2015, Rafal Szabla (rafal.szabla@gmail.com)
#  
#     This program is free software; you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation; either version 2 of the License, or (at your
#     option) any later version.
#  
#     This program is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#     Public License for more details.
#  
#     You should have received a copy of the GNU General Public License along
#     with this program; if not, write to the Free Software Foundation, Inc.,
#     675 Mass Ave, Cambridge, MA 02139, USA.

__author__     = "Rafal Szabla and Robert W. Gora "
__maintainer__ = "Rafal Szabla"
__email__      = "rafal.szabla@gmail.com"
__license__    = "GPL"
__version__    = "1.0.2"

import sys, getopt, math
import numpy as np
import matplotlib.pylab as pl
from scipy.linalg import sqrtm, det
from math import sqrt, exp

global cmRectoHa
cmRectoHa = 4.5563352812122295e-06

def Usage():
    """Print usage information and exit."""
    print __doc__
    #print "Machine epsilon is: ",np.finfo(np.float64).eps,"for float64 type\n"

    sys.exit()

def Main(argv):
    '''Parse commandline and loop throught the logs'''

    DBG=False
    MWC=True

    # Parse commandline
    try:
        opts, args = getopt.getopt(argv, "hdua:s:", ["help",
                                                     "debug",
                                                     "unweighted",
                                                     "aeg",
                                                     "soc",
                                                     ])
    except getopt.GetoptError, error:
        print(error)
        Usage()
    if not argv:
        Usage()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            Usage()
        elif opt in ("-d", "--debug"):
            DBG=True
        elif opt in ("-u", "--unweighted"):
            MWC=False
        elif opt in ("-u", "--unweighted"):
            MWC=False
        elif opt in ("-a", "--aeg="):
            dE=float(arg)
        elif opt in ("-s", "--soc="):
            soc=float(arg)*cmRectoHa**2

    # Parse each data dir
    data_dirs = args

    # Get data for the initial singlet state
    R1, M1, nat = GetCoords(data_dirs[0])
    V1u, V1w = GetNormalModes(data_dirs[0], nat, DBG)
    O1 = GetFreq(data_dirs[0], nat)

    # Print debug information
    if DBG:
        print "Cartesian coordinates:\n", R1
        print "Atomic masses:\n", M1
        print "Mass-unweighted normal modes (dimensionless):\n", V1u
        print "Mass-weighted normal modes [sqrt(a.u.)-1]:\n", V1w
        print "Harmonic frequencies (Omega) [a.u.]:\n", O1
        print "Harmonic frequencies (Omega) [cm-1]:\n", O1/cmRectoHa

    # Get data for the final triplet state
    R2, M2, nat = GetCoords(data_dirs[1])
    V2u, V2w = GetNormalModes(data_dirs[1], nat, DBG)
    O2 = GetFreq(data_dirs[1], nat)

    # Test atomic masses vectors
    if M1.all() != M2.all():
        print "WARNING! Atomic masses vectors are different."

    # Print debug information
    if DBG:
        print "Cartesian coordinates:\n", R2
        print "Atomic masses:\n", M2
        print "Mass-unweighted normal modes (dimensionless):\n", V2u
        print "Mass-weighted normal modes [sqrt(a.u.)-1]:\n", V2w
        print "Harmonic frequencies (Omega) [a.u.]:\n", O2
        print "Harmonic frequencies (Omega) [cm-1]:\n", O2/cmRectoHa

    # Select mass-weighted or unweighted normal modes
    if MWC:
        V1=V1w
        V2=V2w
    else:
        V1=V1u
        V2=V2u

    # Duschinsky Matrix calculation
    # J=np.dot(np.transpose(V2),V1)
    J=V2u.T*V1u
    print "Duschinsky's matrix determinant is:", det(J)
    if DBG:
        print "Duschinsky matrix:\n", J.round(2)
        print "Plot of Duschinsky matrix (absolute values):"
        f = pl.figure()
        pl.imshow(np.absolute(J), interpolation='none', cmap=pl.cm.gray_r)
        pl.colorbar()
        pl.show()

    # Displacement vector
    # dR=np.subtract(R2,R1)
    # D=np.dot(np.dot(np.transpose(V2),M1),dR)
    dR=R1-R2
    if MWC:
        D=V2.T*sqrtm(M1)*dR
    else:
        D=V2.T*dR

    if DBG:
        print "Displacement vector R2-R1 [a.u.]:\n", dR
        print "Displacement in terms of normal coordinates D [a.u.]:\n", D

    # Evaluating the formula for short-time approximation
    k_au, k_s = k_short(J,D,O1,O2,soc,dE,nat,DBG)

def k_short(J,D,O1,O2,soc,dE,nat,DBG,linear=False):
    """Calculate ISC rate in short-time approx."""

    # conversion factor from a.u. to s-1
    au2s=2.4188843265052e-17
    n_coord, n_modes = DegFree(nat,linear)

    # matrices and constants
    # M=np.subtract(np.dot(np.dot(np.transpose(J),np.dot(O2,O2)),J),np.dot(O1,O1))		
    # A=np.dot(np.dot(np.transpose(J),np.dot(O2,O2)),D)
    # C=0.5*(np.dot(np.dot(np.transpose(D),np.dot(O2,O2)),D))
    M=J.T*(O2**2)*J-O1**2
    A=J.T*(O2**2)*D
    C=0.5*D.T*(O2**2)*D

    if DBG:
        print "Matrix M:\n", M
        print "Matrix A:\n", A
        print "C=%f\n" % C[0]

    # Alternative calculation of sum_a
    # (M.T*O1.I).trace()
    # np.sum(M.T*O1.I)
    sum_a=0
    for i in range(n_modes):
        sum_a+=M[i,i]/O1[i,i]
    
    sum_b=0
    for i in range(n_modes):
        for j in range(n_modes):
            sum_b+=(M[i,j])**2/(O1[i,i]*O1[j,j])

    sum_c=0
    for i in range(n_modes):
        sum_c+=(A[i])**2/O1[i,i]
    
    # k^{short}_{ISC} from Eq.(11) from JCP 134, 154105 (2011)
    k_sh=soc*sqrt(math.pi/(0.0625*sum_b + 0.25*sum_c))*exp(-1.0*((0.25*sum_a + C[0] - dE)**2)/(0.25*sum_b+sum_c))
    
    print "k_short = %s [a.u.]" % k_sh
    print "k_short = %e [1/s]"  % (k_sh/au2s)

    return k_sh, k_sh/au2s

def GetCoords(data_dir, linear=False):
    """Extract atomic coordinates."""

    # Conversion from AMU to a.u.
    amu2au = 1822.88853

    masses = {"c" : 12.01115, "o" : 15.99940, "n" : 14.00670, "h" : 1.00797, "s" : 32.066}
    atomic_masses=[]
    r=[]
    nat=0

    coord_list=open(data_dir+'/coord','r').readlines()

    # construct M matrix and read the atomic coordinates
    for c in coord_list[1:]:
        if c.startswith('$'):
            break
        temp_list=c.split()
    	atomic_masses.extend(3*[masses[temp_list[-1]]])
    	for x in temp_list[:-1]:
            r.append(float(x))
    	nat+=1

    print "There are %s atoms in %s." % (nat,data_dir)
    n_coord, n_modes = DegFree(nat,linear)

    # Return M (atomic masses) and R (atomic coordinates) vectors
    M=amu2au*np.diag(np.array(atomic_masses))
    R=np.array(r).reshape(n_coord,1)

    return np.matrix(R), np.matrix(M), nat

def DegFree(nat, linear=False):
    """Return no. of coords and normal modes"""
    n_coord = nat*3
    if not linear:
        n_modes = n_coord-6
    else:
        n_modes = n_coord-5

    return n_coord, n_modes

def GetNormalModes(data_dir, nat, DBG, linear=False):
    """Extract mass unweighted normal modes."""

    # Conversion from AMU to a.u.
    amu2au = 1822.88853

    # Read reduced masses from control
    cnt_file=open(data_dir+'/control','r')
    rdm=[]

    while 1:
        line=cnt_file.readline()
        if line=='':
            break
        if line.find('$vibrational reduced masses') !=-1:
            while 1:
                line=cnt_file.readline()
                if line=='' or line.startswith('$'):
                    break
                rdm.extend(line.split())

    Mr=amu2au*np.matrix(rdm[6:],dtype=np.float).reshape(1,-1)
    if DBG:
        print "Pertinent reduced masses [au]:\n",  Mr
        print "Pertinent reduced masses [amu]:\n", Mr/amu2au

    n_coord, n_modes = DegFree(nat,linear)

    print "Extracting mass unweighted vibrational modes from %s" % data_dir
    vib_list=open(data_dir+'/vib_normal_modes','r').readlines()
    vib_modes=[]

    for line in vib_list[1:]:
        for v in line.split()[2:]:
            vib_modes.append(float(v))

    # Return Vmuw (mass unweighted) and Vmw (mass-weighted) coordinates
    V=np.array(vib_modes).reshape(n_coord,n_modes+6)
    Vmuw=np.matrix(V[:,6:])
    Vmw=Vmuw/np.sqrt(Mr)

    if DBG:
        MrTest=1/np.diag(Vmuw.T*Vmuw)
        print "Pertinent reduced masses of mass-unweighted modes:\n",  MrTest

    # Test mass-weighted coordinates
    MrTest=1/np.diag(Vmw.T*Vmw)
    if MrTest.all() != Mr.all():
        print "WARNING! Mass weighted coordinates give different pertinent reduced masses"

    return Vmuw, Vmw

def GetFreq(data_dir, nat):
    """Extract vibrational frequencies."""

    vib_frequencies=[]

    print "Extracting frequencies from %s" % data_dir
    freq_list=open(data_dir+'/vibspectrum','r').readlines()

    for line in freq_list[9:-1]:
        f=line.split()[2]
        vib_frequencies.append(float(f)*cmRectoHa)

    # Return Omega matrix
    Omega=np.diag(np.array(vib_frequencies))

    return np.matrix(Omega)

#----------------------------------------------------------------------------
# Main routine
#----------------------------------------------------------------------------
if __name__ == "__main__":
    Main(sys.argv[1:])

