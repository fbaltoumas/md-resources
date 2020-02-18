# md-resources
Some resources and analysis scripts for Molecular Dynamics simulations

## Force Field Topologies and Parameters
#### NAMD
- [`MARTINI_CHOLESTEROL_NAMD.zip`](MARTINI_CHOLESTEROL_NAMD.zip): NAMD-formatted topology and parameter files for Coarse-Grained cholesterol, based on the original MARTINI model (Marrink et al, J. Phys. Chem. B 2007, 111 (27),  7812-7824).
- [`MARTINI_PIs_NAMD.zip`](MARTINI_PIs_NAMD.zip): NAMD-formatted files for three PI lipids, intended for use with the MARTINI and PACE force fields.  Based on the equivalent GROMACS parameters, as presented in Lopez et al, JCTC, 9:1694–1708, 2013.

#### GROMACS
- [`Divalent_Metal_ions_GMX..zip`](Divalent_Metal_ions_GMX..zip): Parameters for divalent metal ions (Fe2+, Mn2+, Mg2+, Be2+ etc) in GROMACS format for use with the TIP3P (CHARMM, AMBER) and SPC / SPCE (GROMOS) water models.  Based on the work by Li et al, J Chem Theory Comput 2013, 9(6):2733–2748.

## Analysis Scripts
#### Python
- [`GMX Jarzynski (Python 3.x)`](gmx_jarzynski_py3.py): A Python 3.x script to perform Potential of Mean Force (PMF) calculations using Jarzynski's equality, for pulling simulations generated by GROMACS.
- [`GMX Jarzynski (Python 2.7)`](gmx_jarzynski_py2.7.py): The same script but for the (soon to be obsolete) Python 2.7 version.

#### Perl
- [`Karplus Equation Calculator`](karplus_couplings.pl): A Perl script that uses the Karplus equation to predict J-couplings based on dihedral angle calculations from MD simulation data, as produced by the GROMACS "gmx angle" utility.
