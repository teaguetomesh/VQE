'''
Teague Tomesh - 3/24/2019

Use the OpenFermion package to calculate the actual energies for molecules
using the Full Configuration Interaction (FCI) method.

'''

from openfermion.hamiltonians import MolecularData
import glob
import sys

# set molecule information to create the proper output filename
name = 'H2'
basis = 'sto-3g'

# glob all of the molecule files
molecule_files = glob.glob('molecule_data/*sto-3g_singlet*.hdf5')
molecule_files = sorted(molecule_files)

fci_energies = []
bond_lengths = []
for fn in molecule_files:
    # Load the current molecule
    molecule = MolecularData(filename=fn)

    # the bond length should be the description
    bond_lengths += [molecule.description]
    print(molecule.description)
    
    # Print out some results of calculation.
    print('\nAt bond length of {} angstrom, molecular hydrogen has:'.format(
        molecule.description))
    print('FCI energy of {} Hartree.'.format(molecule.fci_energy))
    fci_energies += [float(molecule.fci_energy)]

print(bond_lengths)
print(fci_energies)

# Save data to file
fname = 'FCI_Energies/{}_{}_FCI_energy.txt'.format(name,basis)
with open(fname, 'w') as fn:
    for r, e in zip(bond_lengths, fci_energies):
        fn.write('{0} {1:.7f}\n'.format(r,e))



