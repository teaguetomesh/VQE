'''
Teague Tomesh - 3/24/2019

Use the OpenFermion package to calculate the actual energies for molecules
using the Full Configuration Interaction (FCI) method.

'''

from openfermion.hamiltonians import MolecularData

# Set molecule parameters.
basis = 'sto-3g'
multiplicity = 1
bond_length_interval = 0.1
n_points = 30
name = 'H2'

# Generate molecule at different bond lengths.
fci_energies = []
bond_lengths = []
for point in range(1, n_points + 1):
    bond_length = bond_length_interval * point
    bond_lengths += [bond_length]
    description = str(round(bond_length,2))
    print(description)
    geometry = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length))]
    molecule = MolecularData(geometry, basis, multiplicity, description=description)
    
    # Load data.
    molecule.load()

    # Print out some results of calculation.
    print('\nAt bond length of {} angstrom, molecular hydrogen has:'.format(
        bond_length))
    print('FCI energy of {} Hartree.'.format(molecule.fci_energy))
    fci_energies += [molecule.fci_energy]

# Save data to file
fname = 'FCI_Energies/'+name+'_'+basis+'_FCI_energy.txt'
with open(fname, 'w') as fn:
    for r, e in zip(bond_lengths, fci_energies):
        fn.write('{0:.1f} {1:.7f}\n'.format(r,e))



