from openfermion.hamiltonians import MolecularData
from openfermion.transforms import get_fermion_operator, jordan_wigner, bravyi_kitaev

# This code taken from the example tutorial at:
# https://github.com/quantumlib/OpenFermion/blob/master/examples/openfermion_tutorial.ipynb

# I have added my own notes and code throughout

'''
# Set parameters to make a simple molecule.
diatomic_bond_length = .7414
geometry = [('H', (0., 0., 0.)), ('H', (0., 0., diatomic_bond_length))]
basis = 'sto-3g'
multiplicity = 1
charge = 0
description = str(diatomic_bond_length)

# Make molecule and print out a few interesting facts about it.
molecule = MolecularData(geometry, basis, multiplicity,
					 charge, description)
print('Molecule has automatically generated name {}'.format(
	    molecule.name))
print('Information about this molecule would be saved at:\n{}\n'.format(
	    molecule.filename))
print('This molecule has {} atoms and {} electrons.'.format(
	    molecule.n_atoms, molecule.n_electrons))
for atom, atomic_number in zip(molecule.atoms, molecule.protons):
	    print('Contains {} atom, which has {} protons.'.format(
			    atom, atomic_number))
'''

# Can also generate molecule at different bond distances
# Set molecule parameters.
basis = 'sto-3g'
multiplicity = 1
bond_length_interval = 0.1
n_points = 25

# Generate molecule at different bond lengths.
hf_energies = []
fci_energies = []
bond_lengths = []
for point in range(3, n_points + 1):
	bond_length = bond_length_interval * point
  bond_lengths += [bond_length]
  description = str(round(bond_length,2))
  print(description)
  geometry = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length))]
  molecule = MolecularData(geometry, basis, multiplicity, description=description)

  # Load data.
  molecule.load()

  # Print out some results of calculation.
  print('\nAt bond length of {} angstrom, molecular hydrogen has:'.format(bond_length))
  #print('Hartree-Fock energy of {} Hartree.'.format(molecule.hf_energy))
  #print('MP2 energy of {} Hartree.'.format(molecule.mp2_energy))
  #print('FCI energy of {} Hartree.'.format(molecule.fci_energy))
  #print('Nuclear repulsion energy between protons is {} Hartree.'.format(molecule.nuclear_repulsion))
  #for orbital in range(molecule.n_orbitals):
    #print('Spatial orbital {} has energy of {} Hartree.'.format(orbital, molecule.orbital_energies[orbital]))

  # Get the Hamiltonian in an active space.
  # From https://arxiv.org/pdf/1902.10679.pdf
  # "Increasing the representation accuracy of quantum simulations of chemistry withoutextra quantum resources"
  # "quantum computers leverage theability to target a subset of degrees of freedom containing the essential 
  #  quantum behavior, sometimes called the active space"

  # Set Hamiltonian parameters.
  active_space_start = 1
  active_space_stop = 2

  # From OpenFermion Documentation:
  # - occupied_indices (list) - A list of spatial orbital indices indication which orbitals should be considered doubly occupied.
  # - active_indices (list) - A list of spatial orbital indices indicating which orbitals should be condidered active.
  molecular_hamiltonian = molecule.get_molecular_hamiltonian(occupied_indices=range(active_space_start),
		  active_indices=range(active_space_start, active_space_stop))

  # Map operator to fermions and qubits.
  # Once we have the hamiltonian in the active space we can map to the Fermionic space
  fermion_hamiltonian = get_fermion_operator(molecular_hamiltonian)
  print('The Fermionic Hamiltonian is:\n{}'.format(fermion_hamiltonian))

  # And then map to the Qubit space using either JW or BK
  # Constructing the hamiltonians in this way solves the problem of calculating
  # and setting the correct electron-electron interaction integrals
  qubit_hamiltonian_jw = jordan_wigner(fermion_hamiltonian)
  qubit_hamiltonian_jw.compress()
  print('\nThe Jordan-Wigner Hamiltonian in canonical basis follows:\n{}'.format(qubit_hamiltonian_jw))

  qubit_hamiltonian_bk = bravyi_kitaev(fermion_hamiltonian)
  qubit_hamiltonian_bk.compress()
  print('\nThe Bravyi-Kitaev Hamiltonian is:\n{}'.format(qubit_hamiltonian_bk))

  hf_energies += [molecule.hf_energy]
  fci_energies += [molecule.fci_energy]

# Plot.
#import matplotlib.pyplot as plt
#plt.figure(0)
#plt.plot(bond_lengths, fci_energies, 'x-')
#plt.plot(bond_lengths, hf_energies, 'o-')
#plt.ylabel('Energy in Hartree')
#plt.xlabel('Bond length in angstrom')
#plt.show()

