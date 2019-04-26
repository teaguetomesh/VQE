'''
Teague Tomesh - 4/24/2019

Recreate some of the figures found in Ryabinkin et. al.

Constrained Variational Quantum Eigensolver: Quantum Computer Search Engine in the Fock Space
Ilya G. Ryabinkin, Scott N. Genin, and Artur F. Izmaylov
Journal of Chemical Theory and Computation 2019 15 (1), 249-255
DOI: 10.1021/acs.jctc.8b00943

'''


import sys
from pathlib import Path
from openfermion.hamiltonians import MolecularData
import openfermion.hamiltonians as oh
from openfermion.transforms import get_fermion_operator, jordan_wigner, bravyi_kitaev
from openfermionpsi4 import run_psi4
import glob
import numpy as np


def generate_and_save(geometry, basis, multiplicity, description, filename):
    
  # Initialize the molecule
  molecule = MolecularData(geometry, basis, multiplicity, description=description,
                          filename=filename)
  print('Molecule filename: ',filename)
  #molecule.save()
  # Compute the active space integrals
  print('-computing integrals')
  molecule = run_psi4(molecule,run_mp2=False,run_cisd=False,run_ccsd=False,run_fci=True)
  #molecule = run_pyscf(molecule,run_mp2=True,run_cisd=True,run_ccsd=True,run_fci=True)
  print('Successful generation')
    

def load_and_transform(filename, orbitals, transform):
  # Load data
  print('--- Loading molecule ---')
  molecule = MolecularData(filename=filename)
  molecule.load()

  print('filename: {}'.format(molecule.filename))
  print('n_atoms: {}'.format(molecule.n_atoms))
  print('n_electrons: {}'.format(molecule.n_electrons))
  print('n_orbitals: {}'.format(molecule.n_orbitals))
  #print('Canonical Orbitals: {}'.format(molecule.canonical_orbitals))
  print('n_qubits: {}'.format(molecule.n_qubits))
  
  # Get the Hamiltonian in an active space.
  # Set Hamiltonian parameters.
  occupied_orbitals, active_orbitals = orbitals

  molecular_hamiltonian = molecule.get_molecular_hamiltonian(
          occupied_indices=range(occupied_orbitals),
          active_indices=range(active_orbitals))

  # Map operator to fermions and qubits.
  fermion_hamiltonian = get_fermion_operator(molecular_hamiltonian)
  #print('Fermionic Hamiltonian is:\n{}'.format(fermion_hamiltonian))

  if transform is 'JW':
    qubit_h = jordan_wigner(fermion_hamiltonian)
    qubit_h.compress()
    print('\nJordan-Wigner Hamiltonian:\n{}'.format(qubit_h))
  elif transform is 'BK':
    qubit_h = bravyi_kitaev(fermion_hamiltonian)
    qubit_h.compress()
    print('\nBravyi-Kitaev Hamiltonian is:\n{}'.format(qubit_h))
  else:
    print('ERROR: Unrecognized qubit transformation: {}'.format(transform))
    sys.exit(2)

  return qubit_h


def write_to_file(filename, name, Ne, hamiltonian, description, orbitals):
  
  # Write the resulting qubit H to file
  print('\n\n~~write Qubit Hamiltonian to file~~\n')
  print(filename)
  with open(filename, 'w') as H_file:
    H_file.write('{} {} {} {}\n'.format(name, Ne, orbitals[0], orbitals[1]))
    #for h in my_Hs:
    hstring = '{}'.format(hamiltonian)
    print(hstring)
    print('')
    terms = hstring.split('\n')
    for t in terms:
      t2 = t.split('[')
      if len(t2) is 2:
        coef = t2[0]
        paul = t2[1].split(']')[0]
        # Check for identity operator
        if paul is '':
          paul = 'I0'
        
        # Write coefficients and operators to file
        H_file.write('{0:17s} {1}\n'.format(coef, paul))

      else:
        print('ERROR: Something went wrong parsing string')
  print('Successful write\n')


def generate_figure_1():
  singlet_glob = glob.glob('molecule_data/*sto-3g_singlet*.hdf5')
  singlet_glob = sorted(singlet_glob)
  triplet_glob = glob.glob('molecule_data/*sto-3g_triplet*.hdf5')
  triplet_glob = sorted(triplet_glob)

  singlet_fci = [float(MolecularData(filename=fn).fci_energy) for fn in singlet_glob]
  triplet_fci = [MolecularData(filename=fn).fci_energy for fn in triplet_glob]
  xx = [float(MolecularData(filename=fn).description) for fn in singlet_glob]

  #print(singlet_fci)
  #print(triplet_fci)
  #print(xx)


  geometry = [('H', (0., 0., 0.)), ('H', (0., 0., 1.))]
  basis='sto-3g'
  multiplicity=3
  description='1.0'
  filename='test_test_test_H2'
  molecule = MolecularData(geometry, basis, multiplicity, description=description,
                          filename=filename)
  molecule = run_psi4(molecule,run_mp2=False,run_cisd=False,run_ccsd=True,run_fci=True)
  
  #molecule = MolecularData(filename=triplet_glob[0])
  print('filename: {}'.format(molecule.filename))
  print('n_atoms: {}'.format(molecule.n_atoms))
  print('n_electrons: {}'.format(molecule.n_electrons))
  print('n_orbitals: {}'.format(molecule.n_orbitals))
  #print('Canonical Orbitals: {}'.format(molecule.canonical_orbitals))
  print('n_qubits: {}'.format(molecule.n_qubits))
  print(molecule.ccsd_energy)
  print(molecule.description)
  print(molecule.two_body_integrals)
  print(molecule.multiplicity)



def main(argv):
  '''
  '''

  # Set molecule parameters.
  name = 'H2'
  basis = 'sto-3g'
  multiplicity = [1,3]
  n_points = 30
  bond_length_interval = 0.1
  num_electrons = 2
  transform = 'BK'

  for mlt in multiplicity:
    # Generate molecule at different bond lengths.
    for point in range(1, n_points + 1):
      bond_length = bond_length_interval * point
      description = str(round(bond_length,2))
      
      geometry = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length))]
      
      # If this molecule has not been generated
      if mlt is 1:
        mult = 'singlet'
      elif mlt is 3:
        mult = 'triplet'
      
      molecule_file = 'molecule_data/{}_{}_{}_{}.hdf5'.format(name,basis,mult,bond_length)
      config = Path(molecule_file)
      
      if True: #not config.is_file():
        # Generate it now
        print('--- Generate Molecule: {}_{}_{:.2f} ---'.format(name,basis,bond_length))
        generate_and_save(geometry, basis, mlt, description, molecule_file)
      
      # Load the molecule and perform qubit transformation
      occupied = 2
      active = 2
      orbitals = (occupied, active)
      qubit_h = load_and_transform(molecule_file, orbitals, transform)

      # Write the qubit hamiltonian to file
      folder = 'Ryabinkin_data/{}/'.format(mult.capitalize())
      fn = 'qubitH_{}_{}_{}_{}.txt'.format(name, basis, transform, description)
      qubit_file = folder + fn
      write_to_file(qubit_file, name, num_electrons, qubit_h, description, orbitals)

  generate_figure_1()


if __name__ == "__main__":
  #main(sys.argv[1:])
  generate_figure_1()





















