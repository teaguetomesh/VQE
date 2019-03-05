from openfermion.hamiltonians import MolecularData
from openfermion.transforms import get_fermion_operator, jordan_wigner, bravyi_kitaev

# Set molecule parameters.
basis = 'sto-3g'
multiplicity = 1
bond_length_interval = 0.1
n_points = 3
name = 'LiH'

# Generate molecule at different bond lengths.
fr_hamiltonians = []
jw_hamiltonians = []
bk_hamiltonians = []
bond_lengths = []
for point in range(3, n_points + 1):
  bond_length = bond_length_interval * point
  bond_length = 1.45
  bond_lengths += [bond_length]
  description = str(round(bond_length,2))
  geometry = [('Li', (0., 0., 0.)), ('H', (0., 0., bond_length))]
  molecule = MolecularData(geometry, basis, multiplicity, description=description)

  # Load data.
  molecule.load()

  print('\nAt bond length of {} angstrom, molecular hydrogen has:'.format(bond_length))
  
  # Get the Hamiltonian in an active space.
  # Set Hamiltonian parameters.
  active_space_start = 1
  active_space_stop = 3

  molecular_hamiltonian = molecule.get_molecular_hamiltonian(
		  occupied_indices=range(active_space_start),
		  active_indices=range(active_space_start, active_space_stop))

  # Map operator to fermions and qubits.
  fermion_hamiltonian = get_fermion_operator(molecular_hamiltonian)
  print('Fermionic Hamiltonian is:\n{}'.format(fermion_hamiltonian))

  qubit_hamiltonian_jw = jordan_wigner(fermion_hamiltonian)
  qubit_hamiltonian_jw.compress()
  print('\nJordan-Wigner Hamiltonian:\n{}'.format(qubit_hamiltonian_jw))

  qubit_hamiltonian_bk = bravyi_kitaev(fermion_hamiltonian)
  qubit_hamiltonian_bk.compress()
  print('\nBravyi-Kitaev Hamiltonian is:\n{}'.format(qubit_hamiltonian_bk))

  fr_hamiltonians += [fermion_hamiltonian]
  jw_hamiltonians += [qubit_hamiltonian_jw]
  bk_hamiltonians += [qubit_hamiltonian_bk]

# Write the resulting qubit H to file
# For now, only considering the Bravyi-Kitaev transformations
print('\n\n~~write Qubit Hamiltonian to file~~\n')
H_file = open('qubitH_{0}_{1}_{2}.txt'.format(name, basis, description), 'w')
for h in bk_hamiltonians:
  hstring = '{}'.format(h)
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
        paul = 'I'
      
      # Write coefficients and operators to file
      H_file.write('{0:17s} {1}\n'.format(coef, paul))

    else:
      print('ERROR: Something went wrong parsing string')
H_file.close()


