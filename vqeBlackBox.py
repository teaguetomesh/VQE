'''
This file contains functions which allow the optimizer to make calls to 
construct and execute QuantumCircuits.
'''

import sys
import time
import getopt
import importlib
from qiskit import QuantumCircuit
import matplotlib
import random as rand
import math
from measureCircuit import genMeasureCircuit
import hamiltonian_averaging as ha

# Parse the Hamiltonian
# TO-DO: I should actually call the generate-H script so that I can generate
# new Hamiltonians on the fly at new interatomic distances.
# For now, I will just rely on manually calling generate-H to get the
# Hamiltonians I need
def parseHamiltonian(myPath):
    ''' Parse the qubit hamiltonian given at the _myPath_, and return
    the coefficients and operators as a list of 2-tuples where:

       [(float:coef1), (str list:['op','op',...]), ...]

    Returns: 
        H (list): list of tuples holding the coefficients and operators in each
                  term
        molecule (str): Name of the current molecule
        numE (int): Number of electrons in this molecule
    '''
    H = []
    molecule = ''
    numE = 0
    with open(myPath) as hFile:
        for i, line in enumerate(hFile):
            line = line.split()
            if i is 0:
                molecule = line[0]
                numE = int(line[1])
            else:
                coef = float(line[0])
                ops = line[1:]
                H += [(coef, ops)]

    return H, molecule, numE


# Generate reference state
def genRefState(refStateModule, numQubits, numElectrons):
    ''' Generate a QuantumCircuit that will produce the reference state 
    from a quantum register initialized to [0,0,0,...]

    Returns: 
        circuit (QuantumCircuit): A qiskit circuit that generates the ref state
    '''
    circuit = refStateModule.generateReferenceState(numQubits, numElectrons)
    #circuit.draw(filename='refState', output='mpl')
    return circuit



# Generate ansatz
def genAnsatz(ansatzModule, numQubits, params):
    ''' Generate a QuantumCircuit that will generate
    the ansatz state from a quantum register initialized to a reference state.

    Returns:
        circuit (QuantumCircuit): A qiskit circuit that generates the ansatz
    '''
    circuit = ansatzModule.genCircuit(numQubits, params)
    return circuit



# Measure the qubits
def genMeasure(H, numQubits):
    ''' Generate a QuantumCircuit that will measure
    each term of the given hamiltonian. May return multiple circuits if the
    same qubit needs to be measured in perpendicular bases.

    Returns: QuantumCircuit List
    '''
    circuit = genMeasureCircuit(H, numQubits)
    return circuit


# Measure the expected energy
def hamiltonianAveraging(circList, H, Nq):
    '''
    '''
    E = ha.run(circList, H, Nq)
    return E


def constructQuantumCircuit(refCircuit, ansCircuit, msrCircuit):
    '''
    '''
    circList = []
    for tup in msrCircuit:
      mC = tup[0]
      name = tup[1]
      fullCirc = refCircuit + ansCircuit + mC
      fullCirc.draw(filename='vqeCircuit', output='mpl')
      circList += [(fullCirc, name)]
    
    return circList

# Minimize the measured energy
def minimizeEnergyObjective(hamiltonian, numQubits, ansatzModule,
                    refCircuit, msrCircuit, optimizeModule):
    '''
    '''
    initialparams = [rand.uniform(0,2*math.pi) for i in range(20)]

    # Generate a circuit for the ansatz
    ansCircuit = genAnsatz(ansatzModule, numQubits, initialparams)

    circList = constructQuantumCircuit(refCircuit, ansCircuit, msrCircuit)

    # Energy integration
    energy = hamiltonianAveraging(circList, hamiltonian, numQubits)
    print(energy)

'''
# Main
def main(argv):
  
  start_time = time.time()

  hamPath = ''
  refState = ''
  ansatz = ''
  optimizer = ''
  numQubits = 4

  try:
   opts, args = getopt.getopt(argv,"p:r:a:q:",["hamiltonian=","reference=",
   	                          "ansatz=","qubits=","optimizer="])
  except getopt.GetoptError:
    print ('Usage: \n ./vqeBlackBox.py -p <hamiltonian> -r <reference> \
    	   -a <ansatz> -q <qubits> -o <optimizer>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print ('Usage: \n ./vqeBlackBox.py -p <hamiltonian> -r <reference>',
             ' -a <ansatz> -q <qubits> -o <optimizer>')
      sys.exit()
    elif opt in ("-p", "--hamiltonian"):
      hamPath = arg
    elif opt in ("-r", "--reference"):
      refState = arg
    elif opt in ("-a", "--ansatz"):
      ansatz = arg
    elif opt in ("-q", "--qubits"):
      numQubits = int(arg)
    elif opt in ("-o", "--optimizer"):
      optimizer = arg

  print('\nHamiltonian: ',hamPath,'\nReference State: ',refState,
  	'\nAnsatz: ',ansatz,'\n')

  # Import the modules for the specified reference state, ansatz, and optimizer
  refStateModule = importlib.import_module('.'+refState, package="ReferenceState")
  ansatzModule = importlib.import_module('.'+ansatz, package="Ansatz")
  optimizeModule = importlib.import_module('.'+optimizer, package="Optimizer")

  # Parse the Hamiltonian file as a list of 2-tuples
  # H = [(coef1, ops1), (coef2, ops2), (coef3, ops3), ...]
  hamiltonian, molecule, numElectrons = parseHamiltonian(hamPath)

  # Check that the number of qubits is sufficient to simulate this molecule
  for term in hamiltonian:
    for op in term[1]:
        qIndex = int(op[1])
        if qIndex + 1 > numQubits:
            print("The number of qubits given is too small to simulate this \
                molecule")
            sys.exit()
  
  print('Current hamiltonian is: ', hamiltonian)

  ### PRECOMPUTATION
  # Generate a circuit for the reference state
  refCircuit = genRefState(refStateModule, numQubits, numElectrons)

  # Generate a circuit for measuring the energy
  msrCircuit = genMeasure(hamiltonian, numQubits)

  ### VQE
  final_energy = minimizeEnergyObjective(hamiltonian, numQubits, ansatzModule,
                    refCircuit, msrCircuit, optimizeModule)


if __name__ == "__main__":
  main(sys.argv[1:])
'''














