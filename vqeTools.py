'''
Teague Tomesh - 3/13/2019

These functions give the optimizer the ability to construct and execute
QuantumCircuits in a programmatic fashion.

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
    for n, tup in enumerate(msrCircuit):
      mC = tup[0]
      name = tup[1]
      fullCirc = refCircuit + ansCircuit + mC
      #fullCirc.draw(scale=0.8, filename='vqeCirc_{0}_{1}_{2}'.format(ansCircuit.name,name,n), 
      #  output='mpl', plot_barriers=False, reverse_bits=True)
      #sys.exit()
      circList += [(fullCirc, name)]
    
    return circList


















