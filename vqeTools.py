'''
Teague Tomesh - 3/13/2019

These functions give the optimizer the ability to construct and execute
QuantumCircuits in a programmatic fashion.

'''

import os
import sys
import re
import time
import getopt
import importlib
from qiskit import QuantumCircuit
import matplotlib
import random as rand
import math
from measureCircuit import genMeasureCircuit
import hamiltonian_averaging as ha
import glob

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


def constructQuantumCircuit(refCircuit, ansCircuit, msrCircuits):
    '''
    Given 3 QuantumCircuits, concatenate them together into a single circuit.
    '''
    print_circuit = False
    circList = []
    for n, tup in enumerate(msrCircuits):
      mC, name = tup
      fullCirc = refCircuit + ansCircuit + mC
      if print_circuit:
          qasmCirc = fullCirc.qasm()
          qasm_dir = 'Results/QASM_circuits/'
      
          if not os.path.isdir(qasm_dir):
              os.mkdir(qasm_dir)
              print('Directory {} created.'.format(qasm_dir))
      
          all_qasms = sorted(glob.glob(qasm_dir+'/*'))
          last_num = len(all_qasms)+1
          qasm_fn = qasm_dir + '{}_{}.qasm'.format(fullCirc.name,last_num)
          with open(qasm_fn, 'w') as qfn:
              print('Writing QASM to {}'.format(qasm_fn))
              qfn.write(qasmCirc)
          #fullCirc.draw(scale=0.8, filename='measure_{0}_{1}_{2}'.format(ansCircuit.name,name,n), 
          #  output='mpl', plot_barriers=False, reverse_bits=True)
          #sys.exit()
      circList += [(fullCirc, name)]
    
    return circList


















