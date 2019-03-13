#!/usr/bin/env python

'''
Teague Tomesh - 3/13/2019

Take in a given Hamiltonian, reference state, ansatz, and optimizer then 
initialize all of the necessary parameters to run the Variational Quantum 
Eigensolver algorithm.

'''

import sys
import time
import getopt
import importlib
import vqeTools


# Main
def main(argv):
  
  start_time = time.time()

  hamPath = ''
  refState = ''
  ansatz = ''
  optimizer = ''
  numQubits = 4

  try:
   opts, args = getopt.getopt(argv,"p:r:a:q:o:",["hamiltonian=","reference=",
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
    '\nAnsatz: ',ansatz,'\nOptimizer: ',optimizer,'\n')

  # Import the modules for the specified reference state, ansatz, and optimizer
  refStateModule = importlib.import_module('.'+refState, package="ReferenceState")
  ansatzModule = importlib.import_module('.'+ansatz, package="Ansatz")
  optimizeModule = importlib.import_module('.'+optimizer, package="Optimizer")

  # Parse the Hamiltonian file as a list of 2-tuples
  # H = [(coef1, ops1), (coef2, ops2), (coef3, ops3), ...]
  hamiltonian, molecule, numElectrons = vqeTools.parseHamiltonian(hamPath)

  # Check that the number of qubits is sufficient to simulate this molecule
  for term in hamiltonian:
    for op in term[1]:
        qIndex = int(op[1])
        if qIndex + 1 > numQubits:
            print("The number of qubits given is too small to simulate this \
                molecule")
            sys.exit()
  
  print('Current hamiltonian is: ', hamiltonian)

  ### PRECOMPUTATION ###
  # Generate a circuit for the reference state
  refCircuit = vqeTools.genRefState(refStateModule, numQubits, numElectrons)

  # Generate a circuit for measuring the energy
  msrCircuit = vqeTools.genMeasure(hamiltonian, numQubits)

  ### VQE ###
  final_energy = optimizeModule.minimizeEnergyObjective(hamiltonian, numQubits,
                                          ansatzModule, refCircuit, msrCircuit)


if __name__ == "__main__":
  main(sys.argv[1:])










