#!/usr/bin/env python

'''
Teague Tomesh - 3/13/2019

Take in a given Hamiltonian, reference state, ansatz, and optimizer then 
initialize all of the necessary parameters to run the Variational Quantum 
Eigensolver algorithm.

'''

import sys
import argparse
import time
import getopt
import importlib
import vqeTools
import glob


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--hamiltonian", type=str,
            default=None, help="Path to Hamiltonian")
    parser.add_argument("-r","--refstate", type=str,
            default=None,
            help="Name of the reference state to use (Look in ReferenceState/)")
    parser.add_argument("-a","--ansatz", type=str,
            default=None, help="Name of the ansatz to use (Look in Ansatz/)")
    parser.add_argument("-q","--qubits", type=int,
            default=None, help="Number of qubits")
    parser.add_argument("-o","--optimizer", type=str,
            default=None,
            help="Name of the optimizer to use (Look in Optimizer/)")
    parser.add_argument("-w","--output", type=str,
            default="DEFAULT_OUTPUT.txt", help="Path to the output file")
    args = parser.parse_args()
    return args


def get_distance(h):
  '''
  Parse the Hamiltonian path to get the current separation distance (Angstroms)
  '''
  for i in range(len(h)-1,0,-1):
    curChar = h[i]
    found1 = False
    if curChar == '.' and not found1:
      ir = i+2
      found1 = True
    if curChar == '_':
      il = i+1
      break
  return float(h[il:ir])


# Main
def main():

  # Parse input args
  args = parse_args()

  print('\nHamiltonian: ',args.hamiltonian,'\nReference State: ',args.refstate,
    '\nAnsatz: ',args.ansatz,'\nOptimizer: ',args.optimizer,'\nOutput: ',args.output)

  # Given the desired ansatz and number of qubits, how many parameters are needed?
  num_parameters = -1
  if 'naive' in args.ansatz:
    num_parameters = 20
  elif 'UCCSD' in args.ansatz:
    if args.qubits == 2:
      num_parameters = 1
    elif args.qubits == 4:
      num_parameters = 7
    elif args.qubits == 6:
      num_parameters = 30
    elif args.qubits == 8:
      num_parameters = 98
  if num_parameters == -1:
    print('ERROR: Current ansatz & qubit number not supported')
    sys.exit(2)

  start_time = time.time()

  # Import the modules for the specified reference state, ansatz, and optimizer
  refStateModule = importlib.import_module('.'+args.refstate, package="ReferenceState")
  ansatzModule = importlib.import_module('.'+args.ansatz, package="Ansatz")
  optimizeModule = importlib.import_module('.'+args.optimizer, package="Optimizer")

  # Get Hamiltonian 
  sortH = []
  if args.hamiltonian[-1] == '/':
    # collect all hamiltonians from a folder
    manyH = glob.glob(args.hamiltonian+'*')
    for h in manyH:
      r = get_distance(h)
      sortH += [(r,h)]
    sortH = sorted(sortH, key=lambda hamil: hamil[0])
  else:
    sortH += [(get_distance(args.hamiltonian),args.hamiltonian)]

  #####################################
  ######## Begin Main VQE Loop ########
  results = []
  prevParam = None
  for pair in sortH:
    rParam = pair[0]
    hPath = pair[1]
    # Parse the Hamiltonian file as a list of 2-tuples
    # H = [(coef1, ops1), (coef2, ops2), (coef3, ops3), ...]
    hamiltonian, molecule, numElectrons = vqeTools.parseHamiltonian(hPath)

    # Check that the number of qubits is sufficient to simulate this molecule
    for term in hamiltonian:
      for op in term[1]:
        qIndex = int(op[1])
        if qIndex + 1 > args.qubits:
          print("The number of qubits given is too small to simulate this",
                "molecule")
          sys.exit()
    print('------------------------------------')
    print('Begin VQE for {} at r separation: {}'.format(molecule, rParam))
    print('Current hamiltonian is: ', hamiltonian)

    ### PRECOMPUTATION ###
    # Generate a circuit for the reference state
    refCircuit = vqeTools.genRefState(refStateModule, args.qubits, numElectrons)
    print('--reference circuit generated--')

    # Generate a circuit for measuring the energy
    msrCircuits = vqeTools.genMeasure(hamiltonian, args.qubits)
    print('--measurement circuit generated--')

    ### VQE ###
    print('--initiate {} optimization--'.format(args.optimizer))
    final_energy = optimizeModule.minimizeEnergyObjective(hamiltonian, args.qubits,
                                      ansatzModule, refCircuit, msrCircuits,
                                      prevParam, num_parameters)

    print('--optimization finished--')
    print('\tAchieved {0:.7f} Ha at {1} angstroms'.format(final_energy[1], rParam))
    results += [(rParam, final_energy)]
    prevParam = final_energy[0]

  end_time = time.time()

  # After VQE is finished - write results to file
  fname = args.output
  with open(fname,'w') as outf:
    for p in results:
      writeStr = '{0} {1}   '.format(str(p[0]).ljust(5),str(p[1][1]).rjust(20))
      for s in p[1][0]:
        writeStr += ' {:.5f}'.format(s)
      writeStr += '\n'
      outf.write(writeStr)

  et = end_time - start_time
  print('Elapsed time (H:M:S): {:02d}:{:02d}:{:02d}'.format(int(et // 3600),
      int(et % 3600 // 60), int(et % 60)))

if __name__ == "__main__":
  main()










