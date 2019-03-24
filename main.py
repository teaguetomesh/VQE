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
import glob


# Main
def main(argv):

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

  start_time = time.time()

  # Import the modules for the specified reference state, ansatz, and optimizer
  refStateModule = importlib.import_module('.'+refState, package="ReferenceState")
  ansatzModule = importlib.import_module('.'+ansatz, package="Ansatz")
  optimizeModule = importlib.import_module('.'+optimizer, package="Optimizer")

  if hamPath[-1] == '/':
    # collect all hamiltonians from a folder
    manyH = glob.glob(hamPath+'*')
    sortH = []
    for h in manyH:
      for i in range(len(h)-1,0,-1):
        curChar = h[i]
        found1 = False
        if curChar == '.' and not found1:
          ir = i+2
          found1 = True
        if curChar == '_':
          il = i+1
          break
      r = float(h[il:ir])
      sortH += [(r,h)]
    sortH = sorted(sortH, key=lambda hamil: hamil[0])

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
        if qIndex + 1 > numQubits:
          print("The number of qubits given is too small to simulate this",
                "molecule")
          sys.exit()
    print('------------------------------------')
    print('Begin VQE for {} at r separation: {}'.format(molecule, rParam))
    print('Current hamiltonian is: ', hamiltonian)

    ### PRECOMPUTATION ###
    # Generate a circuit for the reference state
    refCircuit = vqeTools.genRefState(refStateModule, numQubits, numElectrons)
    print('--reference circuit generated--')

    # Generate a circuit for measuring the energy
    msrCircuit = vqeTools.genMeasure(hamiltonian, numQubits)
    print('--measurement circuit generated--')

    ### VQE ###
    print('--initiate {} optimization--'.format(optimizer))
    final_energy = optimizeModule.minimizeEnergyObjective(hamiltonian, numQubits,
                                           ansatzModule, refCircuit, msrCircuit,
                                           prevParam)

    print('--optimization finished--')
    print('\tAchieved {0:.7f} Ha at {1} angstroms'.format(final_energy[1], rParam))
    results += [(rParam, final_energy)]
    prevParam = final_energy[0]

  end_time = time.time()

  # After VQE is finished - write results to file
  for i, c in enumerate(hamPath):
    if c == '/':
      fname = hamPath[i+1:-1]+'_{}_to_{}'.format(sortH[0][0], sortH[-1][0])
      break
  with open('Results/{}.txt'.format(fname),'w') as outf:
    for p in results:
      writeStr = '{0} {1}   '.format(str(p[0]).ljust(5),str(p[1][1]).rjust(20))
      for s in p[1][0]:
        writeStr += ' {:.5f}'.format(s)
      writeStr += '\n'
      outf.write(writeStr)

  elapsed_time = end_time - start_time
  timestr = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
  print('Elapsed time: ', timestr)

if __name__ == "__main__":
  main(sys.argv[1:])










