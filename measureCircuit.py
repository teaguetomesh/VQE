'''
Teague Tomesh - 3/13/2019

Given a particular qubit Hamiltonian, measuring the expected energy of any
given quantum state will depend only on the individual terms of that 
Hamiltonian. 

measureCircuit.py generates a circuit which will measure a quantum state in the
correct bases to allow the energy to be calculated. This may require generating
multiple circuits if the same qubit needs to be measured in two perpendicular
bases (i.e. Z and X).

'''


import sys
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
import numpy as np


def genMeasureCircuit(H, Nq):
    ''' Take in a given Hamiltonian, H, and produce the minimum number of 
    necessary circuits to measure each term of H.
    Initialize a numpy array with dimensions [Nq,3] with all values equal to 0,
    where each column of the matrix represents a measurement in one of the 
    Pauli bases (Z,X,Y).
    For each term in H, set the bit corresponding to the operator acting on a
    specific qubit to high. 
    Check if any of the rows sum to a value >1, if so, then multiple circuits
    will need to be generated.

    Returns:
        List[QuantumCircuits]
    '''

    opDict = {'Z':0, 'X':1, 'Y':2}
    opList = ('Z','X','Y')

    quOpMatrix = np.zeros(shape=(Nq,3), dtype=int)
    for term in H:
        ops = term[1]
        for op in ops:
            if op[0] is 'I': continue
            quOpMatrix[int(op[1]), opDict[op[0]]] = 1

    numCircs = 0
    for row in quOpMatrix:
        rowSum = np.sum(row)
        if rowSum > numCircs:
            numCircs = rowSum  

    circuitList = []
    
    # TODO: correctly compile each of necessary measurement circuits to 
    # satisfy every term in the Hamiltonian

    # Generate the circuitMatrix
    #   Number of rows = number of qubits
    #   Number of cols = number of circuits
    #   Each entry will have a value v <- {'Z','X','Y'} corresponding to the
    #   basis this qubit must be measured in.
    circuitMatrix = []
    print(quOpMatrix)
    for i in range(Nq):
        cur_qubit_entries = []
        m = 0
        while len(cur_qubit_entries) < numCircs:
            m = m % 3
            if quOpMatrix[i,m] == 1:
                cur_qubit_entries += [opList[m]]
            m += 1
        circuitMatrix += [cur_qubit_entries]
    circuitMatrix = np.array(circuitMatrix)
    print(circuitMatrix)
    sys.exit()

    qr = QuantumRegister(Nq, name='qreg')
    cr = ClassicalRegister(Nq, name='creg')
    circ = QuantumCircuit(qr, cr)
    name = ''
    circ.barrier(qr)
    for n in range(Nq):
      circ.measure(qr[n],cr[n])
      name += 'Z'
    
    circuitList += [(circ, name)]

    return circuitList


if __name__ == "__main__":
  H = [(5.076946850678632, ['I0']), (-0.006811585442824442, ['X0', 'X1', 'Y2', 'Y3']), (0.006811585442824442, ['X0', 'Y1', 'Y2', 'X3']), (0.006811585442824442, ['Y0', 'X1', 'X2', 'Y3']), (-0.006811585442824442, ['Y0', 'Y1', 'X2', 'X3']), (-0.4131939082582367, ['Z0']), (0.24086324970819822, ['Z0', 'Z1']), (0.09808757340260295, ['Z0', 'Z2']), (0.10489915884542617, ['Z0', 'Z3']), (-0.41319390825823665, ['Z1']), (0.10489915884542617, ['Z1', 'Z2']), (0.09808757340260295, ['Z1', 'Z3']), (-0.6600956717966215, ['Z2']), (0.09255914539024543, ['Z2', 'Z3']), (-0.6600956717966215, ['Z3'])]
  Nq = 4 
  genMeasureCircuit(H, Nq)











