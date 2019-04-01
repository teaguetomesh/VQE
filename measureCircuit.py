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
    Check of any of the rows sum to a value >1, if so, then multiple circuits
    will need to be generated.

    Returns:
        circs (List[QuantumCircuits])
    '''

    opDict = {'Z':0, 'X':1, 'Y':2}

    quOpMatrix = np.zeros(shape=(Nq,3), dtype=int)
    for term in H:
        ops = term[1]
        for op in ops:
            if op[0] is 'I': continue
            quOpMatrix[int(op[1]), opDict[op[0]]] = 1

    #print(quOpMatrix)

    numCircs = 0
    for row in quOpMatrix:
        rowSum = np.sum(row)
        if rowSum > numCircs:
            numCircs = rowSum

    circuitList = []
    '''
    TODO: correctly compile each of necessary measurement circuits to 
    satisfy every term in the Hamiltonian
    for i in range(numCircs):
        qr = QuantumRegister(Nq, name='qr')
        cr = ClassicalRegister(Nq, name='cr')
        circ = QuantumCircuit(qr, cr)
        for qubit in range(Nq):
            for op in range(3): 
    '''
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












