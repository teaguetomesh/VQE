'''

Implementation of a naive ansatz for use in the VQE algorithm.
Adapted from [EPiQC VQE tutorial by Pranav Gokhale]
(https://www.youtube.com/watch?v=E947xs9-Mso)

In the tutorial, this ansatz was designed for 2 qubits, here I am extending
it to 4 qubits. 

'''

from qiskit import QuantumCircuit, QuantumRegister
import sys


def genCircuit(M, p):
    '''
    '''
    if M is not 4:
        print('ERROR: The naive ansatz is currently implemented for 4 qubits \
            only')
        sys.exit()

    # Initialize quantum register and circuit
    qr = QuantumRegister(M, name='qreg')
    c  = QuantumCircuit(qr, name='naive_ansatz')

    # Perform initial rotations
    for i in range(4):
        c.rx(p[i], qr[i])
        c.rz(p[i+4], qr[i])

    # Perform entangling operations
    for i in range(3):
        # ladder down
        c.h(qr[i])
        c.cx(qr[i],qr[i+1])
    c.h(qr[3])
    for i in range(3,0,-1):
        # ladder up
        c.cx(qr[i],qr[i-1])

    # Perform final rotations
    for i in range(4):
        c.rx(p[i+8] ,qr[i])
        c.rz(p[i+12],qr[i])
        c.rz(p[i+16],qr[i])

    return c



