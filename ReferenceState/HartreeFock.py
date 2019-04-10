'''
Teague Tomesh - 3/13/2019

Implementation of the Hartree Fock reference state for use in the VQE algorithm.
From [Li, Y. , Hu, J. , Zhang, X. , Song, Z. and Yung, M. (2019), 
Variational Quantum Simulation for Quantum Chemistry. Adv. Theory Simul.. 
doi:10.1002/adts.201800182] Section 5.1 we see that Hartree Fock states can
be easily generated for states using the Jordan-Wigner transformation. 

|HF> = |1>^@(N_e) @ |0>^@(M-N_e)

Where '@' is the tensor-product, M is the number of orbitals, and N_e is the
number of electrons.

'''

from qiskit import QuantumCircuit, QuantumRegister


def generateReferenceState(Nq, Ne):
    '''
    '''

    qr = QuantumRegister(Nq,name='qreg')
    circ = QuantumCircuit(qr)

    for n in range(Ne):
        circ.x(qr[n])

    circ.barrier(qr)

    return circ


if __name__ == "__main__":
  circ = generateReferenceState(4,2)
  circ.draw(scale=0.8, filename='HF_reference_circuit', output='mpl', 
        plot_barriers=False, reverse_bits=True)