'''
Teague Tomesh - 3/25/2019

Implementation of the UCCSD ansatz for use in the VQE algorithm.

Based on the description given in Whitfield et al.
(https://arxiv.org/abs/1001.3855?context=physics.chem-ph)

Adapted from a Scaffold implementation by Pranav Gokhale]
(https://github.com/epiqc/ScaffCC) 

NOTE:
Qiskit orders their circuits increasing from top -> bottom
  0 ---
  1 ---
  2 ---

Both Whitfield et al. and Barkoutsos et al. order increasing from bottom -> top
  p 3 ---
  q 2 ---
  r 1 ---
  s 0 ---

Not a problem. Qubit index is what matters. Set reverse_bits = True when 
drawing Qiskit circuit.

'''

from qiskit import QuantumCircuit, QuantumRegister
import sys
import math


PI = math.pi


def M_d(i, circ, p, q, r, s, dagger=False):
    '''
    See Double Excitation Operator circuit in Table A1 of Whitfield et al 2010

    Y in Table A1 of Whitfield et al 2010 really means Rx(-pi/2)
    '''

    if dagger:
        angle = PI/2
    else:
        angle = -PI/2

    qr = circ.qregs[0]

    if i == 1:
        circ.h(qr[p])
        circ.h(qr[q])
        circ.h(qr[r])
        circ.h(qr[s])
    elif i == 2:
        circ.rx(angle, qr[p])
        circ.rx(angle, qr[q])
        circ.rx(angle, qr[r])
        circ.rx(angle, qr[s])
    elif i == 3:
        circ.h(qr[p])
        circ.rx(angle, qr[q])
        circ.h(qr[r])
        circ.rx(angle, qr[s])
    elif i == 4:
        circ.rx(angle, qr[p])
        circ.h(qr[q])
        circ.rx(angle, qr[r])
        circ.h(qr[s])
    elif i == 5:
        circ.rx(angle, qr[p])
        circ.rx(angle, qr[q])
        circ.h(qr[r])
        circ.h(qr[s])
    elif i == 6:
        circ.h(qr[p])
        circ.h(qr[q])
        circ.rx(angle, qr[r])
        circ.rx(angle, qr[s])
    elif i == 7:
        circ.rx(angle, qr[p])
        circ.h(qr[q])
        circ.h(qr[r])
        circ.rx(angle, qr[s])
    elif i == 8:
        circ.h(qr[p])
        circ.rx(angle, qr[q])
        circ.rx(angle, qr[r])
        circ.h(qr[s])

    return circ


def CNOTLadder(circ, controlStartIndex, controlStopIndex):
    '''
    Applies a ladder of CNOTs, as in the dashed-CNOT notation at bottom of
    Table A1 of Whitfield et al 2010

    Qubit indices increase from bottom to top
    '''
  
    qr = circ.qregs[0]

    if controlStopIndex > controlStartIndex:
        delta = 1
        index = controlStartIndex + 1
        controlStopIndex += 1
    else:
        delta = -1
        index = controlStartIndex

    while index is not controlStopIndex:
        circ.cx(qr[index], qr[index-1])
        index += delta

    return circ


def DoubleExcitationOperator(circ, theta, p, q, r, s):
    # Prerequisite: p > q > r > s

    qr = circ.qregs[0]

    for i in range(1,9):

        circ = M_d(i, circ, p, q, r, s, dagger=False)

        circ = CNOTLadder(circ, p, q)
        circ.cx(qr[q],qr[r])
        circ = CNOTLadder(circ, r, s)

        circ.rz(theta,qr[s]) # Rz(reg[s], Theta_p_q_r_s[p][q][r][s]);

        circ = CNOTLadder(circ, s, r)
        circ.cx(qr[q],qr[r])
        circ = CNOTLadder(circ, q, p)

        circ.barrier(qr)

        circ = M_d(i, circ, p, q, r, s, dagger=True)

    return circ


def SingleExcitationOperator(circ, theta, p, q):
    # Prerequisite: p > q
    # See Single Excitation Operator circuit in Table A1 of Whitfield et al 2010

    qr = circ.qregs[0]

    circ.barrier(qr)

    circ.h(qr[p])
    circ.h(qr[q])
    circ = CNOTLadder(circ, p, q)
    circ.rz(theta,qr[q]) # Rz(reg[q], Theta_p_q[p][q]);
    circ = CNOTLadder(circ, q, p)

    circ.barrier(qr)

    circ.h(qr[p])
    circ.h(qr[q])

    circ.rx(-PI/2, qr[p])
    circ.rx(-PI/2, qr[q])
    circ = CNOTLadder(circ, p, q)
    circ.rz(theta,qr[q]) # Rz(reg[q], Theta_p_q[p][q]);
    circ = CNOTLadder(circ, q, p)

    circ.barrier(qr)

    circ.rx(-PI/2, qr[p])
    circ.rx(-PI/2, qr[q])

    return circ


def genCircuit(Nq, param):
    '''
    '''
    if Nq is not 4:
        print('ERROR: UCCSD_4_Whitfield is currently implemented for 4 qubits only')
        sys.exit()

    # UCCSD ansatz for 4 qubits takes 7 different angles
    # 1 for the double excitation operator
    # 6 for the single excitation operator
    # Map the parameter indices to a dictionary here, indexed by p,q,r,s strings
    pdict = {'3210':0,'10':1,'20':2,'21':3,'30':4,'31':5,'32':6}

    # Initialize quantum register and circuit
    qreg = QuantumRegister(Nq, name='qreg')
    circ  = QuantumCircuit(qreg, name='UCCSD_4_Whitfield')

    # enumerate all Nq > p > q > r > s >= 0 and apply Double Excitation Operator
    for p in range(Nq):
      for q in range(p):
        for r in range(q): 
          for s in range(r): 
            #print(p,q,r,s)
            # For the 4 qubit case this function is called a single time
            pqrs = str(p)+str(q)+str(r)+str(s)
            circ = DoubleExcitationOperator(circ,param[pdict[pqrs]],p,q,r,s)

    # enumerate all Nq > p > q >= 0 and apply Single Excitation Operator
    for p in range(Nq):
      for q in range(p):
        pq = str(p)+str(q)
        SingleExcitationOperator(circ, param[pdict[pq]], p, q);

    return circ

if __name__ == "__main__":
  circ = genCircuit(4,[0,1,2,3,4,5,6])
  circ.draw(scale=0.8, filename='test_circuit', output='mpl', 
        plot_barriers=False, reverse_bits=True)



