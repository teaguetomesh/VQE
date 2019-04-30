'''
Teague Tomesh - 4/10/2019

Implementation of the UCCSD ansatz for use in the VQE algorithm.

Based on the description given in Barkoutsos et al.
(https://arxiv.org/abs/1001.3855?context=physics.chem-ph)

'''

from qiskit import QuantumCircuit, QuantumRegister
import sys
import math


PI = math.pi


def U_d(i, circ, s, r, q, p, dagger=False):
    '''
    See Double Excitation Operator circuit in Fig 1b. of Barkoutsos et al 2018

    Y in Fig 1b of Barkoutsos et al 2018 really means Rx(-pi/2)
    '''

    if dagger:
        angle = PI/2
    else:
        angle = -PI/2

    qr = circ.qregs[0]

    if i == 1:
        circ.h(qr[s])
        circ.h(qr[r])
        circ.rx(angle, qr[q])
        circ.h(qr[p])
    elif i == 2:
        circ.rx(angle, qr[s])
        circ.h(qr[r])
        circ.rx(angle, qr[q])
        circ.rx(angle, qr[p])
    elif i == 3:
        circ.h(qr[s])
        circ.rx(angle, qr[r])
        circ.rx(angle, qr[q])
        circ.rx(angle, qr[p])
    elif i == 4:
        circ.h(qr[s])
        circ.h(qr[r])
        circ.h(qr[q])
        circ.rx(angle, qr[p])
    elif i == 5:
        circ.rx(angle, qr[s])
        circ.h(qr[r])
        circ.h(qr[q])
        circ.h(qr[p])
    elif i == 6:
        circ.h(qr[s])
        circ.rx(angle, qr[r])
        circ.h(qr[q])
        circ.h(qr[p])
    elif i == 7:
        circ.rx(angle, qr[s])
        circ.rx(angle, qr[r])
        circ.rx(angle, qr[q])
        circ.h(qr[p])
    elif i == 8:
        circ.rx(angle, qr[s])
        circ.rx(angle, qr[r])
        circ.h(qr[q])
        circ.rx(angle, qr[p])

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


def DoubleExcitationOperator(circ, theta, s, r, q, p):
    # Prerequisite: s > r > q > p

    qr = circ.qregs[0]

    for i in range(1,9):

        circ = U_d(i, circ, s, r, q, p, dagger=False)

        circ = CNOTLadder(circ, s, r)
        circ.cx(qr[r],qr[q])
        circ = CNOTLadder(circ, q, p)

        circ.rz(theta,qr[p]) # Rz(reg[s], Theta_p_q_r_s[p][q][r][s]);

        circ = CNOTLadder(circ, p, q)
        circ.cx(qr[r],qr[q])
        circ = CNOTLadder(circ, r, s)

        circ.barrier(qr)

        circ = U_d(i, circ, s, r, q, p, dagger=True)

    return circ


def SingleExcitationOperator(circ, theta, r, p):
    # Prerequisite: r > p
    # See Single Excitation Operator circuit in Fig 1a. of Barkoutsos et al 2018

    qr = circ.qregs[0]

    circ.barrier(qr)

    circ.rx(-PI/2, qr[r])
    circ.h(qr[p])
    circ = CNOTLadder(circ, r, p)
    circ.rz(theta,qr[p]) # Rz(reg[q], Theta_p_q[p][q]);
    circ = CNOTLadder(circ, p, r)

    circ.barrier(qr)

    circ.rx(PI/2, qr[r])
    circ.h(qr[p])

    circ.h(qr[r])
    circ.rx(-PI/2, qr[p])
    circ = CNOTLadder(circ, r, p)
    circ.rz(theta,qr[p]) # Rz(reg[q], Theta_p_q[p][q]);
    circ = CNOTLadder(circ, p, r)

    circ.barrier(qr)

    circ.h(qr[r])
    circ.rx(PI/2, qr[p])

    return circ


def genCircuit(Nq, param):
    '''
    '''
    if Nq is not 4:
        print('ERROR: UCCSD_4_Barkoutsos is currently implemented for 4 qubits only')
        sys.exit()

    # UCCSD ansatz for 4 qubits takes 7 different angles
    # 1 for the double excitation operator
    # 6 for the single excitation operator
    # Map the parameter indices to a dictionary here, indexed by p,q,r,s strings
    pdict = {'3210':0,'10':1,'20':2,'21':3,'30':4,'31':5,'32':6}

    # Initialize quantum register and circuit
    qreg = QuantumRegister(Nq, name='qreg')
    circ  = QuantumCircuit(qreg, name='UCCSD_4_Barkoutsos')

    # enumerate all Nq > s > r > q > p >= 0 and apply Double Excitation Operator
    for s in range(Nq):
      for r in range(s):
        for q in range(r): 
          for p in range(q): 
            #print(s,r,q,p)
            # For the 4 qubit case this function is called a single time
            srqp = str(s)+str(r)+str(q)+str(p)
            circ = DoubleExcitationOperator(circ,param[pdict[srqp]],s,r,q,p)

    # enumerate all Nq > r > p >= 0 and apply Single Excitation Operator
    for r in range(Nq):
      for p in range(r):
        rp = str(r)+str(p)
        SingleExcitationOperator(circ, param[pdict[rp]], r, p);

    return circ

if __name__ == "__main__":
  circ = genCircuit(4,[0,1,2,3,4,5,6])
  circ.draw(scale=0.8, filename='test_circuit', output='mpl', 
        plot_barriers=False, reverse_bits=True)



