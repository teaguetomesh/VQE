'''
Teague Tomesh - 3/13/2019

Given a circuit and Hamiltonian, run multiple simulations of that circuit to
determine the expected energy of the given Hamiltonian.

'''

import sys
import time
import multiprocessing as mp
from qiskit import BasicAer, execute
import visualization as vis
import numpy as np


def match(cs, ts):
    '''
    '''
    truthVal = True
    for c_cs, c_ts in zip(cs, ts):
        if c_ts == '*':
            continue
        elif c_cs != c_ts:
            truthVal = False

    return truthVal

def matchTermToCircuit(nameList, H, Nq):
    '''
    Which terms can each circuit measure?

    Input:
        nameList (List [str]): contains the names of the different available
                     circuits in the format: 'ZZZZ', 'ZXZX', ....
                     where the characters specify which basis that qubit
                     is measured in.

        H (List [(float, [str])]): representation of the Hamiltonian

        Nq (int): number of qubits used in this simulation

    Returns:
        matching (dict): {'circName':[list of term indices], ...}

        termStrings (List [(str, int)])
    '''

    # TEST: this H and nameList should result in the matching dictionary:
    #       {'ZZZZ':[2,6], 'XXXX':[1,3], 'XZXZ':[4,5]}
    #H = [(0.9, ['I0']), (0.8, ['X1']), (0.1, ['Z0', 'Z1']), (0.1, ['X0', 'X1']),
    #     (0.1, ['X0', 'Z1']), (0.1, ['X2', 'Z3']), (0.8, ['Z0'])]

    #nameList = ['ZZZZ','XXXX','XZXZ']

    #print(H)

    termStrings = ['' for i in range(len(H))]
    #print('empty TS: ',termStrings)
    # Construct a string which reflects the measurement needs of each term
    for n, t in enumerate(H):
        termString = ''
        ops = t[1]
        if ops[0] == 'I0':
            continue
        for op in ops:
            pauli = op[0]
            index = int(op[1])
            while len(termString) < index:
                # Assuming the pauli operators are given in order of 
                # increasing qubit index, add '*' to the string until
                # we get to the index of the current operator
                termString += '*'
            else:
                # Now the term string should get the operator at this index
                termString += pauli
        else:
            # add '*' to the string until its length equals the number of qubits
            while len(termString) < Nq:
                termString += '*'

        termStrings[n] = termString

    # To see the relation between termStrings and H, uncomment this loop
    #for t, h in zip(termStrings, H):
    #    print(t,h)

    #print('termStrings: ',termStrings)

    matching = {}
    termIndices = np.arange(1,len(H))
    for circName in nameList:
        matchedTerms, delList = [], []
        #nTerms = len(termIndices)
        for i, t in enumerate(termIndices):
            ts = termStrings[t]
            if match(circName, ts):
                #print('it matched!')
                matchedTerms += [t]
                delList += [i]
        matching[circName] = matchedTerms
        termIndices = [e for i, e in enumerate(termIndices) if i not in delList]
        #print('del termIndices: ', termIndices)

    #print('matching: ', matching)

    return matching, termStrings


def computeParity(ts, ss):
    '''
    Input:
        ts (str): termstring, i.e.  XXZY or **ZZY**YZZ
        ss (str): statestring, i.e. 0000 or 0101000110
    '''
    val = 1
    for c_ts, c_ss in zip(ts, ss):
        if c_ts == '*':
            # this particular term does not care about this particular qubit
            continue
        if c_ss == '1':
            # measurements of 1 contribute factor -1 --> measuring the 1 state
            # means this qubit is in the |1> eigenspace with eigenvalue = -1
            val = -1*val
        #if c_ss == '0' then...    
            # measurements of 0 contribute factor +1 --> the |0> eigenspace has
            # eigenvalue = +1
            # (do not need to change anything in this case)
    return val


def energyIntegration(H, terms, counts, termStrings):
    '''
    Use the given simulation results to compute the average measured energy

    Input:
        H (List [(float, [str])]): represenation of the Hamiltonian

        terms (List [int]): which terms to consider for the given counts

        counts (Dict): simulation results in key-value format. The keys are
              strings representing the state measured (i.e. '0000', '0011'), 
              and the values are the number of times that state was measured.

        termStrings (List [(str, int)])

    Returns:
        energySum (float): the average energy contributed by these terms
                           of the form:

                E = coef1*(()+()+...)/#shots + coef2*(()+()+...)/#shots + ...
    '''

    energySum = 0
    numShots = np.sum(list(counts.values()))

    for term in terms:
        coef = H[term][0]
        ts = termStrings[term]

        # for this term, average its energy contribution by counting the 
        # number of times it contributed +1 or -1
        weight = 0
        for state in counts:
            count = counts[state]
            parity = computeParity(ts, state)
            weight += (parity*count)

        # multiply by the term's coefficient
        energySum += (coef*weight)

    return (energySum / numShots)


def fire(queue,tup,simulator,H,matching,termStrings):
    circ = tup[0]
    name = tup[1]
    #print('IN FIRE, CIRCNAME = {}'.format(name))
    result = execute(circ, simulator).result()
    counts = result.get_counts(circ)
    #print(counts)
      
    #vis.state_histogram(counts,'Naive Ansatz Counts: {} shots'.format(shots))

    queue.put(energyIntegration(H, matching[name], counts, termStrings))


def run(circList, H, Nq):
    '''
    Perform Hamiltonian averaging or energy integration to determine the
    expected value of the energy for the given H. This is done in a term by
    term fashion because the given H consists of a sum of Pauli terms.

    Input:
        circList (List (QuantumCircuit)): 

        H (List): [(coef1, ops1), (coef2, ops2), (coef3, ops3), ...]

        Nq (int): number of qubits

    Returns:
        totalE (float): estimation for the expected energy
    '''
    #print('IN RUN, START TIME')
    start_time = time.time()
    nameList = [t[1] for t in circList]
    matching, termStrings = matchTermToCircuit(nameList, H, Nq)

    #TODO: add logic to change the number of shots based on the desired 
    #      precision. Li et al: for desired precision, e, need O(w^2/e^2) shots.
    shots   = 1000
    simulator = BasicAer.get_backend('qasm_simulator')
    # The first term of the Hamiltonian is usually the identity operator with
    # coefficient in front. That coefficient is the starting point for
    # the energy calculation.
    multiproc_output = mp.Queue()
    processes = [mp.Process(target=fire,args=([multiproc_output,t,simulator,H,matching,termStrings])) for t in circList]
    #print('IN HA, NUM PROCESSES = {}'.format(len(processes)))
    #print('NUM CIRCUITS = {}'.format(len(circList)))
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    
    results = [multiproc_output.get() for p in processes]
    
    #print(results)
    #print(np.sum(results))
    #end_time = time.time()
    #elapsed_time = end_time - start_time
    #timestr = time.strftime("%H:%M:%S",time.gmtime(elapsed_time))
    #print('FINISHED, elapsed time: ',timestr)

    return np.sum(results)


if __name__ == "__main__":
  #H = [(0.9, ['I0']), (0.8, ['X1']), (0.1, ['Z0', 'Z1']), (0.1, ['X0', 'X1']),(0.1, ['X0', 'Z1']), (0.1, ['X2', 'Z3']), (0.8, ['Z0'])]
  H = [(1, ['I0']), (1, ['X0', 'X1', 'Y2', 'Y3']), (1, ['X0', 'X1', 'Y2', 'Z3', 'Z4', 'Z5', 'Z6', 'Y7']), (1, ['X0', 'X1', 'X3', 'Z4', 'Z5', 'X6']), (1, ['X0', 'X1', 'Y4', 'Y5']), (1, ['X0', 'X1', 'Y6', 'Y7']), (1, ['X0', 'Y1', 'Y2', 'X3']), (1, ['X0', 'Y1', 'Y2', 'Z3', 'Z4', 'Z5', 'Z6', 'X7']), (1, ['X0', 'Y1', 'Y3', 'Z4', 'Z5', 'X6']), (1, ['X0', 'Y1', 'Y4', 'X5']), (1, ['X0', 'Y1', 'Y6', 'X7']), (1, ['X0', 'Z1', 'X2', 'X3', 'Z4', 'X5']), (1, ['X0', 'Z1', 'X2', 'Y3', 'Z4', 'Y5']), (1, ['X0', 'Z1', 'X2', 'X4', 'Z5', 'X6']), (1, ['X0', 'Z1', 'X2', 'Y4', 'Z5', 'Y6']), (1, ['X0', 'Z1', 'X2', 'X5', 'Z6', 'X7']), (1, ['X0', 'Z1', 'X2', 'Y5', 'Z6', 'Y7']), (1, ['X0', 'Z1', 'Y2', 'Y4', 'Z5', 'X6']), (1, ['X0', 'Z1', 'Z2', 'X3', 'Y4', 'Z5', 'Z6', 'Y7']), (1, ['X0', 'Z1', 'Z2', 'X3', 'X5', 'X6']), (1, ['X0', 'Z1', 'Z2', 'Y3', 'Y4', 'Z5', 'Z6', 'X7']), (1, ['X0', 'Z1', 'Z2', 'Y3', 'Y5', 'X6']), (1, ['X0', 'Z1', 'Z2', 'Z3', 'X4']), (1, ['X0', 'Z1', 'Z2', 'Z3', 'X4', 'Z5']), (1, ['X0', 'Z1', 'Z2', 'Z3', 'X4', 'Z6']), (1, ['X0', 'Z1', 'Z2', 'Z3', 'X4', 'Z7']), (1, ['X0', 'Z1', 'Z2', 'Z3', 'Z4', 'X5', 'Y6', 'Y7']), (1, ['X0', 'Z1', 'Z2', 'Z3', 'Z4', 'Y5', 'Y6', 'X7']), (1, ['X0', 'Z1', 'Z2', 'X4']), (1, ['X0', 'Z1', 'Z3', 'X4']), (1, ['X0', 'Z2', 'Z3', 'X4']), (1, ['Y0', 'X1', 'X2', 'Y3']), (1, ['Y0', 'X1', 'X2', 'Z3', 'Z4', 'Z5', 'Z6', 'Y7']), (1, ['Y0', 'X1', 'X3', 'Z4', 'Z5', 'Y6']), (1, ['Y0', 'X1', 'X4', 'Y5']), (1, ['Y0', 'X1', 'X6', 'Y7']), (1, ['Y0', 'Y1', 'X2', 'X3']), (1, ['Y0', 'Y1', 'X2', 'Z3', 'Z4', 'Z5', 'Z6', 'X7']), (1, ['Y0', 'Y1', 'Y3', 'Z4', 'Z5', 'Y6']), (1, ['Y0', 'Y1', 'X4', 'X5']), (1, ['Y0', 'Y1', 'X6', 'X7']), (1, ['Y0', 'Z1', 'X2', 'X4', 'Z5', 'Y6']), (1, ['Y0', 'Z1', 'Y2', 'X3', 'Z4', 'X5']), (1, ['Y0', 'Z1', 'Y2', 'Y3', 'Z4', 'Y5']), (1, ['Y0', 'Z1', 'Y2', 'X4', 'Z5', 'X6']), (1, ['Y0', 'Z1', 'Y2', 'Y4', 'Z5', 'Y6']), (1, ['Y0', 'Z1', 'Y2', 'X5', 'Z6', 'X7']), (1, ['Y0', 'Z1', 'Y2', 'Y5', 'Z6', 'Y7']), (1, ['Y0', 'Z1', 'Z2', 'X3', 'X4', 'Z5', 'Z6', 'Y7']), (1, ['Y0', 'Z1', 'Z2', 'X3', 'X5', 'Y6']), (1, ['Y0', 'Z1', 'Z2', 'Y3', 'X4', 'Z5', 'Z6', 'X7']), (1, ['Y0', 'Z1', 'Z2', 'Y3', 'Y5', 'Y6']), (1, ['Y0', 'Z1', 'Z2', 'Z3', 'Y4']), (1, ['Y0', 'Z1', 'Z2', 'Z3', 'Y4', 'Z5']), (1, ['Y0', 'Z1', 'Z2', 'Z3', 'Y4', 'Z6']), (1, ['Y0', 'Z1', 'Z2', 'Z3', 'Y4', 'Z7']), (1, ['Y0', 'Z1', 'Z2', 'Z3', 'Z4', 'X5', 'X6', 'Y7']), (1, ['Y0', 'Z1', 'Z2', 'Z3', 'Z4', 'Y5', 'X6', 'X7']), (1, ['Y0', 'Z1', 'Z2', 'Y4']), (1, ['Y0', 'Z1', 'Z3', 'Y4']), (1, ['Y0', 'Z2', 'Z3', 'Y4']), (1, ['Z0']), (1, ['Z0', 'X1', 'Z2', 'Z3', 'Z4', 'X5']), (1, ['Z0', 'Y1', 'Z2', 'Z3', 'Z4', 'Y5']), (1, ['Z0', 'Z1']), (1, ['Z0', 'X2', 'Z3', 'Z4', 'Z5', 'X6']), (1, ['Z0', 'Y2', 'Z3', 'Z4', 'Z5', 'Y6']), (1, ['Z0', 'Z2']), (1, ['Z0', 'X3', 'Z4', 'Z5', 'Z6', 'X7']), (1, ['Z0', 'Y3', 'Z4', 'Z5', 'Z6', 'Y7']), (1, ['Z0', 'Z3']), (1, ['Z0', 'Z4']), (1, ['Z0', 'Z5']), (1, ['Z0', 'Z6']), (1, ['Z0', 'Z7']), (1, ['X1', 'X2', 'Y3', 'Y4']), (1, ['X1', 'X2', 'X4', 'Z5', 'Z6', 'X7']), (1, ['X1', 'X2', 'Y5', 'Y6']), (1, ['X1', 'Y2', 'Y3', 'X4']), (1, ['X1', 'Y2', 'Y4', 'Z5', 'Z6', 'X7']), (1, ['X1', 'Y2', 'Y5', 'X6']), (1, ['X1', 'Z2', 'X3', 'X4', 'Z5', 'X6']), (1, ['X1', 'Z2', 'X3', 'Y4', 'Z5', 'Y6']), (1, ['X1', 'Z2', 'X3', 'X5', 'Z6', 'X7']), (1, ['X1', 'Z2', 'X3', 'Y5', 'Z6', 'Y7']), (1, ['X1', 'Z2', 'Y3', 'Y5', 'Z6', 'X7']), (1, ['X1', 'Z2', 'Z3', 'X4', 'X6', 'X7']), (1, ['X1', 'Z2', 'Z3', 'Y4', 'Y6', 'X7']), (1, ['X1', 'Z2', 'Z3', 'Z4', 'X5']), (1, ['X1', 'Z2', 'Z3', 'Z4', 'X5', 'Z6']), (1, ['X1', 'Z2', 'Z3', 'Z4', 'X5', 'Z7']), (1, ['X1', 'Z2', 'Z3', 'X5']), (1, ['X1', 'Z2', 'Z4', 'X5']), (1, ['X1', 'Z3', 'Z4', 'X5']), (1, ['Y1', 'X2', 'X3', 'Y4']), (1, ['Y1', 'X2', 'X4', 'Z5', 'Z6', 'Y7']), (1, ['Y1', 'X2', 'X5', 'Y6']), (1, ['Y1', 'Y2', 'X3', 'X4']), (1, ['Y1', 'Y2', 'Y4', 'Z5', 'Z6', 'Y7']), (1, ['Y1', 'Y2', 'X5', 'X6']), (1, ['Y1', 'Z2', 'X3', 'X5', 'Z6', 'Y7']), (1, ['Y1', 'Z2', 'Y3', 'X4', 'Z5', 'X6']), (1, ['Y1', 'Z2', 'Y3', 'Y4', 'Z5', 'Y6']), (1, ['Y1', 'Z2', 'Y3', 'X5', 'Z6', 'X7']), (1, ['Y1', 'Z2', 'Y3', 'Y5', 'Z6', 'Y7']), (1, ['Y1', 'Z2', 'Z3', 'X4', 'X6', 'Y7']), (1, ['Y1', 'Z2', 'Z3', 'Y4', 'Y6', 'Y7']), (1, ['Y1', 'Z2', 'Z3', 'Z4', 'Y5']), (1, ['Y1', 'Z2', 'Z3', 'Z4', 'Y5', 'Z6']), (1, ['Y1', 'Z2', 'Z3', 'Z4', 'Y5', 'Z7']), (1, ['Y1', 'Z2', 'Z3', 'Y5']), (1, ['Y1', 'Z2', 'Z4', 'Y5']), (1, ['Y1', 'Z3', 'Z4', 'Y5']), (1, ['Z1']), (1, ['Z1', 'X2', 'Z3', 'Z4', 'Z5', 'X6']), (1, ['Z1', 'Y2', 'Z3', 'Z4', 'Z5', 'Y6']), (1, ['Z1', 'Z2']), (1, ['Z1', 'X3', 'Z4', 'Z5', 'Z6', 'X7']), (1, ['Z1', 'Y3', 'Z4', 'Z5', 'Z6', 'Y7']), (1, ['Z1', 'Z3']), (1, ['Z1', 'Z4']), (1, ['Z1', 'Z5']), (1, ['Z1', 'Z6']), (1, ['Z1', 'Z7']), (1, ['X2', 'X3', 'Y4', 'Y5']), (1, ['X2', 'X3', 'Y6', 'Y7']), (1, ['X2', 'Y3', 'Y4', 'X5']), (1, ['X2', 'Y3', 'Y6', 'X7']), (1, ['X2', 'Z3', 'X4', 'X5', 'Z6', 'X7']), (1, ['X2', 'Z3', 'X4', 'Y5', 'Z6', 'Y7']), (1, ['X2', 'Z3', 'Z4', 'Z5', 'X6']), (1, ['X2', 'Z3', 'Z4', 'Z5', 'X6', 'Z7']), (1, ['X2', 'Z3', 'Z4', 'X6']), (1, ['X2', 'Z3', 'Z5', 'X6']), (1, ['X2', 'Z4', 'Z5', 'X6']), (1, ['Y2', 'X3', 'X4', 'Y5']), (1, ['Y2', 'X3', 'X6', 'Y7']), (1, ['Y2', 'Y3', 'X4', 'X5']), (1, ['Y2', 'Y3', 'X6', 'X7']), (1, ['Y2', 'Z3', 'Y4', 'X5', 'Z6', 'X7']), (1, ['Y2', 'Z3', 'Y4', 'Y5', 'Z6', 'Y7']), (1, ['Y2', 'Z3', 'Z4', 'Z5', 'Y6']), (1, ['Y2', 'Z3', 'Z4', 'Z5', 'Y6', 'Z7']), (1, ['Y2', 'Z3', 'Z4', 'Y6']), (1, ['Y2', 'Z3', 'Z5', 'Y6']), (1, ['Y2', 'Z4', 'Z5', 'Y6']), (1, ['Z2']), (1, ['Z2', 'X3', 'Z4', 'Z5', 'Z6', 'X7']), (1, ['Z2', 'Y3', 'Z4', 'Z5', 'Z6', 'Y7']), (1, ['Z2', 'Z3']), (1, ['Z2', 'Z4']), (1, ['Z2', 'Z5']), (1, ['Z2', 'Z6']), (1, ['Z2', 'Z7']), (1, ['X3', 'X4', 'Y5', 'Y6']), (1, ['X3', 'Y4', 'Y5', 'X6']), (1, ['X3', 'Z4', 'Z5', 'Z6', 'X7']), (1, ['X3', 'Z4', 'Z5', 'X7']), (1, ['X3', 'Z4', 'Z6', 'X7']), (1, ['X3', 'Z5', 'Z6', 'X7']), (1, ['Y3', 'X4', 'X5', 'Y6']), (1, ['Y3', 'Y4', 'X5', 'X6']), (1, ['Y3', 'Z4', 'Z5', 'Z6', 'Y7']), (1, ['Y3', 'Z4', 'Z5', 'Y7']), (1, ['Y3', 'Z4', 'Z6', 'Y7']), (1, ['Y3', 'Z5', 'Z6', 'Y7']), (1, ['Z3']), (1, ['Z3', 'Z4']), (1, ['Z3', 'Z5']), (1, ['Z3', 'Z6']), (1, ['Z3', 'Z7']), (1, ['X4', 'X5', 'Y6', 'Y7']), (1, ['X4', 'Y5', 'Y6', 'X7']), (1, ['Y4', 'X5', 'X6', 'Y7']), (1, ['Y4', 'Y5', 'X6', 'X7']), (1, ['Z4']), (1, ['Z4', 'Z5']), (1, ['Z4', 'Z6']), (1, ['Z4', 'Z7']), (1, ['Z5']), (1, ['Z5', 'Z6']), (1, ['Z5', 'Z7']), (1, ['Z6']), (1, ['Z6', 'Z7']), (1, ['Z7'])]
  #nameList = ['ZZZZ','XXXX','XZXZ']
  nameList = ['ZZZZZZZZ', 'XXYYXXYY', 'YYXXYYXX', 'XYYXXYYX', 'YXXYYXXY', 'XZZZXZZZ', 'YZZZYZZZ', 'ZZZXZZZX', 'ZZZYZZZY', 'XXXYYYYX', 'XYYYYXXX', 'ZXZZZXZZ', 'ZYZZZYZZ', 'ZZXZZZXZ', 'ZZYZZZYZ', 'YXYXXYXY', 'YYXXXXYY', 'XZXZXXZX', 'XZXZXYZY', 'YZYXZXZX', 'YZYZYYZY', 'XXYZZZZY', 'XXZXZZXZ', 'XYYZZZZX', 'XYZYZZXZ', 'XZXXZXZZ', 'XZXYZYZZ', 'XZXZXZXZ', 'XZXZYZYZ', 'XZYZYZXZ', 'XZZXYZZY', 'XZZXZXXZ', 'XZZYYZZX', 'XZZYZYXZ', 'XZZZZXYY', 'XZZZZYYX', 'YXXZZZZY', 'YXZXZZYZ', 'YYXZZZZX', 'YYZYZZYZ', 'YZXZXZYZ', 'YZYYZYZZ', 'YZYZXZXZ', 'YZYZYZYZ', 'YZZXXZZY', 'YZZXZXYZ', 'YZZYXZZX', 'YZZYZYYZ', 'YZZZZXXY', 'YZZZZYXX', 'ZXXZXZZX', 'ZXYZYZZX', 'ZXZXXZXZ', 'ZXZXYZYZ', 'ZXZXZXZX', 'ZXZXZYZY', 'ZXZYZYZX', 'ZXZZXZXX', 'ZXZZYZYX', 'ZYXZXZZY', 'ZYYZYZZY', 'ZYZXZXZY', 'ZYZYXZXZ', 'ZYZYYZYZ', 'ZYZYZXZX', 'ZYZYZYZY', 'ZYZZXZXY', 'ZYZZYZYY', 'ZZYZYXZX']
  Nq = 8
  matching, termstrings = matchTermToCircuit(nameList, H, Nq)
  for key in matching:
    print('-----------------')
    matchedterms = matching[key]
    print('size = {}'.format(len(matchedterms)))
    terms = [termstrings[i] for i in matchedterms]
    print(terms)










