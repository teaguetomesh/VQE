'''
Teague Tomesh - 3/13/2019

Given a circuit and Hamiltonian, run multiple simulations of that circuit to
determine the expected energy of the given Hamiltonian.

'''

from qiskit import Aer, execute
from qiskit.providers.aer import QasmSimulator
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

    print(H)

    termStrings = ['' for i in range(len(H))]
    print('empty TS: ',termStrings)
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

    print('termStrings: ',termStrings)

    matching = {}
    termIndices = np.arange(1,len(H))
    #print(termIndices)
    for circName in nameList:
        matchedTerms, delList = [], []
        #nTerms = len(termIndices)
        for i, t in enumerate(termIndices):
            ts = termStrings[t]
            if match(circName, ts):
                print('it matched!')
                matchedTerms += [t]
                delList += [i]
        matching[circName] = matchedTerms
        termIndices = [e for i, e in enumerate(termIndices) if i not in delList]
        print('del termIndices: ', termIndices)

    print('matching: ', matching)

    return matching, termStrings


def computeContrib(ts, ss):
    '''
    '''
    val = 1
    for c_ts, c_ss in zip(ts, ss):
        if c_ts == '*':
            continue
        # measurements of 1 contribute factor -1
        if c_ss == '1':
            val = -1*val
        # measurements of 0 contribute factor +1
        # so no need to multiply by +1
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

    print('these are the terms: ', terms)

    energySum = 0
    numShots = np.sum(list(counts.values()))

    for term in terms:
        coef = H[term][0]
        ts = termStrings[term]
        print('coef: {}'.format(coef))
        print('termString: {}'.format(ts))

        # for this term, average its energy contribution by counting the 
        # number of times it contributed +1 or -1
        runningSum = 0
        for state in counts:
            count = counts[state]
            contrib = computeContrib(ts, state)
            print(state, count, contrib)
            runningSum = runningSum + (contrib*count)

        # multiply by the term's coefficient
        energySum = energySum + (coef*runningSum)

    return (energySum / numShots)


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
    nameList = [t[1] for t in circList]
    matching, termStrings = matchTermToCircuit(nameList, H, Nq)

    #TODO: add logic to change the number of shots based on the desired 
    #      precision. Li et al: for desired precision, e, need O(w^2/e^2) shots.
    shots   = 1000
    simulator = Aer.get_backend('qasm_simulator')
    totalE = H[0][0]
    print(totalE)
    for tup in circList:
      circ = tup[0]
      name = tup[1]
      result = execute(circ, simulator).result()
      counts = result.get_counts(circ)
      print(counts)
      
      #vis.state_histogram(counts,'Naive Ansatz Counts: {} shots'.format(shots))

      totalE = totalE + energyIntegration(H, matching[name], counts, termStrings)

    return totalE











