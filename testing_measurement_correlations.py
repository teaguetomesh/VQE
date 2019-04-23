'''
Teague Tomesh
04/23/2019
'''


import sys
import numpy as np
import matplotlib.pyplot as plt
import importlib
import math
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
import vqeTools
from qiskit import Aer, execute
from qiskit.providers.aer import QasmSimulator
import time


def state_histogram(inputDict, title):
    '''
    Print a histogram of measurement results produced during simulation of
    a circuit. The histogram of state counts is actually implemented with
    a bar graph in matplotlib.
    '''

    keys = sorted(inputDict)
    values = []
    for key in keys:
        values += [inputDict[key]] 

    fig, ax = plt.subplots(figsize=(10,5))
    x = np.arange(0,len(values))
    rects = ax.bar(x, values, width=1, align='center', color='lightblue', 
                   edgecolor='gray')

    
    plt.xticks(x, keys, rotation=45)
    # Pad margins so that markers don't get clipped by the axes
    plt.margins(y=0.2)
    # Tweak spacing to prevent clipping of tick-labels
    #plt.subplots_adjust(bottom=0.15)
    plt.title(title)
    plt.ylabel('Counts')
    #plt.xlim([0,bin.size])
    
    def autolabel(rs):
      """
      Attach a text label above each bar displaying its height
      """
      for r in rs:
        height = r.get_height()
        ax.text(r.get_x() + r.get_width()/2., 1.05*height, 
                '{}'.format(int(height)),
                ha='center', va='bottom')

    autolabel(rects)

    plt.show()
    plt.close()


def gen_single(H, Nq):
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

    circuitMatrix = []
    #print(quOpMatrix)
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
    #print(circuitMatrix)

    circuitList = []
    for m in range(numCircs):
        qr = QuantumRegister(Nq, name='qreg')
        cr = ClassicalRegister(Nq, name='creg')
        circ = QuantumCircuit(qr, cr)
        name = ''
        circ.barrier(qr)
        for n, qubit_row in enumerate(circuitMatrix):
            basis = qubit_row[m]
            if basis == 'X':
                circ.h(qr[n])
            elif basis == 'Y':
                circ.sdg(qr[n])
                circ.h(qr[n])
            name += basis
        circ.barrier(qr)
        for n in range(Nq):
            circ.measure(qr[n],cr[n])
        #print('m: ',m)
        #print('name: ',name)
        circuitList += [(circ, name)]
        #circ.draw(scale=0.8, filename='measure_{}_single_circ{}'.format(name,m), 
        #    output='mpl', plot_barriers=False, reverse_bits=True)

    return circuitList


def gen_triple(H, Nq):
    opDict = {'Z':0, 'X':1, 'Y':2}
    opList = ('Z','X','Y')

    quOpMatrix = np.zeros(shape=(Nq,3), dtype=int)
    for term in H:
        ops = term[1]
        for op in ops:
            if op[0] is 'I': continue
            quOpMatrix[int(op[1]), opDict[op[0]]] = 1

    # Set the quOpMatrix to the full, to simulate a more complicated H
    quOpMatrix = np.ones((Nq,3), dtype=int)

    numCircs = 0
    for row in quOpMatrix:
        rowSum = np.sum(row)
        if rowSum > numCircs:
            numCircs = rowSum  

    circuitMatrix = []
    #print(quOpMatrix)
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
    #print(circuitMatrix)

    circuitList = []
    for m in range(numCircs):
        qr = QuantumRegister(Nq, name='qreg')
        cr = ClassicalRegister(Nq, name='creg')
        circ = QuantumCircuit(qr, cr)
        name = ''
        circ.barrier(qr)
        for n, qubit_row in enumerate(circuitMatrix):
            basis = qubit_row[m]
            if basis == 'X':
                circ.h(qr[n])
            elif basis == 'Y':
                circ.sdg(qr[n])
                circ.h(qr[n])
            name += basis
        circ.barrier(qr)
        for n in range(Nq):
            circ.measure(qr[n],cr[n])
        #print('m: ',m)
        #print('name: ',name)
        circuitList += [(circ, name)]
        #circ.draw(scale=0.8, filename='measure_{}_triple_circ{}'.format(name,m), 
        #    output='mpl', plot_barriers=False, reverse_bits=True)

    return circuitList


def hamAvg_single(circList, H, Nq):
    nameList = [t[1] for t in circList]
    #matching, termStrings = matchTermToCircuit(nameList, H, Nq)

    shots   = 1000
    simulator = Aer.get_backend('qasm_simulator')
    totalE = 0
    for tup in circList:
      circ, name = tup
      #print(name)
      result = execute(circ, simulator).result()
      counts = result.get_counts(circ)
      #print(counts)
      
      #state_histogram(counts,'Measure {} Counts, {} shots'.format(name, shots))

      energySum = 0
      numShots = np.sum(list(counts.values()))

      coef = H[0][0]
      termString = 'XZYZ'
      runningSum = 0
      for state in counts:
        count = counts[state]
        
        ##############################
        # ComputContrib function #
        val = 1
        for c_ts, c_ss in zip(termString, state):
          if c_ts == '*':
            continue
          # measurements of 1 contribute factor -1
          if c_ss == '1':
            val = -1*val
        ##############################

        contrib = val
        #print(state, count, contrib)
        runningSum = runningSum + (contrib*count)
      
      # multiply by the term's coefficient
      energySum = coef*runningSum
      #print(coef, runningSum, energySum)

      avg_contrib = energySum / numShots

      totalE = totalE + avg_contrib

    return totalE


def hamAvg_triple(circList, H, Nq):
    nameList = [t[1] for t in circList]
    #matching, termStrings = matchTermToCircuit(nameList, H, Nq)

    shots   = 1000
    simulator = Aer.get_backend('qasm_simulator')
    totalE = 0
    count_dict = {}
    for tup in circList:
        circ, name = tup
        #print(name)
        result = execute(circ, simulator).result()
        counts = result.get_counts(circ)
        #print(counts)
        if name == 'ZZZZ':
            count_dict['Z'] = counts
        elif name == 'XXXX':
            count_dict['X'] = counts
        elif name == 'YYYY':
            count_dict['Y'] = counts
      
        #state_histogram(counts,'Measure {} Counts, {} shots'.format(name, shots))

    coef = H[0][0]
    termString = 'XZYZ'
    runningSum = 0
    term_sum = 0
    for i, b in enumerate(termString):
        #print(b)
        counts = count_dict[b]
        numShots = np.sum(list(counts.values()))

        ##############################
        # ComputContrib function #
        num_pos = 0
        num_neg = 0
        for state in counts:
            #print(state, counts[state])
            if state[i] == '0':
                num_pos += counts[state]
                #print('adding {} to num_pos'.format(counts[state]))
            elif state[i] == '1':
                num_neg += counts[state]
                #print('adding {} to num_neg'.format(counts[state]))
        ##############################

        cur_term = (num_pos - num_neg) / numShots
        term_sum = term_sum + cur_term

    avg_contrib = coef * term_sum
    totalE = totalE + avg_contrib

    return totalE


def main(hamiltonian):
    '''
    I need to test whether the measurement of a qubit in a particular basis is
    dependent upon the bases in which the neighboring qubits are measured.

    Example: given a term in a Hamiltonian, 'XZYZ', will the distribution look
    different if the average value of this term in measured with a single
    circuit measuring in the "XZYZ" bases, or with 3 different circuits
    measuring in the "ZZZZ", "XXXX", "YYYY" bases?
    '''

    start_time = time.time()

    # Initialize all the parameters as if we were measuring H = a*XZYZ
    # using the VQE algorithm with HartreeFock reference state, UCCSD_Whitfield
    # ansatz, and Nelder_Mead optimization.
    refStateModule = importlib.import_module('.HartreeFock', package="ReferenceState")
    ansatzModule = importlib.import_module('.UCCSD_4_Whitfield', package="Ansatz")
    optimizeModule = importlib.import_module('.Nelder_Mead', package="Optimizer")

    numQubits = 4
    numElectrons = 2
    
    # Generate a circuit for the reference state
    refCircuit = vqeTools.genRefState(refStateModule, numQubits, numElectrons)
    print('--reference circuit generated--')

    # Generate a circuit for the ansatz
    params = [math.pi/2 for i in range(7)]
    ansCircuit = vqeTools.genAnsatz(ansatzModule, numQubits, params)
    print('--ansatz circuit generated--')

    # Generate a circuit for measuring the energy
    # First, use a single circuit measuring in the "XZYZ" bases
    msrCircuit_single = gen_single(hamiltonian, numQubits)
    print('--measurement circuit (single) generated--')

    # Second, use three circuits measuring in the "ZZZZ", "XXXX", "YYYY" bases
    # and average their results
    msrCircuit_triple = gen_triple(hamiltonian, numQubits)
    print('--measurement circuit (triple) generated--')

    # Construct a full quantum circuit and investigate the distributions made
    # by the single and triple schemes
    cl_single = vqeTools.constructQuantumCircuit(refCircuit, ansCircuit, msrCircuit_single)
    cl_triple = vqeTools.constructQuantumCircuit(refCircuit, ansCircuit, msrCircuit_triple)
    print('--all circuits constructed--')

    single_sum = 0
    triple_sum = 0
    energies = []
    avg_energies = []
    for i in range(6000):
        print('ITERATION {}'.format(i+1))
        # Energy integration
        energy_single = hamAvg_single(cl_single, hamiltonian, numQubits)
        #print(energy_single)
        print('--finished averaging for single circuit--')
        energy_triple = hamAvg_triple(cl_triple, hamiltonian, numQubits)
        #print(energy_triple)
        print('--finished averaging for triple circuit--')

        single_sum += energy_single
        triple_sum += energy_triple

        energies += [(energy_single,energy_triple)]

        single_avg = single_sum / (i+1)
        triple_avg = triple_sum / (i+1)

        avg_energies += [(single_avg, triple_avg)]

    end_time = time.time()
    elapsed_time = end_time - start_time
    timestr = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    print('\nElapsed time: ', timestr)


    ###########################
    ######## PLOTTING #########
    ###########################

    avg_single = [t[0] for t in avg_energies]
    avg_triple = [t[1] for t in avg_energies]
    iteration = np.arange(1,len(avg_energies)+1)

    plt.figure(figsize=[10,10])
    plt.scatter(iteration, avg_single, s=10, c='green', label='single')
    plt.plot(iteration,avg_single,lw=1,c='green')
    plt.scatter(iteration, avg_triple, s=10, c='blue', label='triple')
    plt.plot(iteration,avg_triple,lw=1,c='blue')
    plt.legend()
    plt.xlabel('Iteration')
    plt.ylabel('Average reported energy')
    plt.title('Single vs. Triple measurement circuit average energies')

    plt.savefig('Debugging_Measurement/avg_energy_comparison.png')
    plt.close()

    single = [t[0] for t in energies]
    triple = [t[1] for t in energies]

    plt.figure(figsize=[10,10])
    plt.scatter(single, triple, s=10, c='green')
    plt.xlabel('Energy reported by single measure circuit')
    plt.ylabel('Energy reported by triple measure circuit')
    plt.title('Single vs. Triple reported energies')

    plt.savefig('Debugging_Measurement/direct_comparison.png')
    plt.close()



if __name__ == "__main__":
    H = [(-2, ['X0', 'Z1', 'Y2', 'Z3'])]
    main(H)
