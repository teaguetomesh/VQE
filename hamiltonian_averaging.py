'''



'''

from qiskit import Aer, execute
from qiskit.providers.aer import QasmSimulator
import visualization as vis


def run(circList, H):
    '''
    '''

    shots   = 1000
    simulator = Aer.get_backend('qasm_simulator')

    for circ in circList:
      result = execute(circ, simulator).result()
      counts = result.get_counts(circ)
      
      vis.state_histogram(counts,'Naive Ansatz Counts')

    return 'TESTRUNTESTRUN'