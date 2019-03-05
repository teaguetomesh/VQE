#!/usr/bin/env python
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import qiskit
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit import execute, register, get_backend, compile
from qiskit.tools.visualization import plot_histogram
import random as rand
from scipy.optimize import minimize
import time
import math
import os


def genAnsatzCircuit(cName, p, group):
  # Create the circuit that implements the ansatz and performs a measurement
  # in the basis corresponding to its group

  qr = QuantumRegister(2, name='qreg')
  cr = ClassicalRegister(2, name='creg')

  vqeCircuit = QuantumCircuit(qr, cr, name=cName)

  # Apply gates to evaluate the ansatz
  vqeCircuit.rx(p[0], qr[0])
  vqeCircuit.rx(p[1], qr[1])
  vqeCircuit.rz(p[2], qr[0])
  vqeCircuit.rz(p[3], qr[1])
  vqeCircuit.h(qr[0])
  vqeCircuit.cx(qr[0], qr[1])
  vqeCircuit.h(qr[1])
  vqeCircuit.cx(qr[1], qr[0])
  vqeCircuit.rz(p[4], qr[0])
  vqeCircuit.rz(p[5], qr[1])
  vqeCircuit.rx(p[6], qr[0])
  vqeCircuit.rx(p[7], qr[1])
  vqeCircuit.rz(p[8], qr[0])
  vqeCircuit.rz(p[9], qr[1])

  # Rotate each qubit depending on the group such that the measurement 
  # is effectively a measurement in the X basis
  if (group == 2) or (group == 3):
    vqeCircuit.h(qr[0])
  if (group == 2) or (group == 4):
    vqeCircuit.h(qr[1])

  vqeCircuit.measure(qr[0], cr[0])
  vqeCircuit.measure(qr[1], cr[1])

  print("Generated circuit: {}".format(cName))

  return vqeCircuit


def printCircuit(circ):
  # Print out OPENQASM string
  print(circ.qasm())


#def drawCircuit(circ):
  # Draw the circuit
  #circuit_drawer(vqeCircuit)


def hamiltonianAveraging(amplitudes):
  # Evaluate the ansatz and find the expectation value for the energy
  #print("Available backends:", qiskit.Aer.backends())

  # Generate four circuits
  # The first is for group 1 (Z basis measurements) 
  circ1 = genAnsatzCircuit("VQE_isca_tutorial_1", amplitudes, group=1)
  # The second is for group 2 (X basis measurements)
  circ2 = genAnsatzCircuit("VQE_isca_tutorial_2", amplitudes, group=2)
  # The third is for group 3 ( Z on q1, X on q2)
  circ3 = genAnsatzCircuit("VQE_isca_tutorial_3", amplitudes, group=3)
  # The fourth is for group 4 (X on q1, Z on q2)
  circ4 = genAnsatzCircuit("VQE_isca_tutorial_4", amplitudes, group=4)

  circs = [circ1, circ2, circ3, circ4]

  backend = "qasm_simulator_py"
  shots   = 1000
  jobs, results, counts, circNames = [], [], [], []
  for circ in circs:
    newJob = execute(circ, backend=qiskit.Aer.get_backend(backend),
        shots=shots)
    jobs.append(newJob)
    newResult = newJob.result()
    results.append(newResult)
    counts.append(newResult.get_counts())
    circNames.append(newResult.get_names()[0])

  #plot_histogram(counts)

  # Save the state counts (Z-basis) for later plotting
  try:
    count00 = counts[0]['00']
  except:
    print('No 00 counts; set to 0')
    count00 = 0
  try:
    count01 = counts[0]['01']
  except:
    print('No 01 counts; set to 0')
    count01 = 0
  try:
    count10 = counts[0]['10']
  except:
    print('No 10 counts; set to 0')
    count10 = 0
  try:
    count11 = counts[0]['11']
  except:
    print('No 11 counts; set to 0')
    count11 = 0

  print("Counts: [00:{}], [01:{}], [10:{}], [11:{}]".format(count00, count01,
    count10, count11))

  # TODO: myCounts is only keeping track of the z-basis states: 00, 01, ...
  # Could also track the count of x-basis states: ++, +-, ...
  # from the circuits which measure the X terms in the hamiltonian
  myCounts = [count00, count01, count10, count11]

  # Find the expected energy by Hamiltonian averaging
  # GROUP 1

  # Term 1: I*I
  t1_obs_dict = {'00': 1, '01': 1, '10': 1, '11': 1}
  t1 = results[0].average_data(circNames[0], t1_obs_dict)
  a1 = -3.9734

  # Term 3: I*Z
  t3_obs_dict = {'00': 1, '01': -1, '10': 1, '11': -1}
  t3 = results[0].average_data(circNames[0], t3_obs_dict)
  a3 = -1.0052

  # Term 7: Z*I
  t7_obs_dict = {'00': 1, '01': 1, '10': -1, '11': -1}
  t7 = results[0].average_data(circNames[0], t7_obs_dict)
  a7 = -1.0052

  # Term 9: Z*Z
  t9_obs_dict = {'00': 1, '01': -1, '10': -1, '11': 1}
  t9 = results[0].average_data(circNames[0], t9_obs_dict)
  a9 = 0.2779

  # GROUP 2

  # Term 2: I*X
  t2_obs_dict = {'00': 1, '01': -1, '10': 1, '11': -1}
  t2 = results[1].average_data(circNames[1], t2_obs_dict)
  a2 = -0.2385

  # Term 4: X*I
  t4_obs_dict = {'00': 1, '01': 1, '10': -1, '11': -1}
  t4 = results[1].average_data(circNames[1], t4_obs_dict)
  a4 = -0.2385

  # Term 5: X*X
  t5_obs_dict = {'00': 1, '01': -1, '10': -1, '11': 1}
  t5 = results[1].average_data(circNames[1], t5_obs_dict)
  a5 = 0.2343

  # GROUP 3

  # Term 6: X*Z
  t6_obs_dict = {'00': 1, '01': -1, '10': -1, '11': 1}
  t6 = results[2].average_data(circNames[2], t6_obs_dict)
  a6 = 0.2385

  # GROUP 4

  # Term 8: Z*X
  t8_obs_dict = {'00': 1, '01': -1, '10': -1, '11': 1}
  t8 = results[3].average_data(circNames[3], t8_obs_dict)
  a8 = 0.2385

  #expected_value = (count00 + count01 - count10 - count11) / shots
  #truth = 'True' if (t7avgVal == expected_value) else 'False'
  #print('Do the avg vals match: {}'.format(truth))

  avgE = a1*t1 + a2*t2 + a3*t3 + a4*t4 + a5*t5 + a6*t6 + a7*t7 + a8*t8 + a9*t9

  return avgE, myCounts


def hamiltonianAveraging(amplitudes, alphas):
  # Evaluate the ansatz and find the expectation value for the energy
  #print("Available backends:", qiskit.Aer.backends())

  # Generate four circuits
  # The first is for group 1 (Z basis measurements) 
  circ1 = genAnsatzCircuit("VQE_isca_tutorial_1", amplitudes, group=1)
  # The second is for group 2 (X basis measurements)
  circ2 = genAnsatzCircuit("VQE_isca_tutorial_2", amplitudes, group=2)
  # The third is for group 3 ( Z on q1, X on q2)
  circ3 = genAnsatzCircuit("VQE_isca_tutorial_3", amplitudes, group=3)
  # The fourth is for group 4 (X on q1, Z on q2)
  circ4 = genAnsatzCircuit("VQE_isca_tutorial_4", amplitudes, group=4)

  circs = [circ1, circ2, circ3, circ4]

  backend = "qasm_simulator_py"
  shots   = 1000
  jobs, results, counts, circNames = [], [], [], []
  for circ in circs:
    newJob = execute(circ, backend=qiskit.Aer.get_backend(backend),
        shots=shots)
    jobs.append(newJob)
    newResult = newJob.result()
    results.append(newResult)
    counts.append(newResult.get_counts())
    circNames.append(newResult.get_names()[0])

  #plot_histogram(counts)

  # Save the state counts (Z-basis) for later plotting
  try:
    count00 = counts[0]['00']
  except:
    print('No 00 counts; set to 0')
    count00 = 0
  try:
    count01 = counts[0]['01']
  except:
    print('No 01 counts; set to 0')
    count01 = 0
  try:
    count10 = counts[0]['10']
  except:
    print('No 10 counts; set to 0')
    count10 = 0
  try:
    count11 = counts[0]['11']
  except:
    print('No 11 counts; set to 0')
    count11 = 0

  print("Counts: [00:{}], [01:{}], [10:{}], [11:{}]".format(count00, count01,
    count10, count11))

  myCounts = [count00, count01, count10, count11]

  # Find the expected energy by Hamiltonian averaging
  # GROUP 1

  # Term 1: I*I
  t1_obs_dict = {'00': 1, '01': 1, '10': 1, '11': 1}
  t1 = results[0].average_data(circNames[0], t1_obs_dict)
  a1 = alphas[0]

  # Term 3: I*Z
  t3_obs_dict = {'00': 1, '01': -1, '10': 1, '11': -1}
  t3 = results[0].average_data(circNames[0], t3_obs_dict)
  a3 = alphas[2]

  # Term 7: Z*I
  t7_obs_dict = {'00': 1, '01': 1, '10': -1, '11': -1}
  t7 = results[0].average_data(circNames[0], t7_obs_dict)
  a7 = alphas[6]

  # Term 9: Z*Z
  t9_obs_dict = {'00': 1, '01': -1, '10': -1, '11': 1}
  t9 = results[0].average_data(circNames[0], t9_obs_dict)
  a9 = alphas[8]

  # GROUP 2

  # Term 2: I*X
  t2_obs_dict = {'00': 1, '01': -1, '10': 1, '11': -1}
  t2 = results[1].average_data(circNames[1], t2_obs_dict)
  a2 = alphas[1]

  # Term 4: X*I
  t4_obs_dict = {'00': 1, '01': 1, '10': -1, '11': -1}
  t4 = results[1].average_data(circNames[1], t4_obs_dict)
  a4 = alphas[3]

  # Term 5: X*X
  t5_obs_dict = {'00': 1, '01': -1, '10': -1, '11': 1}
  t5 = results[1].average_data(circNames[1], t5_obs_dict)
  a5 = alphas[4]

  # GROUP 3

  # Term 6: X*Z
  t6_obs_dict = {'00': 1, '01': -1, '10': -1, '11': 1}
  t6 = results[2].average_data(circNames[2], t6_obs_dict)
  a6 = alphas[5]

  # GROUP 4

  # Term 8: Z*X
  t8_obs_dict = {'00': 1, '01': -1, '10': -1, '11': 1}
  t8 = results[3].average_data(circNames[3], t8_obs_dict)
  a8 = alphas[7]

  #expected_value = (count00 + count01 - count10 - count11) / shots
  #truth = 'True' if (t7avgVal == expected_value) else 'False'
  #print('Do the avg vals match: {}'.format(truth))

  avgE = a1*t1 + a2*t2 + a3*t3 + a4*t4 + a5*t5 + a6*t6 + a7*t7 + a8*t8 + a9*t9

  return avgE, myCounts



def energy_objective(amplitudes):
  ''' Evaluate the expected energy 
  Args:
    amplitudes(ndarray): stores the thetas of the gate rotations

  Returns:
    energy(float): Expected energy of the Hamiltonian for the given angles
  '''

  # Evaluate the expected energy
  energy, counts = hamiltonianAveraging(amplitudes)
  #energy = -1 * energy

  # Do some bookkeeping on current counts / energies
  # This code moved here from optCallback because depending on the 
  # optimization algorithm, optCallback may be called 3 times while
  # energy_objective may be called a few hundred times
  global iterationCount
  iterationCount = iterationCount + 1

  dataList = [iterationCount, energy]
  dataList.extend(counts)

  fn = 'vqe_run_data/vqe_opt_out.txt'
  with open(fn, 'a') as sf:
    sf.write('{0:<6d} {1:>7.3f} {2:>6d} {3:>6d} {4:>6d} {5:>6d}\n'.format(
      iterationCount, energy, counts[0], counts[1], counts[2],
      counts[3]))

    print("Wrote to {}".format(fn))

  # return the expected energy <psi|H|psi> for the given amplitudes
  return energy


def energy_objective_PEC(amplitudes, coefs):
  ''' Evaluate the expected energy 
  Args:
    amplitudes(ndarray): stores the thetas of the gate rotations
    coefs(ndarray): stores the current distance, R, and term coefficients

  Returns:
    energy(float): Expected energy of the Hamiltonian for the given angles
  '''

  # Evaluate the expected energy
  alphas = coefs[1:]
  energy, counts = hamiltonianAveraging(amplitudes, alphas)

  # return the expected energy <psi|H|psi> for the given amplitudes
  return energy


def optCallback(curParams):
  ''' Perform some bookkeeping during each iteration '''
  energy, counts = hamiltonianAveraging()
  energy = -1 * energy

  global iterationCount
  iterationCount = iterationCount + 1

  dataList = [iterationCount, energy]
  dataList.extend(counts)

  fn = 'vqe_run_data/vqe_opt_out.txt'
  with open(fn, 'a') as sf:
    sf.write('{0:<6d} {1:>7.3f} {2:>6d} {3:>6d} {4:>6d} {5:>6d}\n'.format(
      iterationCount, energy, counts[0], counts[1], counts[2],
      counts[3]))

    print("Wrote to {}".format(fn))

  return False


def makeFig():
  ''' Plot the results of the VQE algorithm '''
  matplotlib.rcParams.update({'font.size': 16})

  fn = 'vqe_run_data/vqe_opt_out.txt'
  data = np.genfromtxt(fn)

  iteration = data[:,0]
  energies  = data[:,1]
  count00   = data[:,2]
  count01   = data[:,3]
  count10   = data[:,4]
  count11   = data[:,5]

  fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=[16,12])
  ax1.plot(iteration, count00, label='00')
  ax1.plot(iteration, count01, label='01')
  ax1.plot(iteration, count10, label='10')
  ax1.plot(iteration, count11, label='11')
  ax1.legend(loc='upper right')
  titleStr = opt_method + ' VQE, H=II+IX+IZ+XI+XX+XZ+ZI+ZX+ZZ'
  ax1.set_title(titleStr)
  ax1.set_ylabel('Counts')
  ax2.plot(iteration, energies)
  ax2.set_ylabel('Energy')
  ax2.set_xlabel('Iteration')

  ax2Str = 'Opt Energy: {0:.3f}\nFunc evals: {1:d}\nRuntime: {2:.3f}s'.format(
      opt_energy, opt_result.nfev, elapsed_time)
  ax2.text(0.99,0.98, ax2Str, horizontalalignment='right', 
      verticalalignment='top', transform=ax2.transAxes, 
      bbox=dict(facecolor='white', alpha=0.8))
  plt.savefig('vqe_run_data/vqe_opt_out.png')
  plt.show()
  plt.close()


def makePEC(amplitudes, opt_method, coef_fn):
  ''' Produce a full potential energy curve
  Args:
    amplitudes(ndarray): stores the initial thetas for the gate rotations
    opt_method(str): which scipy optimization method to use
    coef_fn(str): filename of the CSV containing the h_i and h_ij
            coefficients at each distance, R

  Returns:
    Nothing, but writes to vqe_run_data
  '''

  # Load the coefficient table
  coefTable = np.loadtxt(coef_fn, delimiter=',')
  nrows, ncols = coefTable.shape

  # Loop calls to scipy minimize
  for i in range(nrows):
    curCoefs = coefTable[i,:]
    opt_result = minimize(energy_objective_PEC, amplitudes, args=(curCoefs),
        method=opt_method, callback=None,
        options={'disp':True})

    fn = 'vqe_run_data/vqe_PEC_out.txt'
    with open(fn, 'a') as sf:
      sf.write('{0:<7.3} {1:>7.3f}\n'.format(curCoefs[0], opt_result.fun))

    print("Wrote to {}".format(fn))

  # Produce a figure of the PEC


#### Begin Main ####

# My api token
APItoken = '762d784bc177a138205d2d7cb326a4521f2f3f7aea8499dde \
    209db7842a18feda8b632434ca2213bf39b08897a351638d010f6da16b50c46e097570377ff4fa6'
qx_config = {
    "APItoken": APItoken,
    "url":"https://quantumexperience.ng.bluemix.net/api"}

start_time = time.time()

# Ansatz Parameters
alpha0   = rand.uniform(0,2*math.pi)
alpha1   = rand.uniform(0,2*math.pi)
beta0    = rand.uniform(0,2*math.pi)
beta1    = rand.uniform(0,2*math.pi)
gamma0   = rand.uniform(0,2*math.pi)
gamma1   = rand.uniform(0,2*math.pi)
delta0   = rand.uniform(0,2*math.pi)
delta1   = rand.uniform(0,2*math.pi)
epsilon0 = rand.uniform(0,2*math.pi)
epsilon1 = rand.uniform(0,2*math.pi)

'''
alpha0   = 0.0
alpha1   = 0.0
beta0    = 0.0
beta1    = 0.0
gamma0   = 0.0
gamma1   = 0.0
delta0   = 0.0
delta1   = 0.0
epsilon0 = 0.0
epsilon1 = 0.0
'''

plist    = [alpha0, alpha1, beta0, beta1, gamma0, gamma1, delta0, delta1,
    epsilon0, epsilon1]

iterationCount = 0

initial_amps = np.ndarray((len(plist),), buffer=np.array(plist))
print("Initial thetas: {}".format(initial_amps))

#initial_energy = energy_objective(initial_amps)

# Check if a data file already exists
fn = 'vqe_run_data/vqe_opt_out.txt'
if os.path.isfile(fn):
  os.rename(fn, 'vqe_run_data/vqe_opt_out_old.txt')

# Produce a full potential energy curve
opt_method = 'Powell'
coef_fn    = 'vqe_run_data/peruzzo_supp_table2.csv'
makePEC(initial_amps, opt_method, coef_fn)

exit()

# Run scipy optimization to find optimal thetas
opt_method = 'Nelder-Mead'
opt_result = minimize(energy_objective, initial_amps,
    method=opt_method, callback=None,
    options={'disp':True})


opt_energy, opt_amps = opt_result.fun, opt_result.x
print("\nOptimal Energy: {}".format(opt_energy))
print("Optimal Thetas: {}".format(opt_amps))
#print("Initial Energy:{}".format(initial_energy))

elapsed_time = time.time() - start_time

print("Runtime: {}".format(elapsed_time))

# Create figure showing counts / energy over iterations
makeFig()

















