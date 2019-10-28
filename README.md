# VariationalQuantumEigensolver
A modular implementation of the VQE algorithm for ground state estimation

# Required Packages

python3.6

Qiskit python package

    pip install qiskit

OpenFermion python package

    pip install openfermion

# Getting Started
The VQE algorithm is run with the main.py file
```
python main.py -h
```
Will provide a list of arguments that can be supplied.

An example execution:
```
python main.py --hamiltonian Hamiltonians/H2_6-31g_JW_OS1/AS3/qubitH_H2_6-31g_JW_0.1.txt --refstate HartreeFock --ansatz UCCSD_Whitfield --qubits 6 --optimizer Nelder_Mead --output Results/output_file.txt
```
