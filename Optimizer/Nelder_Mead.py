'''
Teague Tomesh - 3/13/2019

Implementation of the Nelder-Mead optimization algorithm for VQE
'''

import vqeTools
import random as rand
import math
import sys
import numpy as np


def amoebatry(simplexM, y, psum, ihi, fac, func):
    '''
    Perform the extrapolation by a factor, fac, through the face of the
    simplex opposite the highest scoring point. If the new point is better
    than the highest point, replace the high point.
    '''

    (npts, ndim) = simplexM.shape
    #print(npts, ndim)

    # construct the trial point
    ptry = np.zeros(ndim)
    fac1 = (1.0 - fac) / ndim
    fac2 = fac1 - fac

    ptry = [psum[j]*fac1 - simplexM[ihi,j]*fac2 for j in range(ndim)]
    
    # evaluate the trial point
    ytry = func(ptry)

    # compare
    if ytry < y[ihi]:
        y[ihi] = ytry
        for j in range(ndim):
            psum[j] += ptry[j] - simplexM[ihi,j]
        simplexM[ihi] = ptry

    return ytry


def nelder_mead(func, x_start,
                delta=0.1, tol=10e-6,
                no_improv_break=10, max_iter=5000):
    '''
    '''

    ###Initialize###
    ndim = len(x_start)
    npts = ndim + 1
    nfunc = 0
    niter = 0
    no_improv = 0

    # populate matrix holding points of the simplex
    simplexM = np.zeros((npts,ndim))
    for i in range(npts):
        for j in range(ndim):
            simplexM[i,j] = x_start[j]
        if (i != 0): simplexM[i,i-1] += delta

    y = np.zeros((npts))
    # compute func at each simplex point
    for i, x in enumerate(simplexM):
        score = func(x)
        nfunc += 1
        y[i] = score

    prev_best = y[0]
    psum = np.sum(simplexM, axis=0)

    ### Loop ###
    while True:
        ilo = 0
        # Determine indices of lowest, highest, 2nd highest scores
        if y[0] > y[1]:
            inhi = 1
            ihi  = 0
        else:
            inhi = 0
            ihi  = 1

        for i in range(npts):
            if (y[i] <= y[ilo]): ilo = i
            if (y[i] > y[ihi]):
                inhi = ihi
                ihi  = i
            elif (y[i] > y[inhi] and i != ihi): inhi = i

        best = y[ilo]
        
        # Check for max iterations
        if max_iter and niter >= max_iter:
            print('MAX ITERATIONS REACHED')
            return (simplexM[ilo], y[ilo])
        
        # Return if no improvment after no_improv_break iterations
        print('Current iteration: {}, Num function calls: {}, best so far: {}\
            '.format(niter, nfunc, best))

        if best < (prev_best - tol):
            no_improv = 0
            prev_best = best
        else:
            no_improv += 1

        if no_improv >= no_improv_break:
            print('REACHED MAX ITERATIONS WITH NO IMPROVEMENT')
            return (simplexM[ilo], y[ilo])

        # Begin a new iteration
        niter += 1

        # Reflect the simplex from the high point
        # Extrapolate by a factor -1 through the face of the simplex across
        # from the high point.
        ytry = amoebatry(simplexM, y, psum, ihi, -1.0, func)
        nfunc += 1
        if (ytry <= y[ilo]):
            # Got better result that our best point
            # Try another extrapolation with factor = 2
            ytry = amoebatry(p,y,psum,ihi,2.0,func)
            nfunc += 1
        elif (ytry >= y[inhi]):
            # The reflected point is worse than the 2nd worst point
            # Do a one-dimensional contraction
            ysave = y[ihi]
            ytry  = amoebatry(p,y,psum,ihi,0.5,func)
            nfunc += 1
            if (ytry >= ysave):
                # still not better
                for i in range(npts):
                    # Contract around the best point
                    if i != ilo:
                        for j in range(ndim):
                            val = 0.5*(simplexM[i,j]+simplexM[ilo,j])
                            psum[j] = val
                            simplexM[i,j] = val
                        y[i] = func(psum)
                        nfunc += 1
                # recompute psum
                psum = np.sum(simplexM, axis=0)
        # Return to beginning of loop
    return 0
    


# Minimize the measured energy
def minimizeEnergyObjective(hamiltonian, numQubits, ansatzModule, refCircuit, msrCircuit):
    '''
    Initialize parameters for and then call the Nelder-Mead optimization 
    function

        @param hamiltonian (List): List representation of Hamiltonian
        @param numQubits (int): number of qubits in the simulation
        @param ansatzModule (module): python module containing a script to
                        generate an ansatz circuit given a vector of parameters
        @param refCircuit, msrCircuit (QuantumCircuit): precomputed circuits
            which produce the reference state and the necessary measurements for
            this particular Hamiltonian.

        return: tuple (best parameter vector and best energy found)
    '''

    def f(params):
        '''
        Function handle which can be passed to the Nelder-Mead optimization
        algorithm. Compiles an ansatz circuit using the input parameters,
        combines the 3 different circuit modules together, then measures the
        expected energy.
        '''
        # Generate a circuit for the ansatz
        ansCircuit = vqeTools.genAnsatz(ansatzModule, numQubits, params)

        circList = vqeTools.constructQuantumCircuit(refCircuit, ansCircuit, msrCircuit)

        # Energy integration
        energy = vqeTools.hamiltonianAveraging(circList, hamiltonian, numQubits)
        #print(energy)
        
        return energy


    ### Start of Nelder-Mead simplex optimization ###
    initialparams = [rand.uniform(0,2*math.pi) for i in range(20)]
    final = nelder_mead(f, initialparams)

    return final






    












