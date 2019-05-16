#!/usr/bin/env python

'''

TEAGUE TOMESH
04/03/2019

Python script for calling multiple instances of main.py

'''

import sys
import subprocess
import time


def main(argv):
    '''
    argv contains the arguments to use when calling main.py as well as an int indicating the
    number of trials to run.
    '''
    num_args = 9

    '''
    if len(argv) != num_args:
        print('\n---ERROR---')
        print('Must pass {} arguments to gen_vqe.py:'.format(num_args))
        print('  <hamiltonian> <referenceState> <ansatz> <numQubits> <optimizer>',
              '<pec_output> <terminal_output> <profile_output> <numIterations>\n')
        sys.exit(2)
    '''
    #hamPath   = argv[0]
    #refState  = argv[1]
    refState  = 'HartreeFock'
    #ansatz    = argv[2]
    ansatz    = 'UCCSD_Whitfield'
    #numQ      = argv[3]
    #optimizer = argv[4]
    optimizer = 'Nelder_Mead'
    #pec_output  = argv[5]
    #term_output = argv[6]
    #prof_output = argv[7]
    #num_iters = int(argv[8])

    profile = False

    start_time = time.time()

    # Modify this file to iterate through all the hamiltonians generated for each 
    # occupied space and active space combination - EDIT: skipping AS4 because currently intractable
    #for i in range(num_iters):
    for os in range(5):
        for AS in range(1,4):
            hamPath = 'Hamiltonians/H2_6-31g_JW_OS{}/AS{}/'.format(os,AS)
            #term_output_file = 'Results/'+term_output+'{}.txt'.format(i)
            term_output_file = 'Results/UCCSD_Whitfield_6-31g/OS{}_AS{}/terminal_output.txt'.format(os,AS)
            numQ = str(AS * 2)
            print('Current Iteration: OS{} AS{}'.format(os,AS))
            print('---starting subprocess---')
            with open(term_output_file,'w') as tf:
                #outputPath = pec_output + str(i)
                outputPath = 'UCCSD_Whitfield_6-31g/OS{}_AS{}/pec_output.txt'.format(os,AS)
                if profile:
                    #profile_fn = 'Results/'+prof_output+'{}.txt'.format(i)
                    profile_fn = 'Results/UCCSD_Whitfield_6-31g/OS{}_AS{}/profile_output.txt'.format(os,AS)
                    subprocess.run(['python','-m','cProfile','-o',profile_fn,'main.py',
                           '-p',hamPath,'-r',refState,'-a',ansatz,'-q',numQ,'-o',
                           optimizer,'-t',outputPath],stdout=tf)

                else:
                    subprocess.run(['./main.py','-p',hamPath,'-r',refState,
                           '-a',ansatz,'-q',numQ,'-o',optimizer,'-t',outputPath],
                           stdout=tf)
                print('Finished OS{} AS{}'.format(os,AS))
                cur_time = time.time()
                elapsed_time = cur_time - start_time
                timestr = time.strftime("%d:%H:%M:%S", time.gmtime(elapsed_time))
                print('Elapsed time (D:H:M:S): ',timestr)

if __name__ == "__main__":
    main(sys.argv[1:])













