#!/usr/bin/env python

'''

TEAGUE TOMESH
04/03/2019

Python script for calling multiple instances of main.py

'''

import sys
import subprocess


def main(argv):
    '''
    argv contains the arguments to use when calling main.py as well as an int indicating the
    number of trials to run.
    '''
    num_args = 9

    if len(argv) != num_args:
        print('\n---ERROR---')
        print('Must pass {} arguments to gen_vqe.py:'.format(num_args))
        print('  <hamiltonian> <referenceState> <ansatz> <numQubits> <optimizer>',
              '<pec_output> <terminal_output> <profile_output> <numIterations>\n')
        sys.exit(2)

    hamPath   = argv[0]
    refState  = argv[1]
    ansatz    = argv[2]
    numQ      = argv[3]
    optimizer = argv[4]
    pec_output  = argv[5]
    term_output = argv[6]
    prof_output = argv[7]
    num_iters = int(argv[8])

    profile = True

    for i in range(num_iters):
        term_output_file = 'Results/'+term_output+'{}.txt'.format(i)
        print('Current Iteration: {}'.format(i))
        print('---starting subprocess---')
        with open(term_output_file,'w') as tf:
            outputPath = pec_output + str(i)
            if profile:
                profile_fn = 'Results/'+prof_output+'{}.txt'.format(i)
                subprocess.run(['python','-m','cProfile','-o',profile_fn,'main.py',
                       '-p',hamPath,'-r',refState,'-a',ansatz,'-q',numQ,'-o',
                       optimizer,'-t',outputPath],stdout=tf)

            else:
                subprocess.run(['./main.py','-p',hamPath,'-r',refState,
                       '-a',ansatz,'-q',numQ,'-o',optimizer,'-t',outputPath],
                       stdout=tf)
        print('Finished subprocess {}'.format(i))

if __name__ == "__main__":
    main(sys.argv[1:])













