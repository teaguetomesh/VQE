#!/usr/bin/env python

'''

TEAGUE TOMESH
04/02/2019

Python script for automating calls to visualization.py

'''

import subprocess
import sys

def main(argv):
    '''
    '''
    
    dataPath = argv[0]

    vis = 'python visualization.py -'
    
    #plot 1
    command = vis + 'sn ' + dataPath
    subprocess.run(command.split())

    #plot 2
    command = vis + 'soin ' + dataPath
    subprocess.run(command.split())




if __name__ == "__main__":
  main(sys.argv[1:])