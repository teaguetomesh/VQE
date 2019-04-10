#!/usr/bin/env python

'''

TEAGUE TOMESH
04/02/2019

Python script for automating calls to visualization.py

'''

import subprocess
import sys
import glob

def main(argv):
    '''
    '''
  pec_path = ''
  save_path = ''
  ter_path = ''
  
  try:
   opts, args = getopt.getopt(argv,"p:s:t:",["pec_path=","save_path=",
                              "ter_path"])
  except getopt.GetoptError:
    print ('Usage: \n ./gen_plots.py -p <pec_path> -s <save_path> -t <ter_path>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print ('Usage: \n ./gen_plots.py -p <pec_path> -s <save_path> -t <ter_path>')
      sys.exit()
    elif opt in ("-p", "--pec_path"):
      pec_path = arg
    elif opt in ("-s", "--save_path"):
      save_image = True
      save_path = arg
    elif opt in ("-t", "--ter_path"):
      ter_path = arg

    vis = 'python visualization.py'
    

    #plot 1
    #command = vis + '-e -p ' + pec_path + ' -s ' + save_path + ' -t ' + ter_path
    #subprocess.run(command.split())

    #plot 2
    #command = vis + '-eio -p ' + pec_path + ' -s ' + save_path + ' -t ' + ter_path
    #subprocess.run(command.split())




if __name__ == "__main__":
  main(sys.argv[1:])