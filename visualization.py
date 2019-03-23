'''
Teague Tomesh - 3/13/2019


'''

import matplotlib.pyplot as plt
import numpy as np
import sys
import getopt
from scipy.optimize import curve_fit
import math





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


def parse_terminal_output(path):
    '''
    '''
    info, R, startPts, iterations, funcCalls = [],[],[],[],[]
    with open(path, 'r') as tf:
        for i, line in enumerate(tf):
            if line == '\n': continue
            if i < 5:
                l = line.split(':')[1].split()[0]
                info += [l]
                continue
            #print(line)
            if 'Begin VQE' in line:
                R += [float(line.split(':')[1])]
            if 'Current iteration: 0' in line:
                startE = float(line.split(':')[3])
                startPts += [startE]
            if 'MAX ITERATIONS' in line:
                numIter = int(prevLine.split(':')[1].split(',')[0])
                numFunc = int(prevLine.split(':')[2].split(',')[0])
                iterations += [numIter]
                funcCalls += [numFunc]

            prevLine = line

    return(R, startPts, iterations, funcCalls, info)


def get_names(path):
    '''
    '''
    for i, c in enumerate(path):
        if c == '/':
            ir = i+1
            break
    name = path[ir:-4]+'_PEC'
    print(name)
    molName = ''
    for c in name:
        if c == '_':
            break
        else:
            molName += c
    return name, molName


def LennardJones(r, epsilon, r_m):
    '''
    Use the Lennard-Jones potential to fit the PES curves.

        V = ep((r_m/r)^12-2(r_m/r)^6)

    '''
    #ror = r_m / r
    #t12 = math.pow(ror,12)
    #t6  = math.pow(ror,6)
    #return epsilon*(t12 - 2*t6)
    return epsilon * ((r_m/r)**12 - 2*(r_m/r)**6)


def plotPEC(options, pecPath, terPath):
    '''
    '''
    save = options[0]
    optimize_layer = options[1]
    iteration_layer = options[2]
    fit_curve = options[3]

    # get molecule name and file name from the path
    name, molName = get_names(pecPath)
    titlestr = 'PEC for {}'.format(molName)

    if optimize_layer or iteration_layer:
        print('~parsing terminal output~')
        extraData = parse_terminal_output(terPath)
        titlestr += ' with {} optimization'.format(extraData[4][3])
    else:
        extraData = []

    # get the data from the file
    data = np.genfromtxt(pecPath)
    r = data[:,0]
    E = data[:,1]

    # Create the figure
    fig, ax = plt.subplots(figsize=(8,5))
    ax.axhline(y=0, color='k', lw=0.6)
    
    # Plot the data
    ax.scatter(r,E,c='b')

    if optimize_layer:
        # Plot the initial ansatz guesses
        ax.scatter(extraData[0],extraData[1],s=4,c='r')

    if fit_curve:
        # Fit the data with Lennard Jones potential
        #ax.plot(xx, pecFit, color='r')
        #popt, pcov = curve_fit(LennardJones, r, E, bounds=(0,100))
        #pecFit = [LennardJones(x,popt[0],popt[1]) for x in xx]
        #print(popt)
        xx = np.arange(0.5,2.5,0.01)
        yy = [LennardJones(x, 1, 0.75) for x in xx]
        plt.plot(xx,yy)
        

    plt.xlim(0,3.5)
    plt.ylim(-1.5,2)

    plt.xlabel('Interatomic distance (angstroms)')
    plt.ylabel('Energy (Ha)')
    plt.title(titlestr)
    
    if save: 
        # Save figure to file
        savepath = 'Plots/'+name+'.png'
        plt.savefig(savepath)
        print('Figure saved to: ',savepath)
    
    plt.show()
    plt.close()


# Main
def main(argv):

  pecPath = ''
  terminalOutputPath = ''
  optimize_layer = False
  iteration_layer = False
  save_image = False
  fit_curve = False
  
  try:
   opts, args = getopt.getopt(argv,"siofp:t:",["pec=","terminal="])
  except getopt.GetoptError:
    print ('Usage: \n python visualization.py -p <pec> -t <terminal>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print ('Usage: \n python visualization.py -p <pec> -t <terminal>')
      sys.exit()
    elif opt in ("-p", "--pec"):
      pecPath = arg
    elif opt in ("-t", "--terminal"):
      terminalOutputPath = arg
    elif opt in ("-i"):
      iteration_layer = True
    elif opt in ("-o"):
      optimize_layer = True
    elif opt in ("-s"):
      save_image = True
    elif opt in ("-f"):
      fit_curve = True
    

  print('\nPEC: ',pecPath,'\nTerminal Output: ',terminalOutputPath,'\n')

  options = [save_image, optimize_layer, iteration_layer, fit_curve]
  print('Options:')
  print(' Save: {}'.format(options[0]))
  print(' Optimize_layer: {}'.format(options[1]))
  print(' Iteration_layer: {}'.format(options[2]))
  print(' Fit_curve: {}'.format(options[3]))

  plotPEC(options, pecPath, terminalOutputPath)

  

if __name__ == "__main__":
  main(sys.argv[1:])

















