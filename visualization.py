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
    name = path[ir:-4]
    narr = name.split('.')
    name = narr[0][:-2] + narr[2][1:]
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
    savestr, molName = get_names(pecPath)
    r'$\alpha > \beta$'
    if molName == 'H2':
        mstr = r'$H_2$'
    titlestr = 'PEC for '+mstr

    print('~parsing terminal output~')
    extraData = parse_terminal_output(terPath)
    titlestr += ': {} & {}'.format(extraData[4][2].split('_')[0].capitalize(),
                                   extraData[4][3].replace('_',' '))

    fciData = np.genfromtxt('Hamiltonians/FCI_Energies/H2_sto-3g_FCI_energy.txt')

    # get the data from the file
    data = np.genfromtxt(pecPath)
    r = data[:,0]
    E = data[:,1]

    # Create the figure
    fig, ax = plt.subplots(figsize=(8,5))
    ax.axhline(y=0, color='k', lw=0.6)
    
    # Plot the data
    ax.plot(fciData[:,0],fciData[:,1],c='k',ls=':',label='FCI')
    ax.scatter(r,E,c='b',label='Final Energy')

    if optimize_layer:
        # Plot the initial ansatz guesses
        ax.scatter(extraData[0],extraData[1],s=6,c='r',label='Initial Energy')
        savestr += '_wOpt'

    if iteration_layer:
        for x, y, numi in zip(r, E, extraData[2]):
            ax.text(x, y-0.1,'{}'.format(numi), ha='center', va='top', 
                fontsize='small', color='dimgrey', fontstyle='oblique')
        ax.plot([],[],c='dimgrey',label='# Iterations')
        savestr += '_wItr'

    if fit_curve:
        # Fit the data with Lennard Jones potential
        #ax.plot(xx, pecFit, color='r')
        #popt, pcov = curve_fit(LennardJones, r, E, bounds=(0,100))
        #pecFit = [LennardJones(x,popt[0],popt[1]) for x in xx]
        #print(popt)
        xx = np.arange(0.5,2.5,0.01)
        yy = [LennardJones(x, 1, 0.75) for x in xx]
        plt.plot(xx,yy)
        
    ax.legend(loc='best')

    plt.xlim(0,3.5)
    plt.ylim(-1.5,2)

    plt.xlabel(r'Interatomic distance ($\AA$)')    
    plt.ylabel(r'$\langle H \rangle$ (Ha)')
    plt.title(titlestr)
    
    if save: 
        # Save figure to file
        savepath = 'Plots/'+savestr+'.png'
        plt.savefig(savepath)
        print('Figure saved to: ',savepath)
    
    plt.show()
    plt.close()


def plot_double():
    '''
    '''

    path1 = 'Results/H2_sto-3g_JW_0.1_to_3.0_wPrevParam1.txt'
    path2 = 'Results/H2_sto-3g_JW_0.1_to_3.0_try4.txt'
    fciPath = 'Hamiltonians/FCI_Energies/H2_sto-3g_FCI_energy.txt'

    savestr, molName = get_names(path1)
    savestr += '_vs_random_zoom'
    if molName == 'H2':
        mstr = r'$H_2$'
    titlestr = 'PEC for '+mstr+', Nelder-Mead, PrevParam vs Random ansatz'

    # get the data from the files
    data1 = np.genfromtxt(path1)
    r1 = data1[:,0]
    E1 = data1[:,1]
    data2 = np.genfromtxt(path2)
    r2 = data2[:,0]
    E2 = data2[:,1]
    rcolor = 'SkyBlue'
    pcolor = 'IndianRed'

    # get FCI energy from file
    fciData = np.genfromtxt(fciPath)
    fciR = fciData[:,0]
    fciE = fciData[:,1]

    # Create the figure
    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(8,8))
    pecAx = ax[0]
    errAx = ax[1]

    # Plot PEC
    pecAx.axhline(y=0, color='k', lw=0.6)
    
    # Plot the data
    pecAx.plot(fciR,fciE,c='k',ls=':',label='FCI')
    pecAx.scatter(r2,E2,c=rcolor,s=24,label='Random')
    #ax.plot(r2,E2,lw=1,c='b')
    pecAx.scatter(r1,E1,c=pcolor,s=24,label='PrevParam')
    #ax.plot(r1,E1,lw=1,c='r')

    pecAx.set_xlim(0,3.1)
    pecAx.set_ylim(-1.5,5)

    pecAx.set_ylabel(r'$\langle H \rangle$ (Ha)')
    pecAx.set_title(titlestr)
    pecAx.legend(loc='best')

    # Plot error with respect to FCI calculation
    errAx.axhline(y=0, color='k', lw=0.6)
    prevEDiff = [e1 - e2 for e1, e2 in zip(E1, fciE)]
    randEDiff = [e1 - e2 for e1, e2 in zip(E2, fciE)]
    width = 0.04
    rects1 = errAx.bar(r1 - width/2, prevEDiff, width, color=pcolor)
    rects2 = errAx.bar(r1 + width/2, randEDiff, width, color=rcolor)
    errAx.axhline(y=1.6e-3, color='k', lw=1, ls='-.', label='Chemical Accuracy')
    errAx.set_ylabel('Error wrt FCI energy (Ha)')
    errAx.set_xlabel(r'Interatomic distance ($\AA$)')
    errAx.legend(loc='upper right')
    errAx.set_ylim(0,0.1)

    
    # Save figure to file
    savepath = 'Plots/'+savestr+'.png'
    plt.savefig(savepath)
    print('Figure saved to: ',savepath)
    plt.show()
    plt.close()


# Main
def main(argv):

  fontweight = 'normal'
  plt.rcParams["font.weight"] = fontweight
  plt.rcParams["axes.labelweight"] = fontweight
  plt.rcParams["axes.titleweight"] = fontweight

  # 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'
  fontsize = 12
  plt.rcParams["font.size"] = fontsize
  plt.rcParams["axes.labelsize"] = fontsize
  plt.rcParams["axes.titlesize"] = fontsize
  # 'serif' | 'sans-serif' | 'cursive' | 'fantasy' | 'monospace'
  fontfamily = 'sans-serif'
  plt.rcParams["font.family"] = fontfamily

  optimize_layer = False
  iteration_layer = False
  save_image = False
  fit_curve = False
  vqe_name = ''
  plot_pec = False
  
  try:
   opts, args = getopt.getopt(argv,"dsiofn:",["name="])
  except getopt.GetoptError:
    print ('Usage: \n python visualization.py -dsiof -n <name>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print ('Usage: \n python visualization.py -dsiof -n <name>')
      sys.exit()
    elif opt in ("-n", "--name"):
      vqe_name = arg
      plot_pec = True
    elif opt in ("-i"):
      iteration_layer = True
    elif opt in ("-o"):
      optimize_layer = True
    elif opt in ("-s"):
      save_image = True
    elif opt in ("-f"):
      fit_curve = True
    elif opt in ("-d"):
      plot_double()
    
  if plot_pec:
    pecPath = 'Results/'+vqe_name
    ns = vqe_name.split('_')
    terminalOutputPath = 'Results/'+ns[0]+'_'+ns[1]+'_vqe_output_'+ns[-1]
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

















