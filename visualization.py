'''
Teague Tomesh - 3/13/2019


'''

import matplotlib.pyplot as plt
import numpy as np





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




















