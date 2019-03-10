'''


'''

import matplotlib.pyplot as plt
import numpy as np


def state_histogram(inputDict, title):
    '''
    '''

    keys = inputDict.keys()
    values = inputDict.values()

    fig, ax = plt.subplots()
    counts, bins, patches = ax.hist(values, facecolor='green', edgecolor='gray')

    # Set the ticks to be at the edges of the bins.
    ax.set_xticks(bins)
    
    # Label the bins
    bin_centers = 0.5 * np.diff(bins) + bins[:-1]
    for key, x in zip(keys, bin_centers):
      # Label the bins
      ax.annotate(key, xy=(x, 0), xycoords=('data', 'axes fraction'),
        xytext=(0, -18), textcoords='offset points', va='top', ha='center')

    ax.get_xaxis().set_visible(False)

    # Give ourselves some more room at the bottom of the plot
    plt.subplots_adjust(bottom=0.15)
    plt.show()