import matplotlib as mpl
mpl.use('Agg',force=True)
import numpy as np
import matplotlib.pyplot as plt
import sys
import calculate_haplotype_statistics as chs

mpl.rc('legend', fontsize=13)
mpl.rc('xtick', labelsize=13)
mpl.rc('ytick', labelsize=13)
mpl.rc('axes', labelsize=13)
mpl.rc('axes', labelsize=13)
mpl.rcParams.update({'font.size': 13})
mpl.rc('lines', linewidth=3)

# this function takes data in the form of error_result objects from the calculate_haplotype_statistics class,
# and creates a neat bar chart of genomewide haplotype error rates and completeness metrics.

def plot_sissor(data,labels, outname):
    plt.figure(figsize=(10,5))
    ax = plt.subplot(121)
    N = 2
    errs0 = data[0]
    errs1 = data[1]

    ind = np.arange(N)  # the x locations for the groups
    width = 0.15       # the width of the bars

    plt.bar(ind, (errs0.get_switch_rate(),errs0.get_mismatch_rate()), color='r',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[0])

    plt.bar(ind+width, (errs1.get_switch_rate(),errs1.get_mismatch_rate()), color='b',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[1])

    # add some text for labels, title and axes ticks
    #ax.set_ylim(0,0.025)
    ax.set_ylabel('Error Rate')
    ax.set_xticks(ind+3/2*width)
    ax.set_xticklabels(('Switch','Mismatch'))
    ax.yaxis.grid(True,color='grey')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")

    ax = plt.subplot(122)
    N = 2

    ind = np.arange(N)  # the x locations for the groups
    width = 0.15       # the width of the bars

    plt.bar(ind, (errs0.get_N50()/1000000,errs0.get_AN50()/1000000), color='r',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[0])

    plt.bar(ind+width, (errs1.get_N50()/1000000,errs1.get_AN50()/1000000), color='b',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[1])

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Megabases')
    ax.set_xticks(ind+3/2*width)
    ax.set_xticklabels(('N50','AN50'))
    ax.yaxis.grid(True,color='grey')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outname)
