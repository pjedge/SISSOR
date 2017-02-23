import matplotlib as mpl
mpl.use('Agg',force=True)
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/pedge/git/hapcut2/utilities')

import calculate_haplotype_statistics as chs

mpl.rc('legend', fontsize=13)
mpl.rc('xtick', labelsize=13)
mpl.rc('ytick', labelsize=13)
mpl.rc('axes', labelsize=13)
mpl.rc('axes', labelsize=13)
mpl.rcParams.update({'font.size': 13})
mpl.rc('lines', linewidth=3)

def plot_sissor(data,labels, outname):
    labels[3] = 'Basic Processing,\ncoverage >= 2'
    plt.figure(figsize=(10,5))
    ax = plt.subplot(121)
    N = 2
    errs0 = sum(data[0],chs.error_result())
    errs1 = sum(data[1],chs.error_result())
    errs2 = sum(data[2],chs.error_result())
    errs3 = sum(data[3],chs.error_result())

    ind = np.arange(N)  # the x locations for the groups
    width = 0.15       # the width of the bars

    plt.bar(ind, (errs0.get_switch_rate(),errs0.get_mismatch_rate()), color='k',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[0])

    plt.bar(ind+width, (errs1.get_switch_rate(),errs1.get_mismatch_rate()), color='r',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[1])

    plt.bar(ind+2*width, (errs2.get_switch_rate(),errs2.get_mismatch_rate()), color='b',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[2])

    plt.bar(ind+3*width, (errs3.get_switch_rate(),errs3.get_mismatch_rate()), color='y',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[3])


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

    plt.bar(ind, (errs0.get_N50()/1000000,errs0.get_AN50()/1000000), color='k',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[0])

    plt.bar(ind+width, (errs1.get_N50()/1000000,errs1.get_AN50()/1000000), color='r',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[1])

    plt.bar(ind+2*width, (errs2.get_N50()/1000000,errs2.get_AN50()/1000000), color='b',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[2])

    plt.bar(ind+3*width, (errs3.get_N50()/1000000,errs3.get_AN50()/1000000), color='y',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[3])

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




def plot_sissor_minimal(data,labels, outname):
    plt.figure(figsize=(10,5))
    ax = plt.subplot(121)
    N = 2
    data = [data[0],data[2]]
    labels = ['raw fragments','processed fragments']

    errs0 = sum(data[0],chs.error_result())
    errs1 = sum(data[1],chs.error_result())
    #errs2 = sum(data[2],chs.error_result())

    ind = np.arange(N)  # the x locations for the groups
    width = 0.3       # the width of the bars

    plt.bar(ind, (errs0.get_switch_rate(),errs0.get_mismatch_rate()), color='r',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=width,      # smaller bar width
            align='center',
            label=labels[0])

    plt.bar(ind+width, (errs1.get_switch_rate(),errs1.get_mismatch_rate()), color='b',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=width,      # smaller bar width
            align='center',
            label=labels[1])
    '''
    plt.bar(ind+2*width, (errs2.get_switch_rate(),errs2.get_mismatch_rate()), color='b',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[2])
    '''

    # add some text for labels, title and axes ticks
    #ax.set_ylim(0,0.025)
    ax.set_ylabel('Error Rate')
    ax.set_xticks(ind+1/2*width)
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
    width = 0.3       # the width of the bars

    plt.bar(ind, (errs0.get_N50()/1000000,errs0.get_AN50()/1000000), color='r',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=width,      # smaller bar width
            align='center',
            label=labels[0])

    plt.bar(ind+width, (errs1.get_N50()/1000000,errs1.get_AN50()/1000000), color='b',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=width,      # smaller bar width
            align='center',
            label=labels[1])
    '''
    plt.bar(ind+2*width, (errs2.get_N50()/1000000,errs2.get_AN50()/1000000), color='b',
            ecolor='black', # black error bar color
            alpha=0.5,      # transparency
            width=0.15,      # smaller bar width
            align='center',
            label=labels[2])
    '''

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Megabases')
    ax.set_xticks(ind+1/2*width)
    ax.set_xticklabels(('N50','AN50'))
    ax.yaxis.grid(True,color='grey')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    plt.legend(title='Haplotype assembly using:')
    plt.tight_layout()
    plt.savefig(outname)
