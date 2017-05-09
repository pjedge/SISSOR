import pickle
import sys
from collections import defaultdict
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np

def generate_phased_calls_table(output_file, same_cell_counts, all_cell_counts, cross_cell_counts, ind_same_cell_counts):

    line0 = ['']
    line1 = ['total pos']
    line2 = ['overlap CGI']
    line3 = ['total SNPs']
    line4 = ['different call from CGI/WGS']
    line5 = ['different calls from CGI/WGS+BAC']
    line6 = ['different calls from CGI/WGS+BAC+third_chamber']
    line7 = ['error rate upper bound']
    line8 = ['FDR']

    for mode,f in [('same_cell',same_cell_counts),('all_cell',all_cell_counts),('cross_cell',cross_cell_counts)]:

        line0.append(mode)
        counts = pickle.load(open(f,'rb'))
        err = sum([x[1] for x in counts.items() if x[0].matches_CGI == 0 and x[0].matches_third_chamber == False])/sum([x[1] for x in counts.items() if x[0].matches_CGI != None])

        line1.append(sum(counts.values()))
        overlap_CGI = sum([x[1] for x in counts.items() if x[0].matches_CGI != None])
        line2.append(overlap_CGI)
        num_snp = sum([x[1] for x in counts.items() if x[0].is_SNP])
        line3.append(num_snp)

        CGI = sum([x[1] for x in counts.items() if x[0].matches_CGI == 0])
        CGI_BAC = sum([x[1] for x in counts.items() if x[0].matches_CGI == 0 and not x[0].matches_BAC == 1])
        CGI_BAC_THIRDCH = sum([x[1] for x in counts.items() if x[0].matches_CGI == 0 and not x[0].matches_BAC == 1 and not x[0].matches_third_chamber])

        line4.append(CGI)
        line5.append(CGI_BAC)
        line6.append(CGI_BAC_THIRDCH)
        line7.append(CGI_BAC_THIRDCH / overlap_CGI)

        SNP_err = sum([x[1] for x in counts.items() if x[0].is_SNP and x[0].matches_CGI == 0 and not x[0].matches_BAC == 1 and not x[0].matches_third_chamber])
        line8.append(SNP_err / num_snp)

    # independent same cell
    counts = pickle.load(open(ind_same_cell_counts,'rb'))

    for cell,cellname in zip([0,1,2],['PGP1_21','PGP1_22','PGP1_A1']):
        line0.append(cellname)
        line1.append(sum([x[1] for x in counts.items() if x[0].cell == cell]))
        overlap_CGI = sum([x[1] for x in counts.items() if x[0].cell == cell and x[0].matches_CGI != None])
        line2.append(overlap_CGI)
        num_snp = sum([x[1] for x in counts.items() if x[0].cell == cell and x[0].is_SNP])
        line3.append(num_snp)

        CGI = sum([x[1] for x in counts.items() if x[0].cell == cell and x[0].matches_CGI == 0])
        CGI_BAC = sum([x[1] for x in counts.items() if x[0].cell == cell and x[0].matches_CGI == 0 and not x[0].matches_BAC == 1])
        CGI_BAC_THIRDCH = sum([x[1] for x in counts.items() if x[0].cell == cell and x[0].matches_CGI == 0 and not x[0].matches_BAC == 1 and not x[0].matches_third_chamber])

        line4.append(CGI)
        line5.append(CGI_BAC)
        line6.append(CGI_BAC_THIRDCH)
        line7.append(CGI_BAC_THIRDCH / overlap_CGI)

        SNP_err = sum([x[1] for x in counts.items() if x[0].is_SNP and x[0].matches_CGI == 0 and not x[0].matches_BAC == 1 and not x[0].matches_third_chamber])
        line8.append(SNP_err / num_snp)

    with open(output_file,'w') as opf:
        for line in [line0,line1,line2,line3,line4,line5,line6,line7,line8]:
            print('\t'.join([str(x) for x in line]),file=opf)

def generate_unphased_calls_table(output_file, count_files, cutoffs):

    countlst = []
    lines = [[] for i in range(8)]
    lines[0].append('')
    lines[1].append('Total Calls')
    lines[2].append('Calls Overlapping Ref')
    lines[3].append('Calls Different from Ref')
    lines[4].append('FPR')
    lines[5].append('Error Rate Upper Bound')
    lines[6].append('Matching SNVs')
    lines[7].append('FDR')

    assert(len(count_files) == len(cutoffs))

    for f,cut in zip(count_files,cutoffs):
        counts = pickle.load(open(f,'rb'))
        countlst.append(counts)
        lines[0].append(cut)
        lines[1].append(sum([x[1] for x in counts.items()]))
        have_truth = sum([x[1] for x in counts.items() if x[0].matches_ref != None])
        lines[2].append(have_truth)
        mismatch = sum([x[1] for x in counts.items() if x[0].matches_ref == False])
        snv_mismatch = sum([x[1] for x in counts.items() if x[0].matches_ref == False and x[0].is_SNP])

        lines[3].append(mismatch)

        true_negatives = sum([x[1] for x in counts.items() if x[0].matches_ref == True and x[0].ref_genotype == '0/0'])

        lines[4].append(snv_mismatch/(snv_mismatch + true_negatives))
        lines[5].append(mismatch/have_truth)
        snv_match =sum([x[1] for x in counts.items() if x[0].matches_ref and x[0].is_SNP])
        lines[6].append(snv_match)
        lines[7].append(mismatch / (snv_match + mismatch))

    with open(output_file,'w') as table:
        for line in lines:
            print('\t'.join([str(x) for x in line]),file=table)

def generate_strand_strand_mismatch_table(output_file, same_cell_counts, all_cell_counts, cross_cell_counts):

    lines = [[] for i in range(4)]
    lines[0].append('')
    lines[1].append('total strand-strand mismatch')
    lines[2].append('total strand-strand match')
    lines[3].append('mismatch rate')

    for mode,f in [('same_cell',same_cell_counts),('all_cell',all_cell_counts),('cross_cell',cross_cell_counts)]:

        counts = pickle.load(open(f,'rb'))

        mismatch = sum([x[1] for x in counts.items() if x[0][1] == 'mismatch'])
        match = sum([x[1] for x in counts.items() if x[0][1] == 'match'])
        mismatch_rate = mismatch / (mismatch + match)

        lines[0].append(mode)
        lines[1].append(mismatch)
        lines[2].append(match)
        lines[3].append(mismatch_rate)

    with open(output_file,'w') as table:
        for line in lines:
            print('\t'.join([str(x) for x in line]),file=table)

def generate_nucleotide_substitution_plot(output_file,input_file):

    counts = pickle.load(open(input_file,'rb'))

    sub_dict = defaultdict(int)
    total_sub = 0
    for k,v in counts.items():
        if k[3] != k[4]:
            sub_dict[(k[3],k[4])] += v
            total_sub += v

    sub_list = sorted(list(sub_dict.items()),key=lambda x: x[1])

    ts = 0
    tv = 0

    x_val = []
    y_val = []

    for k,v in sub_list:
        print("{}\t{}\t{}".format(k[0],k[1],v))
        if {k[0],k[1]} <= {'A','G'} or  {k[0],k[1]} <= {'C','T'}:
            ts += v
        else:
            tv += v

        x_val.append("{} -> {}".format(k[0],k[1]))
        y_val.append(v)

    mpl.rc('legend', fontsize=11)
    mpl.rc('xtick', labelsize=11)
    mpl.rc('ytick', labelsize=11)
    mpl.rc('axes', labelsize=11)
    mpl.rc('axes', labelsize=11)
    mpl.rcParams.update({'font.size': 11})
    mpl.rc('lines', linewidth=1.5)
    mpl.rc('mathtext',default='regular')

    alpha1 = 0.8
    alpha2 = 0.5
    width = 0.75


    ind = np.array(list(range(0,12)))
    ax = plt.subplot(111)
    plt.bar(ind, y_val, color='c',
            ecolor='black',
            alpha=alpha1,
            width=width,
            align='center',
            label='subs')

        # add some text for labels, title and axes ticks
    ax.set_ylabel('Mismatch Counts')
    ax.set_xticks(ind)
    ax.set_xticklabels(tuple(x_val))
    ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=90)

    ax.yaxis.grid(True,color='grey')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")
    #plt.legend(loc='upper left')
    plt.xlim(-1,12)
    plt.savefig(output_file)

    print("G->T rate: {}".format(sub_dict[('G','T')]/total_sub))
    print("G<->T rate: {}".format((sub_dict[('G','T')]+sub_dict[('T','G')])/total_sub))
    print("Ts/Tv: {}".format(ts/tv))
