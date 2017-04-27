import pickle

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
        CGI_WGS_BAC = sum([x[1] for x in counts.items() if x[0].matches_CGI == 0 and not x[0].matches_BAC == 1])
        CGI_WGS_BAC_THIRDCH = sum([x[1] for x in counts.items() if x[0].matches_CGI == 0 and not x[0].matches_BAC == 1 and not x[0].matches_third_chamber])

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

    for f in count_files:
        counts = pickle.load(open(f,'rb'))
        countlst.append(counts)
        lines[0].append(cut)
        lines[1].append(sum([x[1] for x in counts.items()]))
        have_truth = sum([x[1] for x in counts.items() if x[0].matches_CGI != None])
        lines[2].append(have_truth)
        mismatch = sum([x[1] for x in counts.items() if x[0].matches_CGI == False])
        snv_mismatch = sum([x[1] for x in counts.items() if x[0].matches_CGI == False and x[0].is_SNP])

        lines[3].append(mismatch)

        true_negatives = sum([x[1] for x in counts.items() if x[0].matches_CGI == True and x[0].CGI_genotype == '0/0'])

        lines[4].append(snv_mismatch/(snv_mismatch + true_negatives))
        lines[5].append(mismatch/have_truth)
        snv_match =sum([x[1] for x in counts.items() if x[0].matches_CGI and x[0].is_SNP])
        lines[6].append(snv_match)
        lines[7].append(mismatch / (snv_match + mismatch))

    with open(output_file,'w') as table:
        for line in lines:
            print('\t'.join([str(x) for x in line]),file=table)
