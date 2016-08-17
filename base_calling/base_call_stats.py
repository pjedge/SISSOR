import numpy as np
import os
from collections import defaultdict
from matplotlib import pyplot as plt
desc = 'Use cross-chamber information to correct base calls in SISSOR data'
sample_rate = 0.0005

chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX','chrY']

cells = ['PGP1_21','PGP1_22','PGP1_A1']
chambers = list(range(1,25))

# next function that returns 0 instead of raising StopIteration
# this is convenient for iterating over file 2 at a time
def safenext(iterator):
    try:
        nextval = next(iterator)
        return nextval
    except StopIteration:
        return 0

short_chrom_names = set([str(x) for x in range(1,23)]+['X','Y'])
chrom_names = set(['chr'+str(x) for x in range(1,23)]+['chrX','chrY'])

def format_chrom(chrom_str):
    if chrom_str in short_chrom_names:
        chrom_str = 'chr' + chrom_str
    assert(chrom_str in chrom_names)
    return chrom_str
        
def base_call_stats():

    # open all 72 files (3 cells x 24 chambers) into a list of handles
    # input files will be read and processed in parallel for equivalent indices
    # output files will contain the same information but will be corrected based on information from other chambers
    
    input_files = []
    output_files = []

    for cell in cells:
        
        outdir = 'pileups_subsample2/{}'.format(cell)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            
        for chamber in chambers:
            
            infile_name = 'pileups_subsample/{}/ch{}.pileup'.format(cell,chamber)
            outfile_name = '{}/ch{}.pileup'.format(outdir,chamber)

            input_files.append(open(infile_name,'r'))
            output_files.append(open(outfile_name,'w'))
        
    # step over files in parallel
    chrom = chroms.pop(0)
    pos = 0
    
    lines = []
    
    for file in input_files:
        
        lines.append(safenext(file))
    
        
    N = len(input_files)

    # print header lines to output files 
    for ix in range(N):
        while 1:
            
            if lines[ix] and lines[ix][0] == '#':
                print(lines[ix],file=output_files[ix],end='')
                lines[ix] = safenext(input_files[ix])
            else:
                break
    
    all_cell_hist = defaultdict(int)
    single_cell_hist = defaultdict(int)
    coverage_cut = 3
    # until all files have reached stopiteration
    while sum([not l for l in lines]) < N:
        
        # "element list"
        # el_lst[x] contains a list with the current line elements for file x (unique cell and chamber)
        el_lst = []
        for l in lines:
            if not l:
                el_lst.append(0)
            else:
                el_lst.append(l.strip().split('\t'))

        # fix inconsistent chromosome names
        for i in range(len(el_lst)):
            if el_lst[i]:
                el_lst[i][0] = format_chrom(el_lst[i][0])
            
        on_chrom = [el and el[0] == chrom for el in el_lst]
        on_chrom_pos = [el and el[0] == chrom and int(el[1])-1 == pos for el in el_lst]

        # if there are no more files on our chromosome, we move onto the next chromosome
        if sum(on_chrom) == 0:
            if chroms == []:
                break
            else:
                chrom = chroms.pop(0)
                pos = 0
                continue
        
        # now we know that some subset of the N files are processing our chromosome
        # if none of the files are on our genomic index we also need to iterate forward
        # until some are
        
        if sum(on_chrom_pos) == 0:
            pos = float('Inf')  # positions for each file
            for el in el_lst:
                if el and el[0] == chrom and int(el[1])-1 < pos:
                    pos = int(el[1])-1  # select minimum position index of all currently considered lines
        
        on_chrom_pos = [el and el[0] == chrom and int(el[1])-1 == pos for el in el_lst]
        assert(sum(on_chrom_pos) > 0)
        
        # now we have a valid chromosome, being considered by some file handles
        # we also have a valid position index, of which some file handles are on it.
        
        # process each infile with a 1 in on_chrom_pos
        # then iterate each of those files
        
        
        # statistics to run: 
        # - histogram of total number of chambers with fragments above a given coverage (start with 5) for a given position
        #      - we expect at most 12
        # - histogram of total number of cells with fragments above a given coverage (start with 5) for a given position
        #      - we expect at most 4
        
        # all cells strand count histogram
        counts = 0
        for i in np.where(on_chrom_pos)[0]:
            el = el_lst[i]
            if not el:
                continue
            
            if int(el[3]) >= coverage_cut:
                counts += 1
            
        all_cell_hist[counts] += 1
        
        # single cells strand count histogram
        counts1 = 0
        counts2 = 0
        counts3 = 0
        for i in np.where(on_chrom_pos)[0]:
            el = el_lst[i]
            if not el:
                continue
            
            if int(el[3]) >= coverage_cut:
                if i < 24:
                    counts1 += 1
                elif i < 48:
                    counts2 += 1
                else:
                    counts3 += 1
            
        single_cell_hist[counts1] += 1
        single_cell_hist[counts2] += 1
        single_cell_hist[counts3] += 1
        
        for i in np.where(on_chrom_pos)[0]:
            
            #el = el_lst[i]
            
            #new_line = '\t'.join(el)
            #print(new_line,file=output_files[i])
            
            #print(lines[i],file=output_files[i],end='')
            lines[i] = safenext(input_files[i])
            
    print("Combined cells position coverage")
    xdata = []
    ydata = []
    for k,v in sorted(list(all_cell_hist.items())):
        print('{}\t{}'.format(k,v))
        xdata.append(k)
        ydata.append(v)

    plt.figure()
    plt.plot(xdata[1:],ydata[1:])
    plt.title("Combined cells position coverage")
    plt.xlabel("Coverage")
    plt.ylabel("Number of positions")
    plt.xlim((1,10))
    print("Separate cells position coverage")
    xdata = []
    ydata = []
    for k,v in sorted(list(single_cell_hist.items())):
        print('{}\t{}'.format(k,v))
        xdata.append(k)
        ydata.append(v)

    plt.figure()
    plt.plot(xdata[1:],ydata[1:])
    plt.title("Separate cells position coverage")
    plt.xlabel("Coverage")
    plt.ylabel("Number of positions")
    plt.xlim((1,5))
    plt.xticks([1,2,3,4,5])

    
    
    
    

    # close all chambers
            
    for handle in input_files:
        handle.close()
    for handle in output_files:
        handle.close()            

if __name__ == '__main__':
    base_call_stats()
