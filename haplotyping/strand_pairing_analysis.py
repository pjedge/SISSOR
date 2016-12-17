# -*- coding: utf-8 -*-
"""
Created on Thu May 19 21:44:48 2016

@author: peter
"""
import sys
sys.path.append('../../HapTools')
#import fileIO
import fragment
import sys
from math import log10
from itertools import combinations

# OBJECTIVE:
# pair haplotype strands and determine three classes
# SPSC = SAME PARENT SAME CELL
# DPSC = DIFFERENT PARENT SAME CELL
# SPDC = SAME PARENT DIFFERENT CELL
# DPDC = DIFFERENT PARENT DIFFERENT CELL

# we want to read in a haplotype fragment file in HAPCUT2 format, find overlaps,
# and filter based on number of het SNVs overlapped and the consistency
# perhaps find a set of SNVs of largest span that has no switch errors?

# output a bed file with chromosoe, start, stop of region (based on het SNV span)
# and the cells and chambers that contain the paired fragments
# might as well also give the fragment IDs too.

# This bed file will allow us to parse a CCF file and mark pairs of chambers
# based on the verification they've undergone, i.e.
# SPSC.1.2.15 = same parent same cell, match, to cell 2, chamber 15
# SPSC.0.2.15 = same parent same cell, MISMATCH, to cell 2, chamber 15


VERBOSE = False
DEBUG = False
ID_LEN_NO_MODS = 5

tinylog = -1e5


def total_SER(f1, f2):

    assert(f1.seq[0][0] <= f2.seq[0][0])

    f2iter = iter(f2.seq) # iterator for fragment2
    a2 = safenext(f2iter)
    if not a2:
        return 0

    switches = 0
    total_pos = 0
    done  = False
    first_pos = True
    switched = False

    for a1 in f1.seq:

        if done:
            break

        while a1[0] > a2[0]:
            a2 = safenext(f2iter)
            if not a2:
                done = True
                break

        if done or a1[0] < a2[0]:
            continue

        assert a1[0] == a2[0]
        if first_pos:
            switched = (a1[2] != a2[2])
            first_pos = False

        total_pos += 1
        if (a1[2] != a2[2] and not switched) or (a1[2] == a2[2] and switched):
            switches += 1
            switched = not switched

        a2 = safenext(f2iter)
        if not a2:
            break

    return switches, total_pos

# next function that returns 0 instead of raising StopIteration
# this is convenient for iterating over file 2 at a time
def safenext(iterator):
    try:
        nextval = next(iterator)
        return nextval
    except StopIteration:
        return 0

# determine if two fragments are consistent with each other
def consistent(f1, f2, i,j, contam_bounds, t):

    answer = True
    assert(f1.seq[0][0] <= f2.seq[0][0])

    f2iter = iter(f2.seq) # iterator for fragment2
    a2 = safenext(f2iter)
    if not a2:
        return answer, contam_bounds

    mismatches = 0
    done  = False
    first_pos = True
    switched = False
    possible_switch_pos = None
    last_pos = None

    for a1 in f1.seq:

        if done:
            break

        while a1[0] > a2[0]:
            a2 = safenext(f2iter)
            if not a2:
                done = True
                break

        if done or a1[0] < a2[0]:
            continue

        assert a1[0] == a2[0]
        if first_pos:
            switched = (a1[2] != a2[2])
            first_pos = False

        if (a1[2] == a2[2] and not switched) or (a1[2] != a2[2] and switched):
            mismatches = 0
            possible_switch_pos = None
        else:
            if possible_switch_pos == None:
                possible_switch_pos = set(range(last_pos+1,a1[0]+1))
            mismatches += 1
            if mismatches >= t:
                contam_bounds[(i,j)] = set.union(contam_bounds[(i,j)], possible_switch_pos)
                answer = False
                switched = not switched
                mismatches = 0
                possible_switch_pos = None

        last_pos = a1[0]

        a2 = safenext(f2iter)
        if not a2:
            break

    contam_bounds[(j,i)] = contam_bounds[(i,j)]

    return answer, contam_bounds


def filter_discordant_fragments_SER(flist, SER_threshold):

    flist = sorted(flist,key=lambda x: x.seq[0][0])

    N = len(flist)


    new_flist = []

    for i in range(N):
        f1 = flist[i]

        errs = []
        e_total = 0
        t_total = 0
        for j in range(i+1,N):

            f2 = flist[j]

            if not overlap(f1,f2,2):

                continue

            e, t = total_SER(f1,f2)
            e_total += e
            t_total += t

            errs.append(e/t)

        if t_total <= 10 or e_total / t_total < SER_threshold: #len(errs) <= 1 or min(errs) < SER_threshold:
            new_flist.append(f1)

    return new_flist

def addlogs(a,b):
    if a > b:
        return a + log10(1 + pow(10, b - a))
    else:
        return b + log10(1 + pow(10, a - b))

CUTOFF = 0.99

def overlap(f1, f2, amt=3):

    if f1.seq[0][0] > f2.seq[0][0]:
        
        temp = f1
        f1 = f2
        f2 = temp

    assert(f1.seq[0][0] <= f2.seq[0][0])

    s1 = set([a[0] for a in f1.seq])
    s2 = set([a[0] for a in f2.seq])
    
    inter = set.intersection(s1,s2)
    
    if len(inter) < amt:
        return None, None
    
    inter_seq1 = [x for x in f1.seq if x[0] in inter]  
    inter_seq2 = [x for x in f2.seq if x[0] in inter] 

    gpos = [x[1] for x in inter_seq1]  
    start = min(gpos)
    end   = max(gpos)
    
    p_samehap = 0
    p_diffhap = 0
    total = 0
    matches = 0
    for (snp_ix1, gen_ix1, a1, q1), (snp_ix2, gen_ix2, a2, q2) in zip(inter_seq1, inter_seq2):
        
        
        q1 = 10**((ord(q1)-33)/-10)        
        q2 = 10**((ord(q2)-33)/-10)   

        total += 1
        if a1 == a2:
            matches += 1
            p_samehap += log10((1-q1)*(1-q2)+q1*q2)
            p_diffhap += log10((1-q1)*q2+(1-q2)*q1)
        else:
            p_diffhap += log10((1-q1)*(1-q2)+q1*q2)
            p_samehap += log10((1-q1)*q2+(1-q2)*q1)
    
    #post_samehap = 10**(p_samehap - addlogs(p_samehap, p_diffhap))
    #print("{} {} {}".format(matches, total, post_samehap))
    
    if matches/total > 0.8:#post_samehap > CUTOFF:
        
        return start, end
        
    else:
        return None, None


def assign_fragments(flist, outputfile):#hapblocks):

    total = 0
    
    with open(outputfile, 'w') as opf:
        for f1,f2 in combinations(flist, 2):
                
            start, end = overlap(f1,f2,4)
            
            if start == None:
                continue
            
            total += end - start
            
            el1 = f1.name.split(':')
            chrom = el1[0]
            cell1 = el1[3]
            chamber1 = el1[4]
            el2 = f2.name.split(':')
            cell2 = el2[3]
            chamber2 = el2[4]
    
            print('{}\t{}\t{}\t{}\t{}'.format(chrom, start, end, cell1, chamber1, cell2, chamber2), file=opf)
            
    print(total)
    
    
cell_name_map = {'PGP1_21':0,'PGP1_22':1,'PGP1_A1':2}    

def compared_paired_strands(chamber_call_file=in1, fragment_assignment_file, output_file=in1):
    
    fragment_boundaries = []
    for i in range(0,n_cells*n_chambers):
        bfile = fragment_boundary_files[i]
        fragment_boundaries.append(parse_bedfile(bfile))
        
    #fragment_list[i][j] is the jth fragment
    fragment_list = [[[]] for i in range(n_chambers*n_cells)]
    
    
    # read in file with spans of paired strands with matching haplotypes
    paired_fragments = [] # list of (start, end, ix1, ix2) tuples
    with open(fragment_assignment_file,'r') as infile:
        for line in infile:
            el = line.strip().split()
            chrom = el[0]
            start = el[1]
            end   = el[2]
            cell1 = cell_name_map[el[3]]
            chamber1 = int(el[4])
            ix1 = cell1 * n_chambers + chamber1
            cell2 = cell_name_map[el[5]]
            chamber2 = int(el[6])
            ix2 = cell2 * n_chambers + chamber2
            
            paired_fragments.append((chrom, start, end, ix1, ix2))

    with open(chamber_call_file,'r') as ccf:
        #print('chr\tpos\tsissor_call\tCGI_allele\tref',file=mof)
        
        snp_ix = -1
        prev_ccf_chrom = None
        for line in ccf:
            ccf_line = line.strip().split('\t')
            ccf_chrom = ccf_line[0]

            if ccf_chrom != prev_ccf_chrom:
                snp_ix = -1
            prev_ccf_chrom = ccf_chrom

            snp_ix += 1

            ccf_pos   = int(ccf_line[1])
            ref_allele = {ccf_line[2]}

            if ccf_chrom in ['chrX','chrY']:
                continue
            
            tags = ccf_line[80].split(';')
            if 'TOO_MANY_ALLELES' in tags or 'TOO_MANY_CHAMBERS' in tags or 'ADJACENT_INDEL_OR_CLIP' in tags:
                continue            
            
            call = ccf_line[3]

            el2 = call.split(';')

            #genotype_prob = -1
            #max_genotype = 'NN'
            #for entry in el2:

            #    genotype, prob = entry.split(':')
            #    prob = float(prob)

            #    if prob > genotype_prob:
            #        genotype_prob = prob
            #        #max_genotype = genotype

            base_call_list = [x for x in ccf_line[5:80] if 'CELL' not in x]
            
            for i,call in enumerate(base_call_list):                
                
                xchamber_calls, basic_calls, pileup = call.split('|')                
                                                
                cell_num = int(i / n_chambers)
                ch_num   = int(i % n_chambers)+1
                cell = cells[cell_num]
                
                
                
                if len(pileup) < MIN_COV:
                    continue
                            
                if XCHAMBER:
                    call = xchamber_calls
                else:
                    call = basic_calls
                
                if call == '*':
                    continue

                el2 = call.split(';')

                max_prob = -1
                max_allele = 'N'
                for entry in el2:

                    a_info, prob = entry.split(':')
                    prob = float(prob)

                    if len(a_info) == 1:
                        allele = {a_info}
                    elif len(a_info) == 2:
                        allele = {a_info[0],a_info[1]}

                    if prob > max_prob:
                        max_prob = prob
                        max_allele = allele
                
                if max_prob < MIN_Q:
                    continue
                
                if len(max_allele) == 2:
                    binary_allele = 'M'
                elif max_allele == ref_allele:
                    binary_allele = '0'
                else:
                    binary_allele = '1'
                    
                p_err = 1 - max_prob
                if p_err < 1e-10:
                    p_err = 1e-10
                qual = int(-10 * log10(p_err))

                q_char = '~' if qual>=93 else chr(33 + qual)                
                
                fragment_list[i][-1].append((snp_ix, ccf_chrom, ccf_pos, binary_allele, q_char, frg_start, frg_end, cell, ch_num))
                    
    
def pair_strands(fragmentfile='../fragmat_analysis/augmented_fragmats/chr1.all', vcf_file='../fragmat_analysis/PGP1_VCFs/chr1.vcf'):
    
    # READ HAPLOTYPE BLOCKS

    # READ FRAGMENT MATRIX
    flist = fragment.read_fragment_matrix(fragmentfile,vcf_file)
    
    # ASSIGN FRAGMENTS TO HAPLOTYPES
    flist = assign_fragments(flist, 'foobar.txt')
    
if __name__ == '__main__':
    pair_strands()
