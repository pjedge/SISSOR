from math import log10
from collections import defaultdict
import os

cells = ['PGP1_21','PGP1_22','PGP1_A1']
chambers = list(range(1,25))
in1 = 'base_calling/het_vcfs/cutoff3/whole_genome.out'
in2 = []
for cell in cells:
    for chamber in chambers:
        in2.append('base_calling/fragment_boundary_beds/{}/ch{}.bed'.format(cell,chamber))
in3 = 'output_dir_CCF_frag'
        
n_cells = 3
n_chambers = 24
XCHAMBER = True
MIN_Q    = 30
MIN_COV  = 5

chroms = ['chr{}'.format(x) for x in range(1,23)] + ['chrX','chrY']
chrom_num = dict()
for i,chrom in enumerate(chroms):
    chrom_num[chrom] = i
    
def parse_bedfile(input_file):

    boundaries = []
    with open(input_file,'r') as inf:
        for line in inf:

            if len(line) < 3:
                continue

            el = line.strip().split('\t')

            chrom = el[0]
            start = int(el[1])
            end = int(el[2])

            boundaries.append((chrom, start, end))

    return boundaries    
    
def create_hapcut_fragment_matrix_CCF(chamber_call_file=in1, fragment_boundary_files=in2, output_dir=in3):
    
    fragment_boundaries = []
    for i in range(0,n_cells*n_chambers):
        bfile = fragment_boundary_files[i]
        fragment_boundaries.append(parse_bedfile(bfile))
        
    #fragment_list[i][j] is the jth fragment
    fragment_list = [[[]] for i in range(n_chambers*n_cells)]
        
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
                
                if fragment_boundaries[i] == []: #no more fragments
                    continue
                
                frg_chrom, frg_start, frg_end = fragment_boundaries[i][0]
                
                # fragment boundaries are behind our positions
                while (chrom_num[frg_chrom] < chrom_num[ccf_chrom]) or (frg_chrom == ccf_chrom and frg_end < ccf_pos):

                    fragment_boundaries[i].pop(0)
                    if fragment_boundaries[i] == []:
                        frg_chrom = None
                        break
                    frg_chrom, frg_start, frg_end = fragment_boundaries[i][0]
                    
                    if fragment_list[i][-1] != []: # if the last fragment isn't empty, start fresh with a new one
                        fragment_list[i].append([])
                    
                # fragment boundary is now either ahead of us or we're inside it
                # if we're not inside boundary then skip ahead
                if not (frg_chrom == ccf_chrom and frg_start <= ccf_pos and frg_end >= ccf_pos):
                    continue
                
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
                    
    lines = defaultdict(list)
    fragcount = 0
        # now take this information and make it into a fragment matrix file line
    for i, frag_snps in enumerate(fragment_list):

        for fs in frag_snps:
            if len(fs) < 2:
                continue
            fragcount += 1

            fragstr = ''
            num_pairs = 0
            prev_snp_ix = -2
            qual = ''
            chrom     = fs[0][1]
            firstpos  = fs[0][2]
            frg_start = fs[0][5]
            frg_end   = fs[0][6]
            cell      = fs[0][7]
            ch_num    = fs[0][8]
            
            name = '{}:{}-{}:{}:CH{}'.format(chrom,frg_start,frg_end,cell,ch_num)

            for snp_ix, chrom, pos, allele, q_char, frg_start, frg_end, cell, ch_num in fs:

                diff = snp_ix - prev_snp_ix

                if diff == 1:
                    fragstr += allele
                else:
                    num_pairs += 1
                    fragstr += ' {} {}'.format(snp_ix+1, allele)

                prev_snp_ix = snp_ix
                qual += q_char

            fragstr += ' ' + ''.join(qual)

            prefix = '{} {}'.format(num_pairs,name)
            fragstr = prefix + fragstr

            lines[chrom].append((firstpos, fragstr))

    # go through the output lines we've accrued for each chromosome.
    # sort them and print them to a fragment file for each chromosome.
    for chrom in lines.keys():
        lines[chrom].sort()

        output_fragment_file = os.path.join(output_dir, chrom)

        with open(output_fragment_file, 'w') as opf:
            for firstpos, line in lines[chrom]:
                print(line, file=opf)

if __name__ == '__main__':
    create_hapcut_fragment_matrix_CCF()