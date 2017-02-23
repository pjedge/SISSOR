from collections import defaultdict, namedtuple
import pysam
import pickle
import re

short_chrom_names = set([str(x) for x in range(1,23)]+['X','Y'])
chrom_names = set(['chr'+str(x) for x in range(1,23)]+['chrX','chrY'])
chroms = ['chr{}'.format(i) for i in range(1,23)] + ['chrX','chrY']
count_key = namedtuple('count_key', 'is_SNP matches_CGI half_match mismatch')
cells = ['PGP1_21','PGP1_22','PGP1_A1']
chambers = list(range(1,25))
bases = {'A','C','G','T'}
n_chambers = 24
n_cells = 3

dp_pat = re.compile("DP=(\d+)")
qr_pat = re.compile(";QR=(\d+)")
type_pat = re.compile("TYPE=([^;]+);")

short_chrom_names = set([str(x) for x in range(1,23)]+['X','Y'])
chrom_names = set(['chr'+str(x) for x in range(1,23)]+['chrX','chrY'])

def format_chrom(chrom_str):
    if chrom_str in short_chrom_names:
        chrom_str = 'chr' + chrom_str
    if chrom_str not in chrom_names:
        return False
    return chrom_str

chr_num = dict()
for i,chrom in enumerate(chroms):
    chr_num[chrom] = i

def split_vcf(input_vcf, chunklist, output_vcfs):

    regions_output = list(zip(chunklist, output_vcfs))
    (chrom, start, stop), outputfile = regions_output.pop(0)
    output = open(outputfile, 'w')
    with open(input_vcf,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')

            vcf_chrom = format_chrom(vcf_line[0])
            if not vcf_chrom:
                continue
            vcf_pos = int(vcf_line[1])

            if vcf_chrom != chrom or vcf_pos > stop:

                done = False
                while not (vcf_chrom == chrom and vcf_pos >= start and vcf_pos <= stop):
                    output.close()
                    if len(regions_output) == 0:
                        done = True
                        break
                    (chrom, start, stop), outputfile = regions_output.pop(0)
                    output = open(outputfile, 'w')
                if done:
                    break

            # write to lifted over vcf
            print(line,end='',file=output)

    while len(regions_output) > 0:
        output.close()
        (chrom, start, stop), outputfile = regions_output.pop(0)
        output = open(outputfile, 'w')

    if not output.closed:
        output.close()

# next function that returns 0 instead of raising StopIteration
# this is convenient for iterating over file 2 at a time
def safenext(iterator):
    try:
        nextval = next(iterator)
        return nextval
    except StopIteration:
        return 0


class gff_reader():
    def __init__(self, GFF_file, fasta_file):

        self.gff_iter = iter(self.read_gff(GFF_file, fasta_file))
        self.chrom = None
        self.pos = None
        self.allele = None

    def get(self, chrom, pos):

        # we're already on the correct value
        if (chrom, pos) == (self.chrom, self.pos):
            return self.allele

        # we're already past the correct value
        if (self.chrom, self.pos) != (None, None) and (
                chr_num[self.chrom] > chr_num[chrom] or (self.chrom == chrom and self.pos > pos)):
            return False

        nxt = safenext(self.gff_iter)
        if not nxt:
            return False
            #raise Exception('GFF Iterator spun out of bounds. Did you try to access a previous genomic index?')

        self.chrom, self.pos, self.allele = nxt
        if (chr_num[self.chrom] > chr_num[chrom] or (self.chrom == chrom and self.pos > pos)):
            return False
        while ((self.chrom, self.pos) != (chrom, pos)):
            nxt = safenext(self.gff_iter)
            if not nxt:
                return False
                #raise Exception('GFF Iterator spun out of bounds. Did you try to access a previous genomic index?')

            self.chrom, self.pos, self.allele = nxt

            if (chr_num[self.chrom] > chr_num[chrom] or (self.chrom == chrom and self.pos > pos)):
                return False

        assert ((self.chrom, self.pos) == (chrom, pos))
        return self.allele

    '''
    def read_vcf(self, ):
        # add heterozygous SNVs observed in our second CGI dataset to the dictionary
        with open(CGI_VCF_file,'r') as vcf:
            for line in vcf:
                if line[0] == '#' or len(line) < 3:
                    continue

                vcf_line = line.strip().split('\t')
                vcf_chrom = vcf_line[0]
                vcf_pos   = int(vcf_line[1])

                if not (vcf_chrom == region_chrom and vcf_pos >= region_start and vcf_pos <= region_end):
                    continue

                alleles = {vcf_line[3],vcf_line[4]}

                vcf_dict[(vcf_chrom,vcf_pos)] = alleles

        return vcf_dict
    '''

    def read_gff(self, GFF_file, fasta_file):
        hg19_fasta = pysam.FastaFile(fasta_file)
        # add heterozygous alleles observed in one CGI dataset to a dictionary
        with open(GFF_file, 'r') as gff:
            for line in gff:
                if line[0] == '#' or len(line) < 3:
                    continue
                gff_line = line.strip().split('\t')

                if gff_line[2] == 'SNP':

                    gff_chrom = gff_line[0]
                    gff_pos = int(gff_line[3])

                    assert (gff_pos == int(gff_line[4]))

                    infoline = gff_line[8]
                    info = infoline.strip().split(';')
                    a_info = info[0]
                    assert ('alleles' in a_info)
                    a_info = a_info[8:]

                    # ref_allele = None
                    # found = re.findall(ref_allele_re, infoline)
                    # if len(found) > 0:
                    #    ref_allele = found[0]

                    if len(a_info) == 1:
                        alleles = {a_info}
                    elif len(a_info) == 3:
                        alleles = {a_info[0], a_info[2]}
                    else:
                        print("Unexpected number of alleles")
                        assert (False)

                    yield (gff_chrom, gff_pos, alleles)
                    # we observed only 1 SNP allele and hasn't been seen in other datasets
                    # add reference base because of this CGI dataset's base false negative rate
                    # if len(alleles) == 1 and ref_allele != None:
                    #    truth_dict[(gff_chrom,gff_pos)] = truth_dict[(gff_chrom,gff_pos)].union({ref_allele})

                elif gff_line[2] == 'REF':

                    gff_chrom = gff_line[0]
                    gff_start = int(gff_line[3])
                    gff_end = int(gff_line[4])

                    for gff_pos in range(gff_start, gff_end + 1):

                        ref_lookup = str.upper(hg19_fasta.fetch(gff_chrom, gff_pos - 1, gff_pos))  # pysam uses 0-index

                        if ref_lookup not in bases:  # e.g. 'N' alleles
                            continue

                        yield (gff_chrom, gff_pos, {ref_lookup})

def format_chrom(chrom_str):
    if chrom_str in short_chrom_names:
        chrom_str = 'chr' + chrom_str
    if chrom_str not in chrom_names:
        return False
    return chrom_str

def where(lst):
    return [i for i in range(len(lst)) if lst[i]]

def freebayes_strand_pair(input_files, fragment_assignment_file, gff_file, fasta_file, count_file, mismatch_file, same_cell_only):

    input_files = [open(inf,'r') for inf in input_files]
    gffr = gff_reader(gff_file,fasta_file)
    print('Reading fragment pair file...')
    counts = defaultdict(int)
    paired_fragments = [] # list of (start, end, ix1, ix2) tuples
    with open(fragment_assignment_file,'r') as infile:
        for line in infile:
            el = line.strip().split()
            chrom = el[0]
            start = int(el[1])
            end   = int(el[2])
            cell1 = int(el[3])
            chamber1 = int(el[4])
            cell2 = int(el[5])
            chamber2 = int(el[6])

            paired_fragments.append((chrom, start, end, cell1, chamber1, cell2, chamber2))

    paired_fragments.sort(key=lambda x: (chr_num[x[0]],x[1]))

    current_paired_fragments = []

    # open all 72 files (3 cells x 24 chambers) into a list of handles
    # input_files should be ordered as cell1.ch1, cell1,ch2, cell1.ch3..., cell2.ch1...

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
                lines[ix] = safenext(input_files[ix])
            else:
                break

    processed = 0
    print('Iterating through VCF file...')

    # until all files have reached stopiteration
    with open(mismatch_file,'w') as mof:
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
            on_chrom_pos = [el and el[0] == chrom and int(el[1]) == pos for el in el_lst]

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
                    if el and el[0] == chrom and int(el[1]) < pos:
                        pos = int(el[1])  # select minimum position index of all currently considered lines

            on_chrom_pos = [el and el[0] == chrom and int(el[1]) == pos for el in el_lst]
            assert(sum(on_chrom_pos) > 0)

            # now we have a valid chromosome, being considered by some file handles
            # we also have a valid position index, of which some file handles are on it.

            # process each infile with a 1 in on_chrom_pos
            # then iterate each of those files

            # do something with line element that is on the current genomic index

            cgi_allele = gffr.get(chrom,pos)
            # pop info for haplotype-paired regions into the "current" list
            while paired_fragments != [] and (chr_num[paired_fragments[0][0]] < chr_num[chrom] or (paired_fragments[0][0] == chrom and paired_fragments[0][1] <= pos)):
                current_paired_fragments.append(paired_fragments.pop(0))

            # filter out paired fragments that we're past
            criteria = lambda pf: (pf[0] == chrom and pf[1] <= pos and pos <= pf[2])
            current_paired_fragments = list(filter(criteria,current_paired_fragments))

            ref_allele = None

            allele_dict = dict()

            for i in where(on_chrom_pos):

                cell_num = int(i / n_chambers)
                ch_num   = i % n_chambers

                el = el_lst[i]
                if not el:
                    continue

                if ref_allele == None:
                    ref_allele = {str.upper(el[3])}
                else:
                    assert(ref_allele == {str.upper(el[3])})

                variant_alleles = set(el[4].split(','))

                if ref_allele == {'N'}:
                    continue

                if len(el) < 10:
                    continue

                vcf_chrom = format_chrom(el[0])
                vcf_pos   = int(el[1])

                assert(vcf_chrom == chrom)
                assert(vcf_pos   == pos)

                fields1 = el[7]
                labels = el[8]
                assert(labels[:2] == 'GT')
                fields2 = el[9]
                genotype = fields2[0]

                # filter out depth < 5
                depth = int(float(re.findall(dp_pat, fields1)[0]))
                if depth < 5:
                    continue

                # filter out qual < 30
                qual = int(float(re.findall(qr_pat, fields1)[0])) if genotype == '0' else float(el[5])
                if qual < 30:
                    continue

                # we only care about reference and SNPs, no indels or MNPs
                pos_type = re.findall(type_pat, line)
                if not (pos_type == [] or pos_type[0] == 'snp'):
                    continue
                if not (len(el[3]) == 1 and (len(el[4]) == 1 or len(el[4]) == 3)):
                    continue

                allele = None
                if genotype == '0':
                    allele = ref_allele
                else:
                    allele = variant_alleles
                # skip cases with N alleles, or with any other unusual cases
                if len(allele.intersection(bases)) != len(allele):
                    continue

                if len(allele) > 2:
                    continue

                allele_dict[(cell_num, ch_num)] = allele

            allele_list = list(allele_dict.items())

            is_SNP      = False
            matches_CGI = True if cgi_allele else None
            half_match  = False
            mismatch    = False
            num_pairs = 0
            for (chrom, start, end, cell1, chamber1, cell2, chamber2) in current_paired_fragments:

                if cell1 != cell2 and same_cell_only:
                    continue

                if (cell1,chamber1) in allele_dict and (cell2,chamber2) in allele_dict:

                    num_pairs += 1
                    allele1 = allele_dict[(cell1,chamber1)]
                    allele2 = allele_dict[(cell2,chamber2)]

                    if allele1 != ref_allele or allele2 != ref_allele:
                        is_SNP = True

                    if cgi_allele and not (allele1 <= cgi_allele and allele2 <= cgi_allele):
                        matches_CGI = False
                    if allele1 != allele2:
                        if allele1 <= allele2 or allele1 >= allele2:
                            half_match = True
                        else:
                            mismatch = True

            if num_pairs > 0:

                counts[count_key(is_SNP,matches_CGI,half_match,mismatch)] += 1

                if (half_match or mismatch) and matches_CGI == False:

                    print('Found a mismatch at {} {}...'.format(chrom, pos))

                    if half_match and mismatch:
                        print('HALF_MATCH and MISMATCH!',file=mof)
                    elif half_match:
                        print('HALF_MATCH!',file=mof)
                    elif mismatch:
                        print('MISMATCH!',file=mof)
                    print('CGI ALLELE: {}'.format(''.join(list(cgi_allele))),file=mof)
                    print('STRAND_PAIRS:',file=mof)
                    for (chrom, start, end, cell1, chamber1, cell2, chamber2) in current_paired_fragments:
                        print("{}:{} paired to {}:{}".format(cell1, chamber1, cell2, chamber2),file=mof)
                    print('VCF INFO:',file=mof)

                    for i in where(on_chrom_pos):

                        cell_num = int(i / n_chambers)
                        ch_num   = i % n_chambers
                        vcf_info = ['{}:{}'.format(cell_num,ch_num)] + el_lst[i]
                        print('\t'.join(vcf_info),file=mof)

                    print('**************************************************',file=mof)

                processed += 1
#            if processed > 500:
#                return
            # iterate to the next lines for each line on the current index
            for i in where(on_chrom_pos):

                lines[i] = safenext(input_files[i])

    print('Dumping counts to pickle file {}...'.format(count_file))
    pickle.dump(counts, open(count_file,'wb'))
    # close all chambers

    for handle in input_files:
        handle.close()

    print('Program complete.')


def accuracy_aggregate(counts_pickle_files, pickle_outfile):

    counts = defaultdict(int)

    for pfile in counts_pickle_files:
        temp_dict = pickle.load(open(pfile,'rb'))
        for k,v in temp_dict.items():
            counts[k] += v

    pickle.dump(counts,open(pickle_outfile,'wb'))
    #print("Strand Mismatch, not in CGI, divided by total:")
    #print(sum([x[1] for x in counts.items() if x[0][1] != None and x[0][2] == True]) / sum([x[1] for x in counts.items() if x[0][1] != None]))
