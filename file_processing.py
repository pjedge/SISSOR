import re
import sys

caret_money = re.compile(r'\^.|\$')
comma_or_period = re.compile(r'\.|\,')
plus_or_minus = re.compile(r'\+|-')
coverage_cut = 2
base_letters = 'ACGTacgt'
numbers = '0123456789'
bases = {'A','T','G','C'}

def parse_mpileup_base_qual(raw_bd, raw_qd, ref_base):

    bd = [] # base data
    qd = [] # qual data

    indel_count = len(re.findall(caret_money, raw_bd))

    if not re.search(plus_or_minus,raw_bd):

        bd = str.upper(re.sub(caret_money,'', raw_bd))
        bd = re.sub(comma_or_period, ref_base, bd)
        qd = [10**((ord(q) - 33) * -0.1) for q in raw_qd]
        paired_bd_qd = [(b,q) for b,q in zip(bd,qd) if b not in ['>','<','*','n','N']]
        if len(paired_bd_qd) < coverage_cut:
            return [],[], 0

        bd, qd = zip(*paired_bd_qd)
        bd = list(bd)
        qd = list(qd)

    else:

        i = 0   # index for base data
        j = 0   # index for qual data

        while i < len(raw_bd):
            if raw_bd[i] == '.' or raw_bd[i] == ',':   # reference base call
                bd.append(ref_base)
                qd.append(10**((ord(raw_qd[j]) - 33) * -0.1))
                i += 1
                j += 1
            elif raw_bd[i] in base_letters:            # variant call
                bd.append(str.upper(raw_bd[i]))
                qd.append(10**((ord(raw_qd[j]) - 33) * -0.1))
                i += 1
                j += 1
            elif raw_bd[i] == '+' or raw_bd[i] == '-': # indel
                indel_count += 1
                num = int(raw_bd[i+1])
                i += 2
                while(raw_bd[i] in numbers):
                    num *= 10
                    num += int(raw_bd[i])
                    i += 1
                i += num
            elif raw_bd[i] == '^':
                i += 2
            elif raw_bd[i] == '$':
                i += 1
            elif raw_bd[i] in '><*nN':                 # reference skip or deletion
                i += 1
                j += 1

        assert(i == len(raw_bd))
        assert(j == len(raw_qd))

    assert(len(bd) == len(qd))

    for b in bd:
        assert(b in bases)

    return bd, qd, indel_count

def parse_bedfile(input_file):

    boundaries = []
    with open(input_file,'r') as inf:
        for line in inf:

            if len(line) < 3:
                continue

            el = line.strip().split('\t')

            chrom = el[0]
            start = int(el[1])
            stop  = int(el[2])

            boundaries.append((chrom, start, stop))

    return boundaries

def split_vcf(input_vcf, chunklist, output_vcfs):

    regions_output = list(zip(chunklist, output_vcfs))
    (chrom, start, stop), outputfile = regions_output.pop(0)
    output = open(outputfile, 'w')
    with open(input_vcf,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')

            vcf_chrom = vcf_line[0]
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

    if not output.closed:
        output.close()

def prune_hapblock_file(hapblock_file, output_file, snp_conf_cutoff, split_conf_cutoff, use_refhap_heuristic):

    with open(hapblock_file,'r') as inf, open(output_file,'w') as of:
        blk_count = 0
        for line in inf:
            if 'BLOCK' in line:
                blk_count = 0
            if len(line) < 3 or 'BLOCK' in line or '****' in line:
                print(line,file=of,end='')
                continue

            el = line.strip().split()
            pruned_refhap_heuristic = int(el[8])
            split_conf = float(el[9]) if el[9] != '.' else 100
            snp_conf   = float(el[10]) if el[10] != '.' else 100

            if split_conf < split_conf_cutoff and blk_count >= 2:
                print('******** ',file=of)
                print('BLOCK: (from split)',file=of)

            if (use_refhap_heuristic and pruned_refhap_heuristic) or (snp_conf < snp_conf_cutoff):
                el[1] = '-'
                el[2] = '-'
                print('\t'.join(el),file=of)
            else:
                print(line,file=of,end='')

            blk_count += 1


def filter_wgs_vcf(in_vcf,out_vcfs,chroms):

    assert(len(chroms) == len(out_vcfs))

    dp_pat = re.compile("DP=(\d+)")
    qr_pat = re.compile(";QR=(\d+)")
    type_pat = re.compile("TYPE=([^;]+);")

    # step through WGS VCF and filter for high-confidence hetereozygous variants
    # if we are filtering on variants in another file (see above), also apply this filter.
    current_chrom = chroms.pop(0)
    current_vcf = out_vcfs.pop(0)
    outfile = open(current_vcf,'w')
    # add SNVs seen in WGS dataset to dictionary
    with open(in_vcf,'r') as vcf:
        for line in vcf:
            if line[0] == '#' or len(line) < 3:
                continue

            vcf_line = line.strip().split('\t')

            vcf_chrom = vcf_line[0]

            if vcf_chrom != current_chrom:
                outfile.close()
                if len(chroms) == 0:
                    break
                current_chrom = chroms.pop(0)
                current_vcf   = out_vcfs.pop(0)
                outfile = open(current_vcf,'w')

            vcf_pos = int(vcf_line[1])

            fields = vcf_line[7]
            pos_type = re.findall(type_pat,line)
            depth = int(float(re.findall(dp_pat,fields)[0]))
            if depth < 10:
                continue

            # we only care about reference and SNPs, no indels or MNPs
            if not(pos_type == [] or pos_type[0] == 'snp'):
                continue
            if not (len(vcf_line[3]) == 1 and len(vcf_line[4]) == 1):
                continue

            if vcf_line[8][0:2] != 'GT':
                continue
            genotype = vcf_line[9][0:3]
            qual = int(float(re.findall(qr_pat,fields)[0])) if genotype == '0/0' else float(vcf_line[5])

            if qual < 30:
                continue

            ref_allele = str.upper(vcf_line[3])
            variant_allele = str.upper(vcf_line[4][0])
            third_allele = 'N'
            if len(vcf_line[4]) == 3:
                third_allele = str.upper(vcf_line[4][2])

            alleles = set()

            if '0' in genotype:
                alleles.add(ref_allele)
            if '1' in genotype:
                alleles.add(variant_allele)
            if '2' in genotype:
                alleles.add(third_allele)

            # skip cases with N alleles, or with any other unusual cases
            if len(alleles.intersection(bases)) != len(alleles):
                continue

            # we only want heterozygous variants for haplotyping
            if genotype not in ['0/1','1/0','0|1','1|0']:
                continue

            # print the line
            print(line.strip(),file=outfile)

    if not outfile.closed:
        outfile.close()
