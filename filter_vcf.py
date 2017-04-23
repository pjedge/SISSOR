
import re

bases = {'A','T','G','C'}
def filter_vcf(in_vcfs,out_vcfs, filter_set):

    vcf_dict = dict()

    with open(filter_set,'r') as infile:
        for line in infile:

            if line[:1] == '#':
                continue
            el = line.strip().split('\t')
            if len(el) < 5:
                continue

            chrom = el[0]
            genomic_pos = int(el[1])
            a1 = el[3]
            a2 = el[4]
            if el[9] not in ['0/1','1/0','1|0','0|1']:
                continue
            if len(a1) > 1 or len(a2) > 1:
                continue
            if a1 not in bases or a2 not in bases:
                continue

            vcf_dict[(chrom,genomic_pos)] = {a1,a2}

    for in_vcf, out_vcf in zip(in_vcfs,out_vcfs):
        with open(in_vcf,'r') as inf, open(out_vcf,'w') as outf:
            for line in inf:
                if line[:1] == '#':
                    continue
                el = line.strip().split('\t')
                if len(el) < 5:
                    continue

                chrom = el[0]
                genomic_pos = int(el[1])
                a1 = el[3]
                a2 = el[4]
                if el[9] not in ['0/1','1/0','1|0','0|1']:
                    continue
                if len(a1) > 1 or len(a2) > 1:
                    continue
                if a1 not in bases or a2 not in bases:
                    continue

                if (chrom,genomic_pos) not in vcf_dict or vcf_dict[(chrom,genomic_pos)] != {a1,a2}:
                    continue

                print(line.strip(),file=outf)

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
