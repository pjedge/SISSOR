
bases = {'A','T','G','C'}
def filter_vcf(in_vcfs,out_vcfs,filter_set):

    vcf_dict = dict()
    with open(vcf_file,'r') as infile:
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

                if (chrom,genomic_pos) not in vcf_dict or vcf_dict[(chrom,genomic_pos)] != {a1,a2}:
                    continue

                print(line.strip(),file=outf)
