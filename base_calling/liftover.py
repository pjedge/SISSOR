from pyliftover import LiftOver

chroms = set(['chr{}'.format(i) for i in range(1,23)])

def liftover(input_vcfs, chunklist, output_vcfs):
    lo = LiftOver('hg38ToHg19.over.chain')

    regions_output = list(zip(chunklist, output_vcfs))
    (chrom, start, stop), outputfile = regions_output.pop(0)
    output = open(outputfile, 'w')

    for input_vcf in input_vcfs:
        with open(input_vcf,'r') as vcf:
            for line in vcf:
                if line[0] == '#' or len(line) < 3:
                    continue

                vcf_line = line.strip().split('\t')

                vcf_chrom = vcf_line[0]
                if vcf_chrom not in chroms:
                    continue
                vcf_pos_old   = int(vcf_line[1])
                converted = lo.convert_coordinate(vcf_chrom,vcf_pos_old)
                if converted == []:
                    continue

                vcf_pos = int(converted[0][1])
                ref_allele = vcf_line[3]
                variant_allele = vcf_line[4][0]
                third_allele = 'N'
                if len(vcf_line[4]) == 3:
                    third_allele = vcf_line[4][2]

                if vcf_chrom != chrom or vcf_pos > stop:
                    output.close()
                    if len(regions_output) == 0:
                        break
                    (chrom, start, stop), outputfile = regions_output.pop(0)
                    assert(vcf_chrom == chrom)
                    assert(vcf_pos >= start)
                    assert(vcf_pos <= stop)
                    output = open(outputfile, 'w')

                # write to lifted over vcf
                hg19_line_el = [vcf_line[0]]+[str(converted[0][1])]+vcf_line[2:]
                hg19_line = '\t'.join(hg19_line_el)
                print(hg19_line,file=output)

    if not output.closed:
        output.close()

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
                output.close()
                done = False
                while not (vcf_chrom == chrom and vcf_pos >= start and vcf_pos <= stop):
                    if len(regions_output) == 0:
                        done = True
                        break
                    (chrom, start, stop), outputfile = regions_output.pop(0)
                if done:
                    break
                assert(vcf_chrom == chrom)
                assert(vcf_pos >= start)
                assert(vcf_pos <= stop)
                output = open(outputfile, 'w')

            # write to lifted over vcf
            print(line,end='',file=output)

    if not output.closed:
        output.close()
