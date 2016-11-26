
# this takes the genome mappability score (GMS) files and splits them into
# convenient 5 Mb chunks, corresponding to regions called by sissorhands
def split_gms(infiles, chunklist, outputlst):

    chunks_output = list(zip(chunklist, outputlst))

    for infile in infiles:
        (chrom, start, stop), outputfile = chunks_output.pop(0)
        assert(start == 1)              # should be beginning of chrom
        assert(chrom+'.gms' in infile)  # file should have chrom name
        output = open(outputfile,'w')
        with open(infile,'r') as gms:
            for line in gms:
                if line[0] == '#' or len(line) < 3:
                    continue
                el = line.strip().split('\t')

                gms_chrom     = el[0]
                gms_pos       = int(el[1])

                if gms_pos > stop:
                    output.close()
                    if chunks_output == []:
                        break
                    (chrom, start, stop), outputfile = chunks_output.pop(0)
                    output = open(outputfile,'w')

                # make sure we're really on the current region
                assert(gms_chrom == chrom)
                assert(gms_pos >= start)
                assert(gms_pos <= stop)

                print(line,end='',file=output)

    if not output.closed:
        output.close()
