import estimate_parameters as ep
infile = 'pileups/test_sample.pileup'
outfile = 'output/test_sample.out'
ep.estimate_parameters(infile,'testsample')
ep.combine_parameters(['testsample'])

import call_chamber_alleles as cca
cca.call_chamber_alleles(infile,outfile)