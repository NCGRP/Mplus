To compile:  g++ m+.cpp -o m+ -fopenmp
Usage: m+ varfile datfile mincoresize maxcoresize samplingfreq reps outputfile
where, 
varfile = path to MSTRAT .var file with unix line breaks
datfile = path to MSTRAT .dat file 
mincoresize maxcoresize = integers specifying minimum and maximum core size. 
          Usually mincoresize = 2, maxcoresize = total number of accessions
samplingfreq = e.g. integer value 5 will cause coresize=2 then coresize=7, then 12, 
          and so on, to be sampled.
reps = number of times to repeat the measurement of diversity for a particular core size 
          before calculating a mean
outputfile = path to output

Options:
-s summaryfile = compute summary statistics and write to an output file
          with the path summaryfile
-k kernelfile = use an MSTRAT .ker file to specify mandatory members of the 
          core.  The number of included accessions must be less than or equal to mincoresize.

example: ./m+ ./beet.var ./beet.dat 3 28 2 3 ./beetout.txt -s ./beetsum.txt -k beet.ker
