To compile:  use "make"
Usage: m+ varfile datfile [-m mincoresize maxcoresize samplingfreq reps outputfile]
		[-k kernelfile] [-a idealcorefile]
where, 
varfile = path to MSTRAT .var file
datfile = path to MSTRAT .dat file

Options:
-m mincoresize maxcoresize samplingfreq reps outputfile = compute the optimal accessions
		for a given core size using the M+ search algorithm. Arguments are as follows:
		mincoresize maxcoresize = integers specifying minimum and maximum core size. 
			Usually mincoresize = 2, maxcoresize = total number of accessions
		samplingfreq = e.g. integer value 5 will cause coresize=2 then coresize=7, then 12, 
			and so on, to be sampled.
		reps = number of replicate core sets to calculate for a particular core size
		outputfile = path to output
-k kernelfile = use an MSTRAT .ker file to specify mandatory members of the 
        core.  The number of mandatory accessions must therefore be less than or equal to 
        mincoresize.  Option must be used with -m.
-a idealcorefile = compute the minimum set of accessions necessary to retain all variation,
		i.e. the "ideal" or "best" core, using the A* search algorithm, write output to 
		bestcorefile.

Notes:  All input files must have Unix line breaks.

example: ./m+ ./beet.var ./beet.dat -m 3 28 2 3 ./beetout.txt -k beet.ker -a beetideal.txt
