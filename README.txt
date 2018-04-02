To compile:  use "make"
Usage: m+1 varfile datfile [-m mincoresize maxcoresize samplingfreq reps outputfile] [-r]
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
-r use rarefaction to correct for differences in sample size of accessions, applies to
        M+ algorithm only.
-k kernelfile = use an MSTRAT .ker file to specify mandatory members of the 
        core.  The number of mandatory accessions must therefore be less than or equal to 
        mincoresize.  Option only applies to -m, and cannot be used with -a.
-a idealcorefile = compute the minimum set of accessions necessary to retain all variation,
        i.e. the "ideal" or "best" core, using the A* search algorithm, write output to 
        idealcorefile.

Notes:  Missing data must be coded as 9999. To validate input files, omit all options. m+1
        uses OpenMP for parallelization on multicore machines. Default is single core M+
        algorithm, and multicore A* algorithm.  To change behavior, modify variable
        'parallelism_enabled' in file m+.cpp.

Examples:
          ./m+1 ./beet.var ./beet.dat -m 3 28 2 3 ./beetout.txt
          ./m+1 ./beet.var ./beet.dat -m 3 28 2 3 ./beetoutk.txt -k beet.ker
          ./m+1 ./beet.var ./beet.dat -a beetideal.txt
          ./m+1 ./orientalis.var ./orientalisIND.dat -m 2 50 1 1 orINDout.txt
          ./m+1 ./orientalis.var ./orientalisIND.dat -a orINDidealout.txt
          ./m+1 ./WheatSNP.var ./WheatSNP.dat -m 20 21 1 20 ./WheatSNPout.txt
          ./m+1 ./At.var ./At.dat -m 2 10 1 1 ./Atout.txt
          ./m+1 ./At.var ./At.dat -m 2 10 1 1 ./Atoutr.txt -r

