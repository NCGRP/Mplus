#include "m+.hpp"


/*
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
*/



/***************FUNCTIONS*****************/


vector<std::string> MySetKernel(char* KerFilePath)
{
	//declare variables
	std::string foo;
	vector<std::string> KernelAccessionList;
	
	//add mandatory accessions to list
	std::ifstream infile;
	infile.open(KerFilePath);
        while( !infile.eof() ) // To get all the lines.
       {
	        std::getline (infile, foo); // Saves the line in foo.
	        if (foo == "") break; //get out of while when foo string is empty
	        KernelAccessionList.push_back(foo);
	    }

	return KernelAccessionList;
}

//split a string on whitespace
vector<string> split(string const &input) 
{ 
    istringstream buffer(input);
    vector<string> ret;

    copy(std::istream_iterator<string>(buffer), 
              std::istream_iterator<string>(),
              back_inserter(ret));
    return ret;
}

//removes duplicate values from a list but does not alter sort order
vector<std::string> unsortedRemoveDuplicates(vector<std::string> numbers)
{
	vector<std::string> uniqvec;
	std::string b;
	for (unsigned int i=0;i<numbers.size();++i)
	{
		b = numbers[i];
		if (std::find( uniqvec.begin(), uniqvec.end(), b) == uniqvec.end() ) uniqvec.push_back(b); //number has not been seen, add it
	}
	return uniqvec;
}

//tabulates the ploidy for each locus
vector<int> GetPloidy(vector<std::string> AllLociNameList, vector<std::string> UniqLociNameList)
{
	int n;
	std::string b;
	vector<int> PloidyList;
	for (unsigned int i=0;i<UniqLociNameList.size();++i)
	{
		b = UniqLociNameList[i];
		n = std::count (AllLociNameList.begin(), AllLociNameList.end(), b);
		PloidyList.push_back(n);
	}
	
	return PloidyList;
}

int MyProcessVarFile(char* VarFilePath, vector<int>& AllColumnIDList, vector<std::string>& AllLociNameList, vector<int>& ActiveColumnIDList, vector<std::string>& ActiveLociNameList, vector<int>& TargetColumnIDList, vector<std::string>& TargetLociNameList, vector<vector<int> >& ColKeyToAllAlleleByPopList, vector<int>& ReferenceOrTargetKey, vector<int>& PloidyList, vector<std::string>& UniqLociNameList)
{
    //declare variables
    std::string foo;
    vector<std::string> foovector;
    unsigned int k;
    
    int i=0; // i is the row number
	std::ifstream infile;
	infile.open(VarFilePath);
    while( !infile.eof() ) // To get all the lines.
    {
	    std::getline (infile, foo); // Saves the line in foo.
	        
	    //split foo on whitespace
		foovector = split(foo);

		//identify active columns with qualitative data, classify those as reference or target
		if (foovector[1] == "2") //column holds qualitative data
		{
			AllColumnIDList.push_back(i);
			AllLociNameList.push_back(foovector[0]);
				
			if ((foovector[2] == "1") && (foovector[3] == "0")) //reference variable
			{
				ActiveColumnIDList.push_back(i);
				ActiveLociNameList.push_back(foovector[0]);
			}
			else if ((foovector[2] == "0") && (foovector[3] == "1"))  //target variable
			{
				TargetColumnIDList.push_back(i);
				TargetLociNameList.push_back(foovector[0]);
			}
		}
	     
		i++;
		foovector.clear();  //zero vector foovector
    }
	
	infile.close();
	
	
	//make a key showing which loci are target and which are reference
	//this will be used later to sort out the AllAlleleByPopList
	std::string uniqloc;
	std::string currloc;
	
	//find unique locus names, retains sort order
	UniqLociNameList = unsortedRemoveDuplicates(AllLociNameList);

	//tabulate the ploidy for each locus
	PloidyList = GetPloidy(AllLociNameList, UniqLociNameList);
	
	//define 2d vector that is key to columns, size to number of unique loci
	ColKeyToAllAlleleByPopList.resize( UniqLociNameList.size() ); //size to number of loci
	int b;
	for (unsigned int i=0;i<UniqLociNameList.size();++i)
	{
		uniqloc = UniqLociNameList[i];
		for (k=0;k<AllLociNameList.size();++k)
		{
			currloc = AllLociNameList[k];
			if (currloc == uniqloc)
			{
				b = AllColumnIDList[k]; //get the corresponding value in the columnID list
				ColKeyToAllAlleleByPopList[i].push_back(b);	
			}	
		}		
	}
	
	//define a 1d vector that describes whether the loci in the ColKey are reference(0) or target(1)
	double rt;
	ReferenceOrTargetKey.resize( ColKeyToAllAlleleByPopList.size() ); //sized to same length as key
	for (unsigned int i=0;i<ColKeyToAllAlleleByPopList.size();++i)
	{
		rt=0;
		//test whether all elements are categorized as reference or as target, if not, raise error
		for (k=0;k<ColKeyToAllAlleleByPopList[i].size();++k)
		{
			b=ColKeyToAllAlleleByPopList[i][k];
			if(std::find(ActiveColumnIDList.begin(), ActiveColumnIDList.end(), b) != ActiveColumnIDList.end())
			{
   			 	//column is reference
   			 	rt=rt+0;
			} 
			else if(std::find(TargetColumnIDList.begin(), TargetColumnIDList.end(), b) != TargetColumnIDList.end())
			{
    			//column is target
    			rt=rt+1;
			}
		}
		
		//test whether columns in key represent reference or target loci
		if (rt == 0) ReferenceOrTargetKey[i] = 0; //it is a reference
		else if (  rt/ColKeyToAllAlleleByPopList[i].size() == 1 ) ReferenceOrTargetKey[i] = 1; //it is a target
		else 
		{
			cout << "ERROR:  Some loci are described as both reference and target.  Please check your var file. Quitting...\n\n";
			exit (EXIT_FAILURE);
		}
	}
	return 0;
}

//quickly reads a large .dat file into a memory buffer
char * MyBigRead(char* DatFilePath)
{
	FILE * pFile;
	unsigned long long lSize;
	char * buffer;
	size_t result;

	pFile = fopen ( DatFilePath , "r" );
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	rewind (pFile);

	// allocate memory to contain the whole file:
	buffer = (char*) malloc (sizeof(char)*lSize);
	if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

	// copy the file into the buffer:
	result = fread (buffer,1,lSize,pFile);
	if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

	fclose (pFile);
	return buffer;
}

//adds an allele to the appropriate locus in the ActiveAlleleList
void MyUpdateActiveAlleleList(vector<vector<std::string> >& ActiveAlleleList, int CurrItemIndex, std::string NewAllele) //updates ActiveAlleleList by reference
{	
	vector<std::string> OldAlleles;
	
	//add to the current array being built or to the appropriate existing array
	OldAlleles = ActiveAlleleList[CurrItemIndex]; //pull out the list of alleles currently held for this locus
	OldAlleles.push_back (NewAllele); //add the new allele to the existing list of alleles for the current locus
	ActiveAlleleList[CurrItemIndex] = OldAlleles; //update ActiveAlleleList by reference
}

//reduces the master vector of all alleles into subsets containing reference or target loci only
int MyReduceToRef(vector<vector<vector<std::string> > > AllAlleleByPopList, vector<int> ReferenceOrTargetKey, vector<vector<vector<std::string> > >& ActiveAlleleByPopList, vector<vector<vector<std::string> > >& TargetAlleleByPopList)
{
	unsigned int r, t, i, j;
	vector<std::string> b;
	
	//resize level 1 of vectors, they contain the same number of populations as AllAlleleByPopList
	ActiveAlleleByPopList.resize(AllAlleleByPopList.size());
	TargetAlleleByPopList.resize(AllAlleleByPopList.size());
	
	//reserve level 2 of vectors, they contain the number of loci in ReferenceOrTargetKey
	r = std::count(ReferenceOrTargetKey.begin(), ReferenceOrTargetKey.end(), 0);
	t = std::count(ReferenceOrTargetKey.begin(), ReferenceOrTargetKey.end(), 1);
	for (i=0;i<AllAlleleByPopList.size();++i)
	{
		ActiveAlleleByPopList[i].reserve(r);
		TargetAlleleByPopList[i].reserve(t);
	}

	
	for (i=0;i<AllAlleleByPopList.size();++i) //iterate thru populations
	{
		for (j=0;j<AllAlleleByPopList[i].size();++j) //iterate thru loci
		{
			b.clear();
			b.resize( AllAlleleByPopList[i][j].size() );
			b = AllAlleleByPopList[i][j];
			
			if (ReferenceOrTargetKey[j] == 0) //it is a reference locus
			{
				ActiveAlleleByPopList[i].push_back(b);
			}
			
			else if (ReferenceOrTargetKey[j] == 1) //it is a target locus
			{
				TargetAlleleByPopList[i].push_back(b);
			}
		}
	}
	//ActiveAlleleByPopList & TargetAlleleByPopList have now been updated
	return 0;
}

int MyProcessDatFileIII(char* DatFilePath, vector<int> AllColumnIDList, vector<std::string> AllLociNameList, vector<vector<int> > ColKeyToAllAlleleByPopList, vector<vector<set<std::string> > >& AllAlleleByPopListSet, vector<std::string>& FullAccessionNameList, vector<std::string>& IndivPerPop, vector<std::string>& AllAlleles)
{
	//declare variables
    std::string foo;
    vector<std::string> foovector;
    std::string OldLocusName;
    std::string CurrLocusName;
    vector<std::string> LocusNames;
    vector<vector<std::string> > ActiveAlleleList;
    vector<std::string> OldAlleles;
    vector<vector<std::string> > TempList2d;
    vector<std::string> FilteredData;
    vector<std::string> ListToFilter;
    std::string IsNewPop = "no";
    vector<std::string>::iterator it;

    unsigned int i,j,k,l;
    vector<std::string> bufvec;
    std::string NewPopID;
    std::string OldPopID = "init*@#rt4"; //use an unlikely population name for the initialization value
    vector<std::string> TempList;
    std::string NewAllele;

	//read the whole file into a buffer using fread
    char * buffer;
	buffer = MyBigRead(DatFilePath);
	stringstream s(buffer); //put giant char array into a stream
	
    //read buffer into a vector, one line per item
    while (getline(s, foo))//get line from s, put in foo, consecutively
    {
	 	bufvec.push_back(foo);  
	}
	
	//sort vector so that individuals from the same population form consecutive elements
	std::sort(bufvec.begin(), bufvec.end()); //no need to use fancy sort, lexicographic should be fine
	
	//break up vector into a 3d vector by population:  { { {pop1ind1elems},{pop1ind2elems},...}, { {pop2ind1elems},{pop1ind2elems},...} }
	vector<vector<vector<std::string> > > ByPop3d;
	for (i=0;i<bufvec.size();++i)
	{
		TempList = split(bufvec[i]); //split line i on whitespace	
		NewPopID = TempList[0];
		IndivPerPop.push_back(NewPopID); //add the pop ID to a list to calc pop sizes later

		if (NewPopID != OldPopID) //then create a new population in ByPop2D
		{
			ByPop3d.resize(ByPop3d.size() + 1);
			
			//add the new population name to the AccessionNameList
			FullAccessionNameList.push_back (NewPopID);
		}
		
		//push vector TempList, containing elements on current line, onto last item of ByPop3d, which might be a new population as added just above
		ByPop3d[ByPop3d.size()-1].push_back(TempList);
		
		OldPopID = NewPopID; 
	}
	
	/*//print out ByPop3d
	for (i=0;i<ByPop3d.size();++i)
	{
		cout << "Pop" << i << "\n";
		for (j=0;j<ByPop3d[i].size();++j)
		{
			cout << " Ind" << j << "  ";
			for (k=0;k<ByPop3d[i][j].size();++k)
			{
				cout << ByPop3d[i][j][k] << ",";
			}
			cout << "\n";
		}
	
	}*/

	//resize AllAlleleByPopListSet
	AllAlleleByPopListSet.resize(ByPop3d.size());//resize number of populations
	for (i=0;i<AllAlleleByPopListSet.size();++i)
	{
		AllAlleleByPopListSet[i].resize(ColKeyToAllAlleleByPopList.size()); //resize number of loci
																		   //the index in ColKey is the locus index in AllAllelesByPopList level 2
																		   //the value of ColKey is the index of the allele in ByPop3d level 3
	}
	
	//condense alleles by locus in AllAllelesByPopList, within each population
	int AlleleIndex;
	for (i=0;i<ByPop3d.size();++i) //go thru pops
	{
		for (j=0;j<ByPop3d[i].size();++j) //go thru indivs
		{
			TempList = ByPop3d[i][j]; //get the list of column entries in an indiv
			for (k=0;k<ColKeyToAllAlleleByPopList.size();++k) //go through each locus
			{
				for (l=0;l<ColKeyToAllAlleleByPopList[k].size();++l) //assign columns to loci
				{
					AlleleIndex = ColKeyToAllAlleleByPopList[k][l];
					NewAllele = ByPop3d[i][j][AlleleIndex]; //get the allele in the specified column
					AllAlleles.push_back(NewAllele); //add the allele to the list of all alleles, missing data included
					if (NewAllele != "9999") //exclude missing data
						AllAlleleByPopListSet[i][k].insert(NewAllele); //add the allele to the set of unique alleles at locus k, pop i	
				}
			}
		}
	}
	return 0;
}

//removes duplicate alleles and missing data (9999) from the supplied vector
vector<std::string> MyFilterDuplicates(vector<std::string> ListToFilter)
{
	//remove duplicates
	sort( ListToFilter.begin(), ListToFilter.end() );
	ListToFilter.erase( std::unique( ListToFilter.begin(), ListToFilter.end() ), ListToFilter.end() );
	
	//remove missing data
	ListToFilter.erase( std::remove( ListToFilter.begin(), ListToFilter.end(), "9999" ), ListToFilter.end() );
	
	return ListToFilter;
}

//removes duplicate numbers from the supplied vector
vector<std::string> MyFilterDuplicatesII(vector<std::string> ListToFilter)
{
	//remove duplicates
	sort( ListToFilter.begin(), ListToFilter.end() );
	ListToFilter.erase( std::unique( ListToFilter.begin(), ListToFilter.end() ), ListToFilter.end() );
	
	return ListToFilter;
}

//returns maximum number of alleles possible at each locus for active and target
vector<int> MyGetMaxs(vector<vector<vector<std::string> > > ActiveAlleleByPopList)
{
	unsigned int i, j, k;
	vector<int> ActiveMaxAllelesList;
	vector<std::string> CurrLoc;
	set<std::string> NewSet;
	
	for (i=0;i<ActiveAlleleByPopList[0].size();++i)
	{
		NewSet.clear();
		for (j=0;j<ActiveAlleleByPopList.size();j++)
		{
			CurrLoc = ActiveAlleleByPopList[j][i]; //you are traversing the locus 'column' of the 3d grid
												   //get locus i for population j
			for (k=0;k<CurrLoc.size();++k)
			{
				NewSet.insert(CurrLoc[k]);	//place alleles into set to eliminate redundancies
			}
		}
		ActiveMaxAllelesList.push_back(NewSet.size()); //the set size after adding all populations for locus i is the maximum number of alleles
	}
	return ActiveMaxAllelesList; 
}

bool fileExists(const char *fileName)
{
    ifstream infile(fileName);
    return infile.good();
}

vector<std::string> MyRemoveTargetAlleles(vector<std::string> AllAlleles, vector<int> AllColumnIDList, vector<int> TargetColumnIDList)
{
	unsigned int checker, j;
	unsigned int i=0;
	vector<std::string> AllRefAlleles;
	while (i<AllAlleles.size())
	{
		j=0;
		while (j<AllColumnIDList.size()) //this routine follows along, mapping column ids to target or, implicitly, reference
		{
			checker = AllColumnIDList[j];		
			if(std::find(TargetColumnIDList.begin(), TargetColumnIDList.end(), checker) != TargetColumnIDList.end()) 
			{
				// TargetColumnIDList contains checker
				++j;
				++i;
			} 
			else 
			{
				//TargetColumnIDList does not contain checker
				AllRefAlleles.push_back(AllAlleles[i]);
				++j;
				++i;
			}
		}
	}
	return AllRefAlleles;
}	

//places a continuous string of alleles into a 2d vector with samples as rows
void MyMakeRefAllelesIntoRows(vector<std::string> AllRefAlleles, vector<std::string> ActiveLociNameList, vector<vector<std::string> >& RefAllelesIntoRows)
{
	unsigned int j; //counts items per row
	unsigned int i=0, k=0; //k is the row number
	while (i<AllRefAlleles.size())
	{
		j=0;
		while (j<ActiveLociNameList.size()) //keeps track of the number of alleles traversed before setting the next row
		{
			RefAllelesIntoRows[k][j] = AllRefAlleles[i];
			++i;
			++j;
		}
		++k;
	}
}	

//remove columns of data to make 2d vector of alleles by locus 
void MyMakeRefAllelesByLocus(vector<vector<std::string> > RefAllelesIntoRows, vector<string> ActiveLociNameList, vector<std::pair<std::string, vector<std::string> > >& RefAllelesByLocus)
{
	vector<std::pair<std::string, vector<std::string> > > lola;
	std::pair<std::string, vector<std::string> > la; //locus name, allele list pair
	
	unsigned int i, j, k;
	int locindex;
	std::string locname, b;
	vector<std::string> foo;
	for (i=0;i<ActiveLociNameList.size();++i)
	{
		locname = ActiveLociNameList[i];
		//test whether the list of locus names/alleles already contains this locus
		locindex = -1; //default value will cause an error if not explicitly set
		for (j=0;j<lola.size();++j)
		{
			b = lola[j].first;
			if (locname == b) //there is already data for that locus
			{
				locindex = j;
			}
		}
		if (locindex == -1)  //locname not found, add a pair
		{
			locindex = lola.size(); 
			lola.push_back(la);//add an empty pair onto the vector of loci
			lola[locindex].first = locname;
		}
		
		//add the column of data defined by locname to the appropriate pair.second, defined by locindex in lola
		for (k=0;k<RefAllelesIntoRows.size();++k)
		{
			b = RefAllelesIntoRows[k][i];
			if (b != "9999") lola[locindex].second.push_back(b); //ignore missing data
		}
	}
	
	//update RefAllelesByLocus
	RefAllelesByLocus = lola;

	//print out each locus name followed by all the alleles found within it
	/*for (i=0;i<lola.size();++i)
	{
		cout << lola[i].first << "\n";
		foo = lola[i].second;
		for (j=0;j<foo.size();++j)
		{
			cout << " " << foo[j];
		}
		cout << "\n";
	}
	*/
}

//calculates allele frequencies for all alleles at all loci, updates vector of struct Alfreq, which contains the relational data
void MyCalculateAlleleFrequencies(vector<std::pair<std::string, vector<std::string> > > RefAllelesByLocus, vector<Alfreq>& AlleleFrequencies)
{
	unsigned int i, j, z;
	double freq;
	std::string b;
	vector<std::string> AllAlleles;
	vector<std::string> UniqAlleles;
	set<std::string> AlleleSet;

	Alfreq laf; //locusname, allelenames, frequencies
	static const struct Alfreq emptylaf; //this will be used to zero struct between loops
	
	for (i=0;i<RefAllelesByLocus.size();++i)
	{
		//empty containers
		laf = emptylaf;
		vector<std::string>().swap(UniqAlleles); //clear UniqAlleles
		AlleleSet.clear(); //clear AlleleSet
		
		//get locus name, add to struct
		b = RefAllelesByLocus[i].first;
		laf.locusname = b;
		
		//compress list of all alleles at this locus into unique alleles
		AllAlleles = RefAllelesByLocus[i].second; //get the vector of all alleles at locus i
		for (j=0;j<AllAlleles.size();++j) 
		{
			AlleleSet.insert( AllAlleles[j] ); //filter out redundant alleles by dumping the vector into a set
		}
		UniqAlleles.assign(AlleleSet.begin(), AlleleSet.end()); //assign the unique alleles to a vector
		
		//for each allele in UniqAlleles count the number of occurrences in AllAlleles, calc frequency, add allele name and freq to struct
		for (j=0;j<UniqAlleles.size();++j)
		{
			b = UniqAlleles[j];
			laf.allelenames.push_back(b); //add allele name to struct
			z = count(AllAlleles.begin(), AllAlleles.end(), b);
			freq = double(z)/double(AllAlleles.size());
			laf.frequencies.push_back(freq); //add allele frequency to struct
		}
		
		AlleleFrequencies.push_back(laf);
	}
}


/***************MAIN*****************/

int main( int argc, char* argv[] )
{
	//get mandatory command line arguments
	char* VarFilePath = argv[1];
	char* DatFilePath = argv[2];
	
	//initialize optional command line arguments
	unsigned int MinCoreSize;
	unsigned int MaxCoreSize;
	int SamplingFreq;
	int NumReplicates;
	char* OutFilePath;
	char* KerFilePath;
	char* IdealFilePath;

	//declare variables
	unsigned int i, j;
	int b;
	string DoM = "no"; //switch to perform M+ optimization
	string Kernel = "no"; //switch to include a mandatory set in the core
	string Ideal = "no"; //switch to compute the ideal core, using A* algorithm
	vector<std::string> KernelAccessionList;
	vector<std::string> BadFiles;
	string bf;

	//parse the command line for options
	for (int i=0;i<argc;i++)
	{
		if ( string(argv[i]) == "-m" ) 
    	{
        	DoM = "yes";
        	MinCoreSize = atoi(argv[i+1]);
			MaxCoreSize = atoi(argv[i+2]);
			SamplingFreq = atoi(argv[i+3]);
			NumReplicates = atoi(argv[i+4]);
			OutFilePath = argv[i+5];
		}

		if ( string(argv[i]) == "-k" ) 
    	{
        	Kernel = "yes";
        	KerFilePath = argv[i+1];
        	KernelAccessionList = MySetKernel(KerFilePath);
			//verify that specified input file actually exists
			if (fileExists(KerFilePath) == 0) 
			{
				bf = "KerFilePath = ";
				bf += KerFilePath;
				BadFiles.push_back(bf);
			}
		}
		
		if ( string(argv[i]) == "-a" ) 
    	{
        	Ideal = "yes";
        	IdealFilePath = argv[i+1];
		}
	}
	
	//test whether all files specified on the command line exist
	if (fileExists(VarFilePath) == 0) 
	{
		bf = "VarFilePath = ";
		bf += VarFilePath;
		BadFiles.push_back(bf);
	}
	if (fileExists(DatFilePath) == 0) 
	{
		bf = "DatFilePath = ";
		bf += DatFilePath;
		BadFiles.push_back(bf);
	}
	
	if (BadFiles.size() > 0)
	{
		cout << "\nThe following variables appear to contain misspecified paths:\n";
		for (i=0;i<BadFiles.size();++i)
		{
			cout << "  " << BadFiles[i] << "\n";
		}
		cout << "\nPlease check the command line.  Quitting...\n\n";
		exit (EXIT_FAILURE);
	}
	
	//print out input variables
	cout << "\nInput variables:\n  VarFilePath = " << VarFilePath << "\n";
	cout << "  DatFilePath = " << DatFilePath << "\n";
	if (DoM == "yes")
	{
		cout << "  -m invoked:\n";
		cout << "    MinCoreSize = " << MinCoreSize << "\n";
		cout << "    MaxCoreSize = " << MaxCoreSize << "\n";
		cout << "    SamplingFreq = " << SamplingFreq << "\n";
		cout << "    NumReplicates = " << NumReplicates << "\n";
		cout << "    OutFilePath = " << OutFilePath << "\n";
	}
	if (Kernel == "yes") 
	{
		cout << "  -k invoked:\n";
		cout << "    KerFilePath = " << KerFilePath << "\n";
	}
	if (Ideal == "yes")
	{
		cout << "  -a invoked:\n";
		cout << "    IdealFilePath = " << IdealFilePath << "\n";
	}
	
	//catch some errors in the command line
	if (MinCoreSize > MaxCoreSize) 
	{
		cout << "ERROR:  The minimum core size ("<<MinCoreSize<<") is greater than the maximum core size.  Please correct the command line arguments.  Quitting...\n\n";
    	exit (EXIT_FAILURE);
    }
    if (MinCoreSize == 1)
    {
    	cout << "ERROR:  A minimum core size of 1 is not allowed.  Please correct the command line.  Quitting...\n\n"; 
		exit (EXIT_FAILURE);
	}
	if (Kernel == "yes" && Ideal == "yes" && DoM == "no")
	{
		cout << "ERROR:  A* search cannot be performed with a kernel file.  Please correct the command line.  Quitting...\n\n";
		exit (EXIT_FAILURE);
	}
    if (KernelAccessionList.size() > MinCoreSize)
    {
		cout << "ERROR:  The number of mandatory accessions ("<<KernelAccessionList.size()
			 <<") is greater than the minimum core size ("<<MinCoreSize
        	 <<").  Please modify the kernel file or the command line argument.  Quitting...\n\n";
			exit (EXIT_FAILURE);
	}
	

	//DETERMINE MACHINE CONFIGURATIION
	int ncpu = sysconf( _SC_NPROCESSORS_ONLN );
	
	//PROCESS INPUT DATA
	
	//start the clock
	time_t starti,endi;
	time (&starti);

	//.var file
	vector<int> AllColumnIDList;
	vector<std::string> AllLociNameList;
	vector<int> ActiveColumnIDList;
	vector<std::string> ActiveLociNameList;
    vector<int> TargetColumnIDList;
	vector<std::string> TargetLociNameList;
	vector<vector<int> > ColKeyToAllAlleleByPopList;
	vector<int> ReferenceOrTargetKey;
	vector<int> PloidyList;
	vector<std::string> UniqLociNamesList;

	MyProcessVarFile(VarFilePath, AllColumnIDList, AllLociNameList, ActiveColumnIDList, ActiveLociNameList, TargetColumnIDList, TargetLociNameList, ColKeyToAllAlleleByPopList, ReferenceOrTargetKey, PloidyList, UniqLociNamesList ); 
	//all but first variable above are updated as references in MyProcessVarFile

		//Print out ActiveColumnIDList and ActiveLociNameList
		if (ActiveColumnIDList.size() == 0)
		{
			cout << "\nThere are no active loci.\n";
		}
		else
		{
			cout << "\nActive Locus	Column\n";
			for (i=0; i < ActiveColumnIDList.size(); i++)
			{
				cout << ActiveLociNameList[i] << "\t" << (ActiveColumnIDList[i] + 1) << "\n";
			}
		}
		
		//Print out TargetColumnIDList and TargetLociNameList
		if (TargetColumnIDList.size() == 0)
		{
			cout << "\nThere are no target loci.\n";
		}
		else
		{
			cout << "\nTarget Locus	Column\n";
			for (i=0; i < TargetColumnIDList.size(); i++)
			{
				cout << TargetLociNameList[i] << "\t" << (TargetColumnIDList[i] + 1) << "\n";
			}
		}

	//process .dat file
	cout << "\nProcessing .dat file...\n";
	
	vector<std::string> FullAccessionNameList;
	vector<std::string> IndivPerPop;
	vector<std::string> AllAlleles;

	//switch for new MyProcessDatFileIII
	vector<vector<set<std::string> > > AllAlleleByPopListSet; //structure of this 3D vector is:
	// { { {pop1,loc1 alleles},{pop1,loc2 alleles},...}, { {pop2,loc1 alleles},{pop2,loc2 alleles},...} } }
	MyProcessDatFileIII(DatFilePath, AllColumnIDList, AllLociNameList, ColKeyToAllAlleleByPopList, AllAlleleByPopListSet, FullAccessionNameList, IndivPerPop, AllAlleles);

	vector<vector<vector<std::string> > > AllAlleleByPopList( AllAlleleByPopListSet.size(), vector<vector<std::string> >(UniqLociNamesList.size()) ); //structure of this 3D vector is:
	// { { {pop1,loc1 alleles},{pop1,loc2 alleles},...}, { {pop2,loc1 alleles},{pop2,loc2 alleles},...} } }
	//sized to number of populations, number of loci.  last level left unsized
	
		/*//Print out lists of unique All alleles from AllAlleleByPopListSet
		set<std::string> si;
		for (i=0;i<AllAlleleByPopListSet.size() ;i++)
		{
			cout << "Population " << FullAccessionNameList[i] << "\n";
			for (j=0;j<AllAlleleByPopListSet[i].size();j++)
			{
				cout << "Locus " << j << "\n";
				si = AllAlleleByPopListSet[i][j];
				for (std::set<std::string>::iterator it=si.begin(); it!=si.end(); ++it)
					cout << *it << ",";
				cout << "\n";
			}
		}*/

	//convert set to vector for further processing
	for (i=0;i<AllAlleleByPopList.size();++i)
	{
		for (j=0;j<AllAlleleByPopList[i].size();++j)
		{
			vector<std::string> ttvec(AllAlleleByPopListSet[i][j].begin(), AllAlleleByPopListSet[i][j].end()); //use constructor to convert set to vector	
			AllAlleleByPopList[i][j] = ttvec;
		}
	}
	
	//sort AllAlleleByPopList into reference and target loci lists
	vector<vector<vector<std::string> > > ActiveAlleleByPopList; 
	vector<vector<vector<std::string> > > TargetAlleleByPopList; 
	MyReduceToRef(AllAlleleByPopList, ReferenceOrTargetKey, ActiveAlleleByPopList, TargetAlleleByPopList); //latter 2 variables updated as reference
			
		/*	
		//Print out FullAccessionNameList
		cout << "\n\nPopulation names\n";
		for (i=0; i<FullAccessionNameList.size();i++) cout << FullAccessionNameList[i] << "\n";
		
		//Print out lists of unique reference alleles from ActiveAlleleByPopList
		cout << "ActiveAlleleByPopList:\n";
		for (i=0; i<ActiveAlleleByPopList.size() ;i++)
		{
			cout << "Population " << i << "\n";
			for (j=0;j<ActiveAlleleByPopList[i].size();j++)
			{
				cout << "Locus " << j << "\n";
				for (k=0;k<ActiveAlleleByPopList[i][j].size();k++)
				{
					cout << ActiveAlleleByPopList[i][j][k] << ",";
				}
				cout << "\n";
			}
		
		}
		
		//Print out lists of unique target alleles from TargetAlleleByPopList
		cout << "TargetAlleleByPopList:\n";
		for (i=0; i<TargetAlleleByPopList.size() ;i++)
		{
			cout << "Population " << i << "\n";
			for (j=0;j<TargetAlleleByPopList[i].size();j++)
			{
				cout << "Locus " << j << "\n";
				for (k=0;k<TargetAlleleByPopList[i][j].size();k++)
				{
					cout << TargetAlleleByPopList[i][j][k] << ",";
				}
				cout << "\n";
			}
		}
		*/
		


	//CALCULATE SOME USEFUL VARIABLES
	//get total number of accessions
	unsigned int NumberOfAccessions = ActiveAlleleByPopList.size();
	
	//catch a possible error in the command line
	if (MaxCoreSize > NumberOfAccessions)
	{
		cout << "ERROR:  Maximum core size of "<<MaxCoreSize<< " is greater than the number of accessions.  Please set to a lower value. Quitting...\n\n";
		exit (EXIT_FAILURE);
	}

	//get number of loci
	int NumLoci = ActiveAlleleByPopList[1].size();
	
	//calculate population sizes
	vector<int> PopSizes;
	for (i=0;i<FullAccessionNameList.size();++i)
	{
		bf = FullAccessionNameList[i];
		b = std::count (IndivPerPop.begin(), IndivPerPop.end(), bf);
		PopSizes.push_back(b); //synchronized with FullAccessionNameList
		//cout << "b="<<b<<"\n";	
	}
	
	//calculate allele frequencies
	
	//1. remove target alleles from AllAlleles
	vector<std::string> AllRefAlleles = MyRemoveTargetAlleles(AllAlleles, AllColumnIDList, TargetColumnIDList);
	
	//2. put AllRefAlleles into a 2d vector, samples are rows
	vector<vector<std::string> > RefAllelesIntoRows(IndivPerPop.size(), vector<std::string> ( ActiveLociNameList.size() ));
	MyMakeRefAllelesIntoRows(AllRefAlleles, ActiveLociNameList, RefAllelesIntoRows);
	
	//3. extract alleles into vector of pairs (ignores missing data 9999)
	vector<std::pair<std::string, vector<std::string> > > RefAllelesByLocus; // first = locus name, second = vector of all alleles present, updated as reference below
	MyMakeRefAllelesByLocus(RefAllelesIntoRows, ActiveLociNameList, RefAllelesByLocus);
	
	//4. calculate allele frequencies, finally (ignores missing data 9999)
	vector<Alfreq> AlleleFrequencies;//(RefAllelesByLocus.size());  //declare vector of struct Alfreq
	MyCalculateAlleleFrequencies(RefAllelesByLocus, AlleleFrequencies);
	
		/*
		//print out structs containing allele frequencies
		Alfreq laf;
		vector<std::string> anames;
		vector<double> frqs;
		for (i=0;i<AlleleFrequencies.size();++i)
		{
			cout << AlleleFrequencies[i].locusname << "\n";
		
			laf = AlleleFrequencies[i];
			anames = AlleleFrequencies[i].allelenames;
			frqs = AlleleFrequencies[i].frequencies;
			for (j=0;j<anames.size();++j)
			{
				cout << " " << anames[j];
				cout << " " << frqs[j] << "\n";
			}
		}
		*/
	
	//from list of all accession names, generate an index that is used later to locate alleles belonging to specific accessions
	vector<int> AccessionNameList;
	for (i=0;i<NumberOfAccessions;i++) AccessionNameList.push_back (i);
	//make a similar list for only those accessions contained in the kernel
	vector<int> KernelAccessionIndex;
	int breaker;//contains a switch value to exit the inner loop
	if (Kernel == "yes") 
	{
		for (i=0;i<KernelAccessionList.size();i++)
		{
			breaker = 0;
			for (j=0;j<FullAccessionNameList.size();j++)
			{
				//find the name in the KernelAccessionList in the FullAccessionNameList
				if (FullAccessionNameList[j] == KernelAccessionList[i])
				{
					//add the index value in the coordinated AccessionNameList to the KernelAccessionIndex
					KernelAccessionIndex.push_back(AccessionNameList[j]);
					breaker = 1;
				}
				if (breaker == 1) break;
			}
		}
	}
	//reverse sort KernelAccessionIndex, necessary later
	std::sort(KernelAccessionIndex.begin(), KernelAccessionIndex.end(), std::greater<int>());
		
	//get maximum number of alleles possible at each locus for active and target
	vector<int> ActiveMaxAllelesList, TargetMaxAllelesList;
	ActiveMaxAllelesList = MyGetMaxs(ActiveAlleleByPopList);
	TargetMaxAllelesList = MyGetMaxs(TargetAlleleByPopList);
		
	//stop the clock
	time (&endi);
	double dif = difftime (endi,starti);
	cout << "Input files processed.  Elapsed time = "<< dif << " seconds.\n\n";
	
	cout << "Number of accessions = " << NumberOfAccessions << ", Number of loci = " << NumLoci;
	if (Kernel == "yes") cout << ", Number of kernel accessions = " << KernelAccessionList.size() << "\n\n";
	else cout << "\n\n";
	
	
	//PERFORM A*
	if (Ideal == "yes")
	{
		int parallelism_enabled = 1; //0=no, not 0 = yes
		if (parallelism_enabled == 0) cout << "\nBeginning serial A* search...\n\n";
		else cout << "\nBeginning parallel A* search...\n\n";
				
		//start the clock
		time_t start1,end1;
		time (&start1);

		//run A*
		aStar
		(
			IdealFilePath, 
			ActiveAlleleByPopList, 
			ActiveMaxAllelesList, 
			UniqLociNamesList, 
			ReferenceOrTargetKey, 
			FullAccessionNameList, 
			PloidyList, 
			PopSizes, 
			AlleleFrequencies, 
			parallelism_enabled,
			ncpu
		);
		
		//stop the clock
		time (&end1);
		double dif = difftime (end1,start1);
		cout << "\nA* search complete.  Elapsed time = "<< dif << " seconds.\n\n";
	}	
	
	//PERFORM M+
	if (DoM == "yes")
	{
		//compile as parallel or not?
		int parallelism_enabled = 1; //0=no, not 0 = yes
		if (parallelism_enabled == 0) cout << "\nBeginning serial M+ search...\n\n";
		else cout << "\nBeginning parallel M+ search...\n\n";

		//start the clock
		time_t startm,endm;
		time (&startm);
		
		//run M+
		mp
		(
		MinCoreSize,
		MaxCoreSize,
		SamplingFreq,
		NumReplicates,
		OutFilePath,
		Kernel,
		KernelAccessionIndex,
		AccessionNameList,
		ActiveAlleleByPopList,
		TargetAlleleByPopList,
		ActiveMaxAllelesList,
		TargetMaxAllelesList,
		FullAccessionNameList,
		parallelism_enabled,
		ncpu
		);
		
		//stop the clock
		time (&endm);
		dif = difftime (endm,startm);
		cout << "\nM+ search complete.  Elapsed time = "<< dif << " seconds.\n\n";
	
	}
	
	return 0;
}
