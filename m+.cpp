/*#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include <numeric>
#include <time.h>
#include <omp.h>
#include <algorithm>
*/
#include "m+.hpp"
//using namespace std;


/*
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
-b idealcorefile = compute the minimum set of accessions necessary to retain all variation, i.e.
		  the "ideal" or "best" core, write output to bestcorefile.  This option ignores the kernel file.
-k kernelfile = use an MSTRAT .ker file to specify mandatory members of the 
          core.  The number of included accessions must be less than or equal to mincoresize.

example: ./m+ ./beet.var ./beet.dat 3 28 2 3 ./beetout.txt -s ./beetsum.txt -b beetideal.txt -k beet.ker
*/



/***************FUNCTIONS*****************/


vector<std::string> MySetKernel(char* KerFilePath)
{
	//declare variables
	int i=0;
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
	for (int i=0;i<numbers.size();++i)
	{
		b = numbers[i];
		if (std::find( uniqvec.begin(), uniqvec.end(), b) == uniqvec.end() ) uniqvec.push_back(b); //number has not been seen, add it
	}
	return uniqvec;
}

//tabulates the ploidy for each locus
vector<int> GetPloidy(vector<std::string> AllLociNameList, vector<std::string> UniqLociNameList)
{
	int i, c, n;
	std::string b;
	vector<int> PloidyList;
	for (i=0;i<UniqLociNameList.size();++i)
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
    int k;
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
	for (i=0;i<UniqLociNameList.size();++i)
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
	for (i=0;i<ColKeyToAllAlleleByPopList.size();++i)
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
	//print out ColKeyToAllAlleleByPopList
	cout << "ColKeyToAllAlleleByPopList\n";
	for (i=0;i<ColKeyToAllAlleleByPopList.size();++i)
	{
		cout << "Locus "<<UniqLociNameList[i]<< " {";
		for (k=0;k<ColKeyToAllAlleleByPopList[i].size();++k)
		{
			cout<<ColKeyToAllAlleleByPopList[i][k]<<",";
		}
		cout << "\n";
	}
	/*
	
	//print out ReferenceOrTargetKey
	for (i=0;i<ReferenceOrTargetKey.size();++i)
	{
		cout << ReferenceOrTargetKey[i] << "\n";
	}	
	
	getchar();
	*/
	return 0;
}

//quickly reads a large .dat file into a memory buffer
char * MyBigRead(char* DatFilePath)
{
	FILE * pFile;
	long lSize;
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

	/* the whole file is now loaded in the memory buffer. */
	//foo(buffer, result);
	// terminate

	fclose (pFile);
	//free (buffer);
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
	
	/*cout << "OldAlleles\n";
	for (j=0;j<OldAlleles.size();j++)
	{
		cout << OldAlleles[j] << ",";
	}
	cout << "\n";
	*/
}

int MyProcessDatFile(char* DatFilePath, vector<int> ActiveColumnIDList, vector<std::string> ActiveLociNameList, vector<vector<vector<std::string> > >& ActiveAlleleByPopList, vector<std::string>& FullAccessionNameList, vector<std::string>& IndivPerPop, vector<std::string>& AllAlleles)
{


//declare variables
    std::string foo;
    vector<std::string> foovector;
    std::string NewPopID;
    std::string OldPopID="0";
    std::string OldLocusName;
    std::string NewAllele;
    std::string CurrLocusName;
    vector<std::string> LocusNames;
    int i=0; //i is the row number
    int j=0;
    int q=0;
    int CurrItemIndex;
    vector<vector<std::string> > ActiveAlleleList;
    vector<std::string> OldAlleles;
    vector<std::string> TempList;
    vector<vector<std::string> > TempList2d;
    vector<std::string> FilteredData;
    vector<std::string> ListToFilter;
    std::string IsNewPop = "no";
    vector<std::string>::iterator it;



	//read in the dat file line by line	
	std::ifstream infile;
	infile.open(DatFilePath);
    while( !infile.eof() ) // To get all the lines.
    {
	    if (IsNewPop == "no") //get a new line
	    {
			vector<std::string>().swap(foovector); //clear foovector
			foo.clear(); //clear foo
			
			std::getline (infile, foo); // Saves the line in foo.
			if (infile.eof()) break;//this is a crappy workaround to exit before extra loop after eof is encountered

			//split on whitespace
			foovector = split(foo);

			//extract population identifier from first column
			NewPopID = foovector[0];
			IndivPerPop.push_back(foovector[0]); //add the pop ID to a list to calc pop sizes later
			
		}
		else if (IsNewPop == "yes") IsNewPop = "no"; //reset it to get a new line unless it is a new allele
	
		
		//determine whether the current line belongs to the current population, or is the first from a new population
		if (NewPopID == OldPopID)
		{
			//extract the alleles, by locus, for a single population into a list
			OldLocusName = "starting value";
			for (i=0;i<ActiveColumnIDList.size();i++)
			{
				//get the current allele
				NewAllele = foovector[ActiveColumnIDList[i]];
				
				AllAlleles.push_back(NewAllele); //add the allele to the list of all alleles
								
				//get the name of the current locus
				CurrLocusName = ActiveLociNameList[i];
				
				/*
				cout << "ActiveColumnID=" << ActiveColumnIDList[i] << "\t";
				cout << "NewAllele=" << NewAllele << "\t";
				cout << "CurrLocusName=" << CurrLocusName << "\n";
				*/
				
				//test whether you already have a list item for this locus, 
				//count command tests whether CurrLocusName appears in LocusNames
				it = std::find(LocusNames.begin(), LocusNames.end(), CurrLocusName);
				if (it != LocusNames.end()) //you have found a match, convert iterator to index, update
				{
					CurrItemIndex = std::distance(LocusNames.begin(), it); //convert iterator to index
					MyUpdateActiveAlleleList(ActiveAlleleList, CurrItemIndex, NewAllele); //update ActiveAlleleList with the new allele by reference
				}
				
				/*
				if ( std::count (LocusNames.begin(), LocusNames.end(), CurrLocusName) >= 1) //you are adding to an existing locus
				{
					//figure out what gene you are in, then add NewAllele to the appropriate item in the ActiveAlleleList
					for (j=0; j<LocusNames.size(); j++)
					{
						if (LocusNames[j] == CurrLocusName)
						{
							CurrItemIndex = j;
							break;	
						}
					}
					MyUpdateActiveAlleleList(ActiveAlleleList, CurrItemIndex, NewAllele); //update ActiveAlleleList with the new allele by reference
				}
				
				*/
				else //it is new locus
				{
					//create a new list item, and add current allele to it. other alleles will be added later in first part of 'if'
					TempList.push_back (NewAllele); //add the single allele to a vector
					ActiveAlleleList.push_back (TempList);  //add the vector as a new list item
					vector<std::string>().swap(TempList); //clear TempList
					CurrItemIndex = ActiveAlleleList.size();
					
					//build a list of the locus names
					LocusNames.push_back (ActiveLociNameList[i]);
				}
				
				//pass current locus name to old locus name
				OldLocusName = CurrLocusName;
				CurrLocusName = ""; //zero it out
			}
			
			//update values
			OldPopID = NewPopID;
			NewPopID = "0";
		}
		else //you have encountered a new population
		{
			//add this population name to the AccessionNameList
			FullAccessionNameList.push_back (NewPopID);
			
			//process the list of alleles for the prior population, eliminate redundancies and missing data
			vector<std::string>().swap(TempList); //clear TempList
			for (q=0;q<ActiveAlleleList.size();q++)
			{
				ListToFilter = ActiveAlleleList[q];
				FilteredData = MyFilterDuplicates(ListToFilter);	
				TempList2d.push_back (FilteredData);
			}
			
			ActiveAlleleList = TempList2d;
			
			//put the allele list you've just created for the prior population into the final list
			ActiveAlleleByPopList.push_back (ActiveAlleleList);
			
			//zero out old variables so you can fill them again with the data from the new population
			vector<vector<std::string> >().swap(TempList2d); //clear TempList2d
			vector<std::string>().swap(ListToFilter); //clear ListToFilter
			vector<std::string>().swap(LocusNames); //clear LocusNames
			vector<vector<std::string> >().swap(ActiveAlleleList); //clear ActiveAlleleList
			OldPopID = NewPopID;
			
			//set up a switch that prevents acquiring a new line when you pop back up to the start of the 'while'
			IsNewPop = "yes";
		}
	}


	//process the list of alleles for the very last population, eliminate redundancies and missing data
	vector<std::string>().swap(TempList); //clear TempList
	for (q=0;q<ActiveAlleleList.size();q++)
	{
		ListToFilter = ActiveAlleleList[q];
		FilteredData = MyFilterDuplicates(ListToFilter);	
		TempList2d.push_back (FilteredData);
	}
			
	ActiveAlleleList = TempList2d;
			
	//add the final set of alleles for the very last population into the main list
	ActiveAlleleByPopList.push_back (ActiveAlleleList);

	//correct the length of the main list because the first item ends up empty
	ActiveAlleleByPopList.erase( ActiveAlleleByPopList.begin() );

	
	//ActiveAlleleByPopList and FullAccessionNameList references have now been completely updated
	return 0;
}


//reduces the master vector of all alleles into subsets containing reference or target loci only
int MyReduceToRef(vector<vector<vector<std::string> > > AllAlleleByPopList, vector<int> ReferenceOrTargetKey, vector<vector<vector<std::string> > >& ActiveAlleleByPopList, vector<vector<vector<std::string> > >& TargetAlleleByPopList)
{
	int r, t, i, j, k;
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
    int q=0;
    int CurrItemIndex;
    vector<vector<std::string> > ActiveAlleleList;
    vector<std::string> OldAlleles;
    vector<vector<std::string> > TempList2d;
    vector<std::string> FilteredData;
    vector<std::string> ListToFilter;
    std::string IsNewPop = "no";
    vector<std::string>::iterator it;

    int i,j,k,l;
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
	
	int i, j, k;
	vector<int> ActiveMaxAllelesList;
	vector<std::string> CurrLoc;
	set<std::string> NewSet;
	
//new way with set
	for (i=0;i<ActiveAlleleByPopList[0].size();++i)
	{
		NewSet.clear();
		for (j=0;j<ActiveAlleleByPopList.size();j++)
		{
			CurrLoc = ActiveAlleleByPopList[j][i]; //you are traversing the locus 'column' of the 3d grid
												   //get locus i for population j
			for (k=0;k<CurrLoc.size();++k)
			{
				NewSet.insert(CurrLoc[k]);	
			}
		}
		ActiveMaxAllelesList.push_back(NewSet.size()); //the set size after adding all populations for locus i is the maximum number of alleles
	}
	
	//for (i=0;i<ActiveMaxAllelesList.size();i++) cout << "ActiveMaxAllelesList[" << i << "]=" << ActiveMaxAllelesList[i] << "\n";
	
	return ActiveMaxAllelesList; 
}

int MyCalculateDiversity(vector<vector<vector<std::string> > > AlleleList, vector<int> ActiveMaxAllelesList, std::string Standardize, double& RandomActiveDiversity, double& AltRandomActiveDiversity)
{
	/*AlleleList structure:
		  Pop1..r
			  locusarray1..n		*/
	int CoreSize = AlleleList.size();
	int NumLoci = AlleleList[0].size();
	int i, j, k, M;
	vector<std::string> NewArray;
	vector<std::string> CurrLoc;
	vector<std::string> ListToFilter;
	vector<int> Mlist(NumLoci);
	set<std::string> AlleleSet;


	if (NumLoci == 0)
	{
		RandomActiveDiversity = -1; 
		AltRandomActiveDiversity = -1;
	}
	else
	{
		
		/*//slow using sort/erase
		for (i=0;i<NumLoci;i++)
		{
			//3. fuse alleles from the same locus into a single array, for all populations in core
			NewArray.clear();
			for (j=0;j<CoreSize;j++)
			{
				CurrLoc = AlleleList[j][i];
				NewArray.insert(NewArray.end(), CurrLoc.begin(), CurrLoc.end());
			}
			
			for (k=0;k<NewArray.size();k++)
			{
				cout << "NewArray["<<k<<"]="<<NewArray[k]<<"\n";
			}
			
			//4. assemble a list of diversity (M) for each locus separately
			//ListToFilter starts with a sorted list of all unique alleles found at a locus.  If there are no alleles (i.e. all missing data), M is set to 0.  Otherwise, M is set to the number of unique alleles at the locus.
			ListToFilter.clear();
			//remove duplicates
			ListToFilter = NewArray;
			sort( ListToFilter.begin(), ListToFilter.end() );
			ListToFilter.erase( std::unique( ListToFilter.begin(), ListToFilter.end() ), ListToFilter.end() );

			
			cout << "Locus "<<i<<"\n";
			for (k=0;k<ListToFilter.size();k++)
			{
				cout << "ListToFilter["<<k<<"]="<<ListToFilter[k]<<"\n";
			}

			if (ListToFilter.size() == 0)
			{
				M=0;
			}
			else
			{
				M=ListToFilter.size(); //M=the number of alleles present
			}
			
			//Mlist contains number of alleles present at each locus
			Mlist[i] = M;
		}
		*/
		
		//fast using set
		for (i=0;i<NumLoci;i++)
		{

			//3. pass alleles from the same locus into a single set, for all populations in core, to remove redundancies
			AlleleSet.clear(); //clear AlleleSet
			for (j=0;j<CoreSize;j++)
			{
				CurrLoc = AlleleList[j][i];
				for (k=0;k<CurrLoc.size();++k)
				{
					AlleleSet.insert(CurrLoc[k]); //locus i for all population j
				}
			}
			
			if (AlleleSet.size() == 0) M=0;
			else M=AlleleSet.size();
			
			Mlist[i] = M;
		}

		
		/*for (k=0;k<Mlist.size();++k)
		{
			cout << Mlist[k] << ",";
		}
		cout << "\n";
		*/
		
		
		
		//5. standardize the M values to the maximum possible number of alleles at that locus, 
		//and add them up to get final estimate of standardized allelic diversity in the core.
		//Then divide by the number of loci to get a number that is comparable across data sets.
		
		//calculate the standardized allelic diversity
		double SAD;
		double SADtemp = 0; //SADtemp is the summed value of standardized M across all loci for a single subcore
		for (i=0;i<NumLoci;i++)
		{
			SADtemp = SADtemp + ( (double) Mlist[i] / (double) ActiveMaxAllelesList[i] );
		}
		SAD = SADtemp / NumLoci;
		
		
		//calculate M by simply adding up the Mlist
		//M = std::accumulate(Mlist.begin(),Mlist.end(),0);  <--std::accumulate is unreliable so it has been removed
		M = 0;
		for (i=0;i<Mlist.size();++i)
		{
			M = M + Mlist[i];
		}
		
		//6. Determine the way variables are updated so that the value for the desired optimality 
		//criterion is RandomActiveDiversity, while the other is AltRandomActiveDiversity.  
		//This way, the output can contain information on both M+ and M criteria, 
		//although only one will be used for optimization.
		if (Standardize == "yes")
		{
			RandomActiveDiversity = SAD;
			AltRandomActiveDiversity = M;
		}
		else if (Standardize == "no")
		{
			RandomActiveDiversity = M;
			AltRandomActiveDiversity = SAD;
		}
	}
	return 0;
}

//write current results for each thread to a recovery file
void WriteRecoveryFile (const char* RecoveryFilePath, int r, double StartingRandomActiveDiversity, double best, double RandomTargetDiversity, double OptimizedTargetDiversity, double StartingAltRandomActiveDiversity, double AltOptimizedActiveDiversity, double AltRandomTargetDiversity, double AltOptimizedTargetDiversity, vector<std::string> TempListStr)
{
	stringstream ss;
	ss << omp_get_thread_num();
	string foo = ss.str();
	ofstream RecoveryFile;
	RecoveryFile.open(RecoveryFilePath, ios::out | ios::app); //open file in append mode
	RecoveryFile <<	r << "\t" << StartingRandomActiveDiversity << "\t" << best << "\t" << RandomTargetDiversity << "\t" << OptimizedTargetDiversity << "\t" << StartingAltRandomActiveDiversity << "\t" << AltOptimizedActiveDiversity << "\t" << AltRandomTargetDiversity << "\t" << AltOptimizedTargetDiversity << "\t(";
	for (int i=0;i<TempListStr.size();++i)
	{
		if (i == (TempListStr.size() - 1) )
		{
			RecoveryFile << TempListStr[i] << ")\n";
		}
		else RecoveryFile << TempListStr[i] << ",";	
	}
	RecoveryFile.close();
}

void printProgBar( int percent)
{
  std::string bar;

  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
}

  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "% complete     thread "<<omp_get_thread_num() <<" "<< std::flush;
}

bool fileExists(const char *fileName)
{
    ifstream infile(fileName);
    return infile.good();
}

vector<std::string> MyRemoveTargetAlleles(vector<std::string> AllAlleles, vector<int> AllColumnIDList, vector<int> TargetColumnIDList)
{
	int checker, j;
	int i=0;
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
	int j; //counts items per row
	int i=0, k=0; //k is the row number
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
	
	int i, j, locindex, k;
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
	int i, j, NumAlleles, z;
	double freq;
	std::string b;
	vector<std::string> AllAlleles;
	vector<std::string> UniqAlleles;
	set<std::string> AlleleSet;

	Alfreq laf; //locusname, allelenames, frequencies
	//vector<vector<double> > freqvec(RefAllelesByLocus.size()); //to set size of struct frequencies
	//laf.frequencies = freqvec;
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
		//cout << laf.locusname << "\n";
		
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
			
			//cout << " b="<< b<<" z="<<z<<"\n";
			//cout << " laf.frequencies["<<j<<"]="<<laf.frequencies[j]<<"\n";
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
	const int MinCoreSize = atoi(argv[3]);
	const int MaxCoreSize = atoi(argv[4]);
	const int SamplingFreq = atoi(argv[5]);
	const int NumReplicates = atoi(argv[6]);
	char* OutFilePath = argv[7];
	char* SumFilePath = NULL;
	char* KerFilePath = NULL;
	char* IdealFilePath = NULL;
	
	//declare variables
	int i, j, k, l, b;
	string DoSummary = "no"; //switch to turn on summary function
	string Kernel = "no"; //switch to include a mandatory set in the core
	string Ideal = "no"; //switch to compute the ideal core, using A* algorithm
	vector<std::string> KernelAccessionList;
	vector<std::string> BadFiles;
	string bf;

	//parse the command line for options
	for (i=0;i<argc;i++)
	{
		if ( string(argv[i]) == "-s" ) 
    	{
        	DoSummary = "yes";
        	SumFilePath = argv[i+1];
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
		
		if ( string(argv[i]) == "-b" ) 
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
	cout << "Input variables:\n  VarFilePath = " << VarFilePath << "\n";
	cout << "  DatFilePath = " << DatFilePath << "\n";
	if (KerFilePath != NULL) cout << "  KerFilePath = " << KerFilePath << "\n";
	if (SumFilePath != NULL) cout << "  SumFilePath = " << SumFilePath << "\n";
	cout << "  MinCoreSize = " << MinCoreSize << "\n";
	cout << "  MaxCoreSize = " << MaxCoreSize << "\n";
	cout << "  SamplingFreq = " << SamplingFreq << "\n";
	cout << "  NumReplicates = " << NumReplicates << "\n";
	cout << "  OutFilePath = " << OutFilePath << "\n";
	
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
    if (KernelAccessionList.size() > MinCoreSize)
    {
		cout << "ERROR:  The number of mandatory accessions ("<<KernelAccessionList.size()
			 <<") is greater than the minimum core size ("<<MinCoreSize
        	 <<").  Please modify the kernel file or the command line argument.  Quitting...\n\n";
				exit (EXIT_FAILURE);
	}

	

	
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
			cout << "\n\nActive Locus	Column\n";
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
	
	//.dat file, all loci
	
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
	
	
/*	//switch for old MyProcessDatFile
	vector<vector<vector<std::string> > > AllAlleleByPopList; //structure of this 3D vector is:
	// { { {pop1,loc1 alleles},{pop1,loc2 alleles},...}, { {pop2,loc1 alleles},{pop2,loc2 alleles},...} } }
	MyProcessDatFile(DatFilePath, AllColumnIDList, AllLociNameList, AllAlleleByPopList, FullAccessionNameList, IndivPerPop, AllAlleles);
	//AllAlleleByPopList, FullAccessionNameList, IndivPerPop, AllAlleles updated as reference
*/	
		
		/*//Print out FullAccessionNameList
		cout << "\n\nPopulation names\n";
		for (i=0; i<FullAccessionNameList.size();i++) cout << FullAccessionNameList[i] << "\n";
		
		//Print out lists of unique All alleles from AllAlleleByPopList
		for (i=0;i<AllAlleleByPopList.size() ;i++)
		{
			cout << "Population " << FullAccessionNameList[i] << "\n";
			for (j=0;j<AllAlleleByPopList[i].size();j++)
			{
				cout << "Locus " << j << "\n";
				for (k=0;k<AllAlleleByPopList[i][j].size();k++)
				{
					cout << AllAlleleByPopList[i][j][k] << ",";
				}
				cout << "\n";
			}
		
		}*/
	
	
	
	
	
	
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
	int NumberOfAccessions = ActiveAlleleByPopList.size();
	
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
	




    /*cout << "ActiveMaxAllelesList:\n";
	for (i=0; i<ActiveMaxAllelesList.size(); i++) cout << ActiveMaxAllelesList[i] << ",\n";
    cout << "\nTargetMaxAllelesList:\n";
	for (i=0; i<TargetMaxAllelesList.size(); i++) cout << TargetMaxAllelesList[i] << ",\n";
	*/
	
	
	//PERFORM A*
	if (Ideal == "yes")
	{
		int parallelism_enabled = 1; //0=no, not 0 = yes
		if (parallelism_enabled == 0) cout << "\nBeginning serial A* search...\n\n";
		else cout << "\nBeginning parallel A* search...\n\n";
				
		//start the clock
		time_t start1,end1;
		time (&start1);

	aStar(IdealFilePath, ActiveAlleleByPopList, ActiveMaxAllelesList, UniqLociNamesList, ReferenceOrTargetKey, FullAccessionNameList, PloidyList, PopSizes, AlleleFrequencies, parallelism_enabled);
		
		//stop the clock
		time (&end1);
		double dif = difftime (end1,start1);
		cout << "\nA* search complete.  Elapsed time = "<< dif << " seconds.\n\n";
	}	
	
	
	
	//PERFORM M+
	
	//start the clock
	time_t startm,endm;
	time (&startm);

	
	//set up variables for monitoring progress
	int percent; //percent of analysis completed
	int progindex = 0;  //index to monitor progress, percent = 100*(progindex/l)
	//below is a stupid way to calculate the number of rows in the output file, value l (which = V1) 
	//used to monitor progress and as the maximum vector index for shared output vectors
	l=0;
	for (i=MinCoreSize;i<MaxCoreSize+1;i=i+SamplingFreq)
	{
		for (int j=0;j<NumReplicates;j++)
		{
			l++;
		}
	}
		

	//set up vectors to fill with results
	double V1 = l; //(MaxCoreSize - MinCoreSize + 1)*NumReplicates; //number of rows in output vectors
	vector<vector<double> > Results(V1, vector<double>(9)); //will contain numerical results
	vector<vector<string> > Members(V1); //will contain core set members
		
	
	


	
	//compile as parallel or not?
	int parallelism_enabled = 1; //0=no, not 0 = yes
	if (parallelism_enabled == 0) cout << "\nBeginning serial M+ search...\n\n";
	else cout << "\nBeginning parallel M+ search...\n\n";

	#pragma omp parallel if(parallelism_enabled) 
	{		
		int r, nr, RandAcc, b, row, bsc, plateau; //r = core size, nr = controller to repeat NumReplicates times
									//row = result vector row number, bsc = holds best sub core member, and other indexed accessions
									//plateau = index of the number of reps in optimization loop with same diversity value
									
		double RandomActiveDiversity;
		double AltRandomActiveDiversity;
		double StartingRandomActiveDiversity;
		double StartingAltRandomActiveDiversity;
		double RandomTargetDiversity;
		double AltRandomTargetDiversity;
		double StartingDiversity;
		double TempAltOptimizedActiveDiversity;
		double AltOptimizedActiveDiversity;
		double OptimizedTargetDiversity;
		double AltOptimizedTargetDiversity;
		double best;
		double nnew;
		vector<vector<vector<std::string> > > AlleleList;
		vector<vector<vector<std::string> > > CoreAlleles;
		vector<vector<vector<std::string> > > TdTempList;
		vector<vector<vector<std::string> > > BestSubCoreAlleles;
		std::string Standardize = "yes";  //a run that mimics the MSTRAT approach can be accomplished by setting Standardize="no", and setting up the var file so that each column in the .dat file is treated as a single locus, rather than two (or more) adjacent columns being treated as a single codominant locus.
		vector<int> AccessionsInCore;
		vector<int> AccessionsInSubCore;
		vector<int> BestSubCore;
		vector<int> BestSubCoreRevSorted;
		vector<int> TempList;
		vector<int> TempList2;
		vector<int> bestcore;
		vector<std::string> TempListStr;
	
		//seed the random number generator for each thread
		int tt;
		tt = (time(NULL));
		srand ( tt ^ omp_get_thread_num() ); //initialize
		
		//set up a recovery file for each thread that saves progress as program runs
		stringstream ss;
		ss << OutFilePath << ".t" << omp_get_thread_num() << ".tmp"; 
		string rfp = ss.str();
		const char* RecoveryFilePath = rfp.c_str();
		ofstream RecoveryFile; 
		RecoveryFile.open(RecoveryFilePath);
		RecoveryFile.close(); //quick open close done to clear any existing file each time program is run
		RecoveryFile.open(RecoveryFilePath, ios::out | ios::app); //open file in append mode
		RecoveryFile << "core size	random active diversity	optimized active diversity	random target diversity	optimized target diversity	alt random active diversity	alt optimized active diversity	alt random target diversity	alt optimized target diversity	core members" << "\n";
		RecoveryFile.close();



		//do parallelization so that each rep by core size combo can be
		//handled by a distinct thread.  this involves figuring out the total
		//number of reps*coresizes taking into account the SamplingFreq

		int rsteps = 1 + floor( (MaxCoreSize - MinCoreSize) / SamplingFreq ); //number of steps from MinCoreSize to MaxCoreSize

		#pragma omp for
			for (int rnr = 0; rnr<rsteps*NumReplicates;++rnr)
			{
				r = MinCoreSize + ((rnr / NumReplicates) * SamplingFreq); //int rounds to floor
				nr = rnr % NumReplicates; // modulo
				
				

		/*old parallelization, outer loop only, if reinstated, make sure closing braces are corrected too
		#pragma omp for
			for (int r=MinCoreSize;r<MaxCoreSize+1;r=r+SamplingFreq)
			{
				//cout << "  CoreSize = " << r << "\n";

				for (int nr=0;nr<NumReplicates;nr++)
				{
		*/		


					//cout << "    rep " << nr + 1 << "\n";
					//develop random starting core set
					//clear AccessionsInCore and set size
					AccessionsInCore.clear();
					AccessionsInCore.resize(r);
					
					//add kernel accessions to core, if necessary
					if (Kernel == "yes")
					{
						for (int i=0;i<KernelAccessionIndex.size();i++)
						{
							AccessionsInCore[i] = KernelAccessionIndex[i];
						}
					}

					//clear TempList and set size					
					TempList.clear();
					TempList.resize( AccessionNameList.size() );
					//vector<int>().swap(AccessionsInCore); //clear AccessionsInCore
					//vector<int>().swap(TempList); //clear TempList
					
					//set list of available accessions in TempList, by erasing those already in the core
					TempList = AccessionNameList;
					//expunge the kernel accessions, so they are not available for random addition below
					//KernelAccessionIndex has been reverse sorted so you don't go outside range after automatic resize by .erase
					for (int i=0;i<KernelAccessionIndex.size();i++)
					{
						b = KernelAccessionIndex[i];
						TempList.erase(TempList.begin()+b);
					}
				
					//randomly add accessions until r accessions are in the core. if there is a kernel, include those (done above)
					//plus additional, randomly selected accessions, until you get r accessions
					//for (int i=0;i<r;i++)
					for (int i=KernelAccessionIndex.size();i<r;i++)
					{
						//choose an accession randomly from those available
						RandAcc = rand() % TempList.size();
						//add it to the list
						AccessionsInCore[i] = TempList[RandAcc];
						//AccessionsInCore.push_back(TempList[RandAcc]);
						
						//remove it from the list of available accessions
						TempList.erase(TempList.begin()+RandAcc);
					}
			
					
					//assemble genotypes for random core and calculate diversity
					//1. put together initial list of active alleles
					//vector<vector<vector<std::string> > >().swap(CoreAlleles); //clear CoreAlleles
					CoreAlleles.clear();
					CoreAlleles.resize( AccessionsInCore.size() );
					for (int i=0;i<AccessionsInCore.size();i++)
					{
						b = AccessionsInCore[i];
						//CoreAlleles.push_back(ActiveAlleleByPopList[b]);
						CoreAlleles[i] = ActiveAlleleByPopList[b];
					}

					//2. calculate diversity from random selection at active loci
					//vector<vector<vector<std::string> > >().swap(AlleleList); //clear AlleleList
					AlleleList.clear();
					AlleleList = CoreAlleles;
			
					MyCalculateDiversity(AlleleList, ActiveMaxAllelesList, Standardize, RandomActiveDiversity, AltRandomActiveDiversity);
					//in MyCalculateDiversity, latter two variables are updated as references
					//save them away in non-updated variables
					StartingRandomActiveDiversity = RandomActiveDiversity;
					StartingAltRandomActiveDiversity = AltRandomActiveDiversity;
					
		//cout << "RandomActiveDiversity="<<RandomActiveDiversity<<"\n";
		//cout << "AltRandomActiveDiversity="<<AltRandomActiveDiversity<<"\n";


					//3. calculate diversity from random selection at target loci
					//vector<vector<vector<std::string> > >().swap(AlleleList); //clear AlleleList
					AlleleList.clear();
					AlleleList.resize( AccessionsInCore.size() );
					for (int j=0;j<AccessionsInCore.size();j++)
					{
						b = AccessionsInCore[j];
						//AlleleList.push_back(TargetAlleleByPopList[b]);
						AlleleList[j] = TargetAlleleByPopList[b];
					}
					MyCalculateDiversity(AlleleList, TargetMaxAllelesList, Standardize, RandomTargetDiversity, AltRandomTargetDiversity);

		
		
		//cout << "RandomTargetDiversity="<<RandomTargetDiversity<<"\n";
		//cout << "AltRandomTargetDiversity="<<AltRandomTargetDiversity<<"\n";
		

		
					//BEGIN OPTIMIZATION
					StartingDiversity = 0; //this is the diversity recovered during the prior iteration.
					plateau = 0; //count of the number of times you have found the best value, evaluates when you are
								 //stuck on a plateau, assuming acceptance criterion allows downhill steps
					//this is the iterations step, now an indefinite loop that is broken when 
					//no improvement is made during the course of the optimization algorithm
					//If r = kernel size = MinCoreSize then do no optimization but still calculate all variables.
					if (KernelAccessionIndex.size() == r)
					{
							//assemble genotypes for core
							//1. put together initial list
							CoreAlleles.clear();
							CoreAlleles.resize(r);
							for (int i=0;i<r;i++)
							{
								b = AccessionsInCore[i];
								CoreAlleles[i] = ActiveAlleleByPopList[b];
							}
							
							AlleleList = CoreAlleles;
							
							MyCalculateDiversity(AlleleList, ActiveMaxAllelesList, Standardize, RandomActiveDiversity, AltRandomActiveDiversity);
							best = RandomActiveDiversity; //best is equivalent to OptimizedActiveDiversity
							AltOptimizedActiveDiversity = AltRandomActiveDiversity;
					}
					else
					{
						//do optimization
						while ( true )
						{
							//assemble genotypes for core
							//1. put together initial list
							CoreAlleles.clear();
							CoreAlleles.resize(r);
							for (int i=0;i<r;i++)
							{
								b = AccessionsInCore[i];
								CoreAlleles[i] = ActiveAlleleByPopList[b];
							}
				
							//2. go through all possible subsets of size r-1, one at a time, noting which is best.
							//If there is a kernel, do not swap out any of those accessions (they are retained as the
							//first KernelAccessionIndex.size() items in CoreAlleles).  Accomplished by starting for loop
							//at KernelAccessionIndex.size().
							best=0;
							//for (int i=0;i<CoreAlleles.size();i++)
							for (int i=KernelAccessionIndex.size();i<CoreAlleles.size();i++)
							{
								//remove each item consecutively from the list of all populations in the core
								AlleleList.clear();
								TdTempList.clear();
							
								TdTempList = CoreAlleles; //swap to temporary vector
								TdTempList.erase( TdTempList.begin() + i);
								AlleleList = TdTempList;
					
								TempList2.clear();
								TempList2 = AccessionsInCore;
								TempList2.erase(TempList2.begin() + i);
								AccessionsInSubCore = TempList2;


								/*Data structure for SubCoreAlleles:
								SubCore 1..r
									Population 1..(r-1)
										AlleleArray 1..NumLoci		
					
								--3. fuse alleles from the same locus into a single array, for all accessions, for the current subcore
								--4. assemble a list of diversity (M) for each locus separately
								--5. standardize the M values to the maximum possible number of alleles at that locus, and add them up to get final estimate of standardized allelic diversity in the core.  then divide by the number of loci to get a number that is comparable across data sets.
								--5.5. simultaneous to the calculation, keep track of which subcore is best
								*/
					
								MyCalculateDiversity(AlleleList, ActiveMaxAllelesList, Standardize, RandomActiveDiversity, AltRandomActiveDiversity);
								nnew = RandomActiveDiversity;
		
								if (nnew >= best) // >= allows sideways movement during hill climbing
								{
									best = nnew;

									BestSubCore.clear();
									BestSubCore = AccessionsInSubCore;
									BestSubCoreAlleles.clear();
									BestSubCoreAlleles = AlleleList;
								}
							}  //for loop cycles thru all subcores

							//reverse sort BestSubCore to support easy assembly of pared TempList below
							BestSubCoreRevSorted = BestSubCore;
							std::sort(BestSubCoreRevSorted.begin(), BestSubCoreRevSorted.end(), std::greater<int>());
			
							/*
							6. take the subcore with greatest diversity and consecutively add each 
							possible additional accession from the base collection.  find the core of size r 
							(not r-1 subcore) that has the greatest diversity.

							suppress the IDs of those accessions found in the BestSubCore from the 
							list of all accessions to get a list of remaining accessions.*/
							TempList = AccessionNameList;
							for (int k=0;k<BestSubCoreRevSorted.size();k++)
							{
								bsc = BestSubCoreRevSorted[k];
								TempList.erase( TempList.begin() + bsc );
							}
						
							//shuffle the list of remaining accessions, so addition order is not predictable
							std::random_shuffle (TempList.begin(), TempList.end());
							
							//add each remaining accession consecutively, calculate diversity, test 
							//whether it is better than the prior one
							best = 0;
							for (int k=0;k<TempList.size();k++)
							{
								bsc = TempList[k];
					
								//define the core
								TempList2 = BestSubCore;
								TempList2.resize( TempList2.size() + 1 );
								//TempList2.push_back(i);
								TempList2[TempList2.size()-1] = bsc; //add new accession to last vector element
								AccessionsInCore = TempList2;
					
								//assemble the allelelist for the core
								TdTempList = BestSubCoreAlleles;
								TdTempList.resize( TdTempList.size() + 1 );
								//TdTempList.push_back( ActiveAlleleByPopList[i] );
								TdTempList[TdTempList.size()-1] = ActiveAlleleByPopList[bsc];
								AlleleList = TdTempList;
					
								//calculate diversity
								MyCalculateDiversity(AlleleList, ActiveMaxAllelesList, Standardize, nnew, TempAltOptimizedActiveDiversity); 


//cout << nnew << "\n";

				
								//test whether current diversity is higher than the best diversity found so far
								if (nnew >= best) // >= allows sideways movement during hill climbing
								{
									best = nnew;
									bestcore = AccessionsInCore;
									//save the alternative diversity value for the best core
									AltOptimizedActiveDiversity = TempAltOptimizedActiveDiversity;
								}
				
							}

							AccessionsInCore = bestcore; //define starting variable for next MSTRAT iteration
				
							//if there has been no improvement from the prior iteration, you have reached
							// the plateau and should exit the repeat
							if (best == StartingDiversity) 
							{
								plateau++;
								if (plateau > 0) break;
							}
							//update starting value and repeat
							else if (best > StartingDiversity) StartingDiversity = best;
							
							/*//standard break
							if (best == StartingDiversity) break;
							else if (best > StartingDiversity) StartingDiversity = best;
							*/
								
						} //while(true) endless loop
					}
		
					//7. Calculate diversity at target loci
					//assemble the target loci allelelist for the accessions in the best core
					//vector<vector<vector<std::string> > >().swap(AlleleList); //clear AlleleList
					AlleleList.clear();
					AlleleList.resize( AccessionsInCore.size() );
					for (int j=0;j<AccessionsInCore.size();j++)
					{
						b = AccessionsInCore[j];
						//AlleleList.push_back(TargetAlleleByPopList[b]);
						AlleleList[j] = TargetAlleleByPopList[b];
					}
			
			
					//calculate diversity at target loci based upon the optimized core selection
					MyCalculateDiversity(AlleleList, TargetMaxAllelesList, Standardize, OptimizedTargetDiversity, AltOptimizedTargetDiversity);

			
					//8. Assemble stats for optimized core and add to output vectors
					//create a list of accession names from the list of accession ID's in AccessionsInCore
					sort( AccessionsInCore.begin(), AccessionsInCore.end() );
					
					TempListStr.clear();
					TempListStr.resize(r);
					for (int i=0;i<AccessionsInCore.size();i++)
					{
						b = AccessionsInCore[i];
						TempListStr[i] = FullAccessionNameList[b];
					}

					//load the variables onto the results vectors
					//numerical results
					row = ((r - MinCoreSize)*NumReplicates) + nr - ( (NumReplicates*(SamplingFreq-1))*( (r-MinCoreSize)/SamplingFreq ) );
					// (r - MinCoreSize)*NumReplicates) + nr specifies row number if SamplingFreq=1
					// (NumReplicates*(SamplingFreq-1)) specifies a step value to correct when SamplingFreq>1
					// ( (r-MinCoreSize)/SamplingFreq ) specifies the replicate on core size, accounting for SamplingFreq
					// see file Calculation of row value.xlsx for development of the 'row' index
					Results[row][0] = r;
					Results[row][1] = StartingRandomActiveDiversity;//RandomActiveDiversity;
					Results[row][2] = best; //equivalent to OptimizedActiveDiversity
					Results[row][3] = RandomTargetDiversity;
					Results[row][4] = OptimizedTargetDiversity;
					Results[row][5] = StartingAltRandomActiveDiversity;//AltRandomActiveDiversity;
					Results[row][6] = AltOptimizedActiveDiversity;
					Results[row][7] = AltRandomTargetDiversity;
					Results[row][8] = AltOptimizedTargetDiversity;
			
					//core set members
					Members[row] = TempListStr; 
			
					//write the results onto the recovery files
					WriteRecoveryFile(RecoveryFilePath, r, StartingRandomActiveDiversity, best, RandomTargetDiversity, OptimizedTargetDiversity, StartingAltRandomActiveDiversity, AltOptimizedActiveDiversity, AltRandomTargetDiversity, AltOptimizedTargetDiversity, TempListStr);
					
					//display progress
					progindex = progindex + 1;
					percent = 100*(progindex/V1);
					printProgBar(percent); 
					//cout << progindex << "\n";
					
				//} //end NumReplicates
			//} //end   for (r=MinCoreSize;r<MaxCoreSize+1;r=r+SamplingFreq)
			} //end   for (int rnr = 0; rnr<rsteps*NumReplicates;++rnr)
					
	} //end #pragma omp parallel	
	

	//set up file stream for output file
	ofstream output; 
	output.open(OutFilePath);
	output.close(); //quick open close done to clear any existing file each time program is run
	output.open(OutFilePath, ios::out | ios::app); //open file in append mode
	output << "core size	random active diversity	optimized active diversity	random target diversity	optimized target diversity	alt random active diversity	alt optimized active diversity	alt random target diversity	alt optimized target diversity	core members" << "\n";
			
	//write out results row by row
	for (i=0;i<V1;i++)
	{
		//write variables
		output 	<< Results[i][0] 
				<< "	" << Results[i][1] 
				<< "	" << Results[i][2] 
				<< "	" << Results[i][3] 
				<< "	" << Results[i][4] 
				<< "	" << Results[i][5] 
				<< "	" << Results[i][6] 
				<< "	" << Results[i][7] 
				<< "	" << Results[i][8] 
				<< "	" << "(";
		//write Accessions retained
		for (j=0;j<Members[i].size();j++)
		{
			if ( j==(Members[i].size() - 1) )
			{
				//add trailing parentheses and move to next row
				output << Members[i][j] << ")\n";
			}
			else
			{
				output << Members[i][j] << ",";
			}
		}
	}
	
	//wrap up write step
	output.close();

	//delete all recovery files
	#pragma omp parallel if(parallelism_enabled) 
	{		
		stringstream ss;
		ss << OutFilePath << ".t" << omp_get_thread_num() << ".tmp"; 
		string rfp = ss.str();
		const char* RecoveryFilePath = rfp.c_str();
		remove(RecoveryFilePath);
	}
	
	
	//stop the clock
	time (&endm);
	dif = difftime (endm,startm);
	cout << "\nM+ search complete.  Elapsed time = "<< dif << " seconds.\n\n";





}
