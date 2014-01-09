#include <sstream>
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
#include "m+.hpp"
using namespace std;


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
-k kernelfile = use an MSTRAT .ker file to specify mandatory members of the 
          core.  The number of included accessions must be less than or equal to mincoresize.

example: ./m+ ./beet.var ./beet.dat 3 28 2 3 ./beetout.txt -s ./beetsum.txt -k beet.ker
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


int MyProcessVarFile(char* VarFilePath, vector<int>& AllColumnIDList, vector<std::string>& AllLociNameList, vector<int>& ActiveColumnIDList, vector<std::string>& ActiveLociNameList, vector<int>& TargetColumnIDList, vector<std::string>& TargetLociNameList, vector<vector<int> >& ColKeyToAllAllelesByPopList, vector<int>& ReferenceOrTargetKey)
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
	vector<std::string> UniqLociNameList;
	
	//find unique locus names
	//easy way, but messes up sort order
	/*
	UniqLociNameList = AllLociNameList;
	sort( UniqLociNameList.begin(), UniqLociNameList.end() );
	UniqLociNameList.erase( std::unique( UniqLociNameList.begin(), UniqLociNameList.end() ), UniqLociNameList.end() );
	*/
	//hard way, retains sort order
	UniqLociNameList = unsortedRemoveDuplicates(AllLociNameList);
	
	//define 2d vector that is key to columns, size to number of unique loci
	ColKeyToAllAllelesByPopList.resize( UniqLociNameList.size() ); //size to number of loci
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
				ColKeyToAllAllelesByPopList[i].push_back(b);	
			}	
		}		
	}
	
	//define a 1d vector that describes whether the loci in the ColKey are reference(0) or target(1)
	double rt;
	ReferenceOrTargetKey.resize( ColKeyToAllAllelesByPopList.size() ); //sized to same length as key
	for (i=0;i<ColKeyToAllAllelesByPopList.size();++i)
	{
		rt=0;
		//test whether all elements are categorized as reference or as target, if not, raise error
		for (k=0;k<ColKeyToAllAllelesByPopList[i].size();++k)
		{
			b=ColKeyToAllAllelesByPopList[i][k];
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
		else if (  rt/ColKeyToAllAllelesByPopList[i].size() == 1 ) ReferenceOrTargetKey[i] = 1; //it is a target
		else 
		{
			cout << "ERROR:  Some loci are described as both reference and target.  Please check your var file. Quitting...\n\n";
			exit (EXIT_FAILURE);
		}
		
	}
	/*
	//print out ColKeyToAllAllelesByPopList
	for (i=0;i<ColKeyToAllAllelesByPopList.size();++i)
	{
		cout << "Locus "<<UniqLociNameList[i]<< " {";
		for (k=0;k<ColKeyToAllAllelesByPopList[i].size();++k)
		{
			cout<<ColKeyToAllAllelesByPopList[i][k]<<",";
		}
		cout << "\n";
	}
	
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

int MyProcessDatFile(char* DatFilePath, vector<int> ActiveColumnIDList, vector<std::string> ActiveLociNameList, vector<vector<vector<std::string> > >& ActiveAlleleByPopList, vector<std::string>& FullAccessionNameList)
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


	//read in the dat file line by line	
	std::ifstream infile;
	infile.open(DatFilePath);
    while( !infile.eof() ) // To get all the lines.
    {
	    if (IsNewPop == "no") //get a new line
	    {
			std::getline (infile, foo); // Saves the line in foo.

			//split on whitespace
			foovector = split(foo);

			//extract population identifier from first column
			NewPopID = foovector[0];
			
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
								
				//get the name of the current locus
				CurrLocusName = ActiveLociNameList[i];
				
				/*
				cout << "ActiveColumnID=" << ActiveColumnIDList[i] << "\t";
				cout << "NewAllele=" << NewAllele << "\t";
				cout << "CurrLocusName=" << CurrLocusName << "\n";
				*/
				
				//test whether you already have a list item for this locus, 
				//note this works for both building the initial list of arrays (using the first part of "or") or adding to the existing list of arrays (after "or")
				//count command tests whether CurrLocusName appears in LocusNames
				if (CurrLocusName == OldLocusName || std::count (LocusNames.begin(), LocusNames.end(), CurrLocusName) >= 1 )
				{
					//determine whether you are building the list or adding to it, then figure out what gene you are in
					if (std::count (LocusNames.begin(), LocusNames.end(), CurrLocusName) >= 1)
					{
						for (j=0; j<LocusNames.size(); j++)
						{
							if (LocusNames[j] == CurrLocusName)
							{
								CurrItemIndex = j;
								break;	
							}
						}
					
					}
					
					//add to the current array being built or to the appropriate existing array
					OldAlleles = ActiveAlleleList[CurrItemIndex]; //pull out the list of alleles currently held for this locus
					OldAlleles.push_back (NewAllele); //add the new allele to the existing list of alleles for the current locus
					ActiveAlleleList[CurrItemIndex] = OldAlleles; //update ActiveAlleleList
					
					/*cout << "OldAlleles\n";
					for (j=0;j<OldAlleles.size();j++)
					{
						cout << OldAlleles[j] << ",";
					}
					cout << "\n";
					*/
				}
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
	
	//resize level 2 of vectors, they contain the number of loci in ReferenceOrTargetKey
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


//MyProcessDatFile using fread and parallelization--in dev--this can't be done because current line
//is categorized contingent on previous line: is it a new population?
int MyProcessDatFileII(char* DatFilePath, vector<int> ActiveColumnIDList, vector<std::string> ActiveLociNameList, vector<vector<vector<std::string> > >& ActiveAlleleByPopList, vector<std::string>& FullAccessionNameList)
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

/*
	//read the whole file into a buffer using fread
    char * buffer;
	cout << "gettingfile" << "\n";
	buffer = MyBigRead(DatFilePath);
	
	cout << "gotfile" << "\n";
*/

	//read in the dat file line by line	
	// code using file stream
	std::ifstream infile;
	infile.open(DatFilePath);
	
	int zz = 0;
	
    while( !infile.eof() ) // To get all the lines.
    {
	    if (IsNewPop == "no") //get a new line
	    {
			std::getline (infile, foo); // Saves the line in foo.

			//split on whitespace
			foovector = split(foo);
			
			//extract population identifier from first column
			NewPopID = foovector[0];
			
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
								
				//get the name of the current locus
				CurrLocusName = ActiveLociNameList[i];

				/*
				cout << "ActiveColumnID=" << ActiveColumnIDList[i] << "\t";
				cout << "NewAllele=" << NewAllele << "\t";
				cout << "CurrLocusName=" << CurrLocusName << "\n";
				*/
				
				//test whether you already have a list item for this locus, note this works for both building the initial list of arrays (using the first part of "or") or adding to the existing list of arrays (after "or")
				//count command tests whether CurrLocusName appears in LocusNames
				if (CurrLocusName == OldLocusName || std::count (LocusNames.begin(), LocusNames.end(), CurrLocusName) >= 1 )
				{
					//determine whether you are building the list or adding to it, then figure out what gene you are in
					if (std::count (LocusNames.begin(), LocusNames.end(), CurrLocusName) >= 1)
					{
						for (j=0; j<LocusNames.size(); j++)
						{
							if (LocusNames[j] == CurrLocusName)
							{
								CurrItemIndex = j;
								break;	
							}
						}
					
					}
					
					//add to the current array being built or to the appropriate existing array
					OldAlleles = ActiveAlleleList[CurrItemIndex]; //pull out the list of alleles currently held for this locus
					OldAlleles.push_back (NewAllele); //add the new allele to the existing list of alleles for the current locus
					ActiveAlleleList[CurrItemIndex] = OldAlleles; //update ActiveAlleleList
					
					/*cout << "OldAlleles\n";
					for (j=0;j<OldAlleles.size();j++)
					{
						cout << OldAlleles[j] << ",";
					}
					cout << "\n";
					*/
				}
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
	int i, j;
	vector<vector<std::string> > TempList;
	vector<vector<std::string> > MaxAlleles;
	vector<int> ActiveMaxAllelesList;
	vector<std::string> NewArray;
	vector<std::string> CurrLoc;
	vector<std::string> ListToFilter;
	vector<std::string> TempList2;
	int Mtemp;
	
	//1. fuse alleles from the same locus into a single array, for all accessions
	for (i=0;i<ActiveAlleleByPopList[1].size();i++) //item 1 is as good as any for defining the number of loci
	{
		vector<std::string>().swap(NewArray); //clear NewArray
		for (j=0;j<ActiveAlleleByPopList.size();j++)
		{
			CurrLoc = ActiveAlleleByPopList[j][i]; //you are traversing the locus 'column' of the 3d grid
			//cout << "j=" << j << ", CurrLoc.size()=" << CurrLoc.size() << "\n";
			NewArray.insert(NewArray.end(), CurrLoc.begin(), CurrLoc.end()); //concatenate the new alleles to the growing list for this locus
		}
		TempList.push_back(NewArray);
	}
	MaxAlleles = TempList;  //the maximum number of alleles possible at each locus, for standardization
	
	//3. calculate diversity (M) for each locus separately, meanwhile adding it up across all loci
	
	//cout << "MaxAlleles.size()=" << MaxAlleles.size() << "\n";
	
	for (i=0;i<MaxAlleles.size();i++)
	{
		ListToFilter = MaxAlleles[i];
		TempList2 = MyFilterDuplicatesII(ListToFilter);
		Mtemp =  TempList2.size();
		ActiveMaxAllelesList.push_back(Mtemp);
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


	if (NumLoci == 0)
	{
		RandomActiveDiversity = -1; 
		AltRandomActiveDiversity = -1;
	}
	else
	{
		for (i=0;i<NumLoci;i++)
		{
			//3. fuse alleles from the same locus into a single array, for all populations in core
			//vector<std::string>().swap(NewArray); //clear NewArray
			NewArray.clear();
			for (j=0;j<CoreSize;j++)
			{
				CurrLoc = AlleleList[j][i];
				NewArray.insert(NewArray.end(), CurrLoc.begin(), CurrLoc.end());
			}
			
			/*for (k=0;k<NewArray.size();k++)
			{
				cout << "NewArray["<<k<<"]="<<NewArray[k]<<"\n";
			}*/
			
			//4. assemble a list of diversity (M) for each locus separately
			//ListToFilter starts with a sorted list of all unique alleles found at a locus.  If there are no alleles (i.e. all missing data), M is set to 0.  Otherwise, M is set to the number of unique alleles at the locus.
			//vector<std::string>().swap(ListToFilter); //clear ListToFilter
			ListToFilter.clear();
			//remove duplicates
			ListToFilter = NewArray;
			sort( ListToFilter.begin(), ListToFilter.end() );
			ListToFilter.erase( std::unique( ListToFilter.begin(), ListToFilter.end() ), ListToFilter.end() );

			
			/*cout << "Locus "<<i<<"\n";
			for (k=0;k<ListToFilter.size();k++)
			{
				cout << "ListToFilter["<<k<<"]="<<ListToFilter[k]<<"\n";
			}*/

			if (ListToFilter.size() == 0)
			{
				M=0;
			}
			else
			{
				M=ListToFilter.size(); //M=the number of alleles present
			}
			
			//Mlist.push_back(M); //Mlist contains number of alleles present at each locus
			Mlist[i] = M;
		}
		
		
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
	
	//declare variables
	int i, j, k, l;
	string DoSummary = "no"; //switch to turn on summary function
	string Kernel = "no"; //switch to include a mandatory set in the core
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
			if (fileExists(SumFilePath) == 0) 
			{
				bf = "SumFilePath = ";
				bf += SumFilePath;
				BadFiles.push_back(bf);
			}
		} 
		
		if ( string(argv[i]) == "-k" ) 
    	{
        	Kernel = "yes";
        	KerFilePath = argv[i+1];
        	KernelAccessionList = MySetKernel(KerFilePath);
			if (fileExists(KerFilePath) == 0) 
			{
				bf = "KerFilePath = ";
				bf += KerFilePath;
				BadFiles.push_back(bf);
			}
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
	
	
	
	//start the clock
	time_t start,end;
	time (&start);

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
	
	//PROCESS INPUT DATA
	//.var file
	vector<int> AllColumnIDList;
	vector<std::string> AllLociNameList;
	vector<int> ActiveColumnIDList;
	vector<std::string> ActiveLociNameList;
    vector<int> TargetColumnIDList;
	vector<std::string> TargetLociNameList;
	vector<vector<int> > ColKeyToAllAllelesByPopList;
	vector<int> ReferenceOrTargetKey;

	
	MyProcessVarFile(VarFilePath, AllColumnIDList, AllLociNameList, ActiveColumnIDList, ActiveLociNameList, TargetColumnIDList, TargetLociNameList, ColKeyToAllAllelesByPopList, ReferenceOrTargetKey ); 
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
	vector<vector<vector<std::string> > > AllAlleleByPopList; //structure of this 3D vector is:
	// { { {pop1,loc1 alleles},{pop1,loc2 alleles},...}, { {pop2,loc1 alleles},{pop2,loc2 alleles},...} } }
	vector<std::string> FullAccessionNameList;
	
	MyProcessDatFile(DatFilePath, AllColumnIDList, AllLociNameList, AllAlleleByPopList, FullAccessionNameList);
	//latter 3 variables updated as reference
	
		/*
		//Print out FullAccessionNameList
		cout << "\n\nPopulation names\n";
		for (i=0; i<FullAccessionNameList.size();i++) cout << FullAccessionNameList[i] << "\n";
		
		//Print out lists of unique All alleles from AllAlleleByPopList
		for (i=0;i<AllAlleleByPopList.size() ;i++)
		{
			cout << "Population " << i << "\n";
			for (j=0;j<AllAlleleByPopList[i].size();j++)
			{
				cout << "Locus " << j << "\n";
				for (k=0;k<AllAlleleByPopList[i][j].size();k++)
				{
					cout << AllAlleleByPopList[i][j][k] << ",";
				}
				cout << "\n";
			}
		
		}
		*/
	
	
	
	
	
	
	//sort AllAlleleByPopList into reference and target loci lists
	vector<vector<vector<std::string> > > ActiveAlleleByPopList; 
	vector<vector<vector<std::string> > > TargetAlleleByPopList; 
	//reference
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

	
	
/*	THIS IS THE OLD ROUTINE, WHICH DIDN'T ALLOW POLYPLOIDY AND RAN THROUGH MYPROCESSDATFILE TWICE (ONCE FOR REFERENCE, ONCE FOR TARGET)

	//.dat file, active loci
	vector<vector<vector<std::string> > > ActiveAlleleByPopList; //structure of this 3D vector is:
	// { { {pop1,loc1 alleles},{pop1,loc2 alleles},...}, { {pop2,loc1 alleles},{pop2,loc2 alleles},...} } }
	vector<std::string> FullAccessionNameList;
	
	MyProcessDatFile(DatFilePath, ActiveColumnIDList, ActiveLociNameList, ActiveAlleleByPopList, FullAccessionNameList);
	//MyProcessDatFileII(DatFilePath, ActiveColumnIDList, ActiveLociNameList, ActiveAlleleByPopList, FullAccessionNameList);
	//latter two variables above are updated as references
		
		//Print out FullAccessionNameList
		cout << "\n\nPopulation names\n";
		for (i=0; i<FullAccessionNameList.size();i++) cout << FullAccessionNameList[i] << "\n";
		
		//Print out lists of unique active alleles from ActiveAlleleByPopList
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
	
	//.dat file, target loci
	vector<vector<vector<std::string> > > TargetAlleleByPopList; //structure of this 3D vector is:
	// { { {pop1,loc1 alleles},{pop1,loc2 alleles},...}, { {pop2,loc1 alleles},{pop2,loc2 alleles},...} } }
	vector<std::string> TargetAccessionNameList;
	
	MyProcessDatFile(DatFilePath, TargetColumnIDList, TargetLociNameList, TargetAlleleByPopList, TargetAccessionNameList);
	//latter two variables above are updated as references. TargetAccessionNameList is never used,
	//just created here and sent as a dummy so the function doesn't update FullAccessionNameList
		
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
	
	//get number of loci
	int NumLoci = ActiveAlleleByPopList[1].size();
	
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
		
	
	cout << "\n\n" << "Number of accessions = " << NumberOfAccessions << ", Number of loci = " << NumLoci;
	if (Kernel == "yes") cout << ", Number of kernel accessions = " << KernelAccessionList.size() << "\n\n";
	else cout << "\n\n";
	

	//catch some errors
	if (MaxCoreSize > NumberOfAccessions)
	{
		cout << "ERROR:  Maximum core size of "<<MaxCoreSize<< " is greater than the number of accessions.  Please set to a lower value. Quitting...\n\n";
		exit (EXIT_FAILURE);
	}
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

    /*cout << "ActiveMaxAllelesList:\n";
	for (i=0; i<ActiveMaxAllelesList.size(); i++) cout << ActiveMaxAllelesList[i] << ",\n";
    cout << "\nTargetMaxAllelesList:\n";
	for (i=0; i<TargetMaxAllelesList.size(); i++) cout << TargetMaxAllelesList[i] << ",\n";
	*/
		
	

	//COMMENCE M+
	
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
		
	
	
	cout << "Optimizing...\n";


	
	//compile as parallel or not?
	int parallelism_enabled = 1; //0=no, not 0 = yes
	
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
	time (&end);
	double dif = difftime (end,start);
	cout << "\nElapsed time = "<< dif << " seconds.\n\n";





}
