#include "m+.hpp"

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

int MyProcessVarFile(char* VarFileBuffer, vector<int>& AllColumnIDList, vector<std::string>& AllLociNameList, vector<int>& ActiveColumnIDList, vector<std::string>& ActiveLociNameList, vector<int>& TargetColumnIDList, vector<std::string>& TargetLociNameList, vector<vector<int> >& ColKeyToAllAlleleByPopList, vector<int>& ReferenceOrTargetKey, vector<int>& PloidyList, vector<std::string>& UniqLociNameList)
{
    //declare variables
    std::string foo;
    vector<std::string> foovector;
    unsigned int k;

    int i=0; // i is the row number

	std::string vfb(VarFileBuffer); //convert char* to string
    std::istringstream iss(vfb); //open string as a stream
    while (std::getline(iss, foo)) //cycle line by line, placing line in foo
    {
	    //split foo on whitespace
		foovector = split(foo);
		
		if (foovector.size() != 0) //omit the hanging last line, which has no information
		{
		
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
    }
	
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
int MyReduceToRef(vector<vector<vector<int> > > AllAlleleByPopList, vector<int> ReferenceOrTargetKey, vector<vector<vector<int> > >& ActiveAlleleByPopList, vector<vector<vector<int> > >& TargetAlleleByPopList)
{
	unsigned int r, t, i, j;
	vector<int> b;
	
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
			vector<int>().swap(b);
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

int MyProcessDatFileIII(char* DatFileBuffer, vector<int> AllColumnIDList, vector<vector<int> > ColKeyToAllAlleleByPopList, vector<vector<vector<int> > >& AllAlleleByPopList, vector<std::string>& FullAccessionNameList, vector<std::string>& IndivPerPop, vector<int>& AllAlleles, std::string Rarify)
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
	
	//initialize important data structures, unsized
	vector<vector<set<int> > > AllAlleleByPopListSet; //structure of this 3D vector is:
	// { { {pop1,loc1 alleles},{pop1,loc2 alleles},...}, { {pop2,loc1 alleles},{pop2,loc2 alleles},...} } }
	// this is the same structure as the vector<vector<vector<int> > > AllAlleleByPopList
	// if Rarify == "yes", AllAlleleByPopList just becomes ByPop3d.  if Rarify == "no", ByPop3d is reduced
	// to unique alleles at each locus using AllAlleleByPopListSet, then converted back into the vector
	// AllAlleleByPopList, which is updated as a reference


	unsigned int i,j,k,l;
	vector<std::string> bufvec;
	std::string NewPopID;
	std::string OldPopID = "init*@#rt4"; //use an unlikely population name for the initialization value

	cout << "  Reading data matrix...\n";
	//start the clock
	time_t startm,endm;
	time (&startm);

	//put giant char array buffer into a stream, then delete the buffer to save memory 
	stringstream s(DatFileBuffer);
	strcpy(DatFileBuffer,""); //clear char* 
	
	
	//read buffer into a vector, one line per item
	while (getline(s, foo))//get line from s, put in foo, consecutively
	{
		bufvec.push_back(foo);  
	}
	
	s.str(""); //clear stringstream
	int row = bufvec.size();
		
	//sort vector so that individuals from the same population form consecutive elements
	std::sort(bufvec.begin(), bufvec.end()); //no need to use fancy sort, lexicographic should be fine
		
	//split lines of bufvec into 2d vector
	vector<vector<std::string> > bufvec2d(bufvec.size()); 
	for (i=0;i<bufvec.size();++i) bufvec2d[i] = split(bufvec[i]);
	vector<std::string>().swap(bufvec); //clear bufvec
	
	//stop the clock
	time (&endm);
	double dif = difftime (endm,startm);
	if (dif==1) cout << "    " << dif << " second.\n";	
	else cout << "    " << dif << " seconds.\n";
	
	//convert alleles to integer coding to save memory, vector access order 2
	cout << "  Recoding data...\n";
	time (&startm);

	vector<vector<int> > bufvec2dint(bufvec2d.size(), vector<int>(bufvec2d[0].size())); //declare and size vector to hold new integer coded alleles
	unsigned int iz = ColKeyToAllAlleleByPopList.size();
	for (unsigned int i=0;i<iz;++i) //go thru each locus
	{
		//get all alleles at the locus
		vector<std::string> AllelesEncountered; //will contain the unique set of alleles at the locus
		unsigned int kz = bufvec2d.size();
		for (unsigned int k=0;k<kz;++k) //go thru all individuals
		{
			unsigned int jz = ColKeyToAllAlleleByPopList[i].size();
			for (unsigned int j=0;j<jz;++j) //get all alleles for an individual at this locus, then move to next indiv
			{
				int ColIndex = ColKeyToAllAlleleByPopList[i][j];
				std::string a = bufvec2d[k][ColIndex];  
				if (a == "9999") bufvec2dint[k][ColIndex] = -9999; //add the missing data value
				else
				{
					int AlleleInt; //the new, integerized, name of the allele
					std::vector<std::string>::iterator itr = std::find(AllelesEncountered.begin(), AllelesEncountered.end(), a);
					if (itr != AllelesEncountered.end()) //the allele has been found before
					{
						AlleleInt = itr - AllelesEncountered.begin(); //convert itr to index, the index is the integerized allele name
						bufvec2dint[k][ColIndex] = AlleleInt; //add the new name
					}
					else // you have a new allele
					{
						AllelesEncountered.push_back(a); //add new allele to list of those encountered
						AlleleInt = AllelesEncountered.size() - 1;  //calculate integerized allele name, starts at 0
						bufvec2dint[k][ColIndex] = AlleleInt;
					}
				}
			}
		}
	}
	
	//stop the clock
	time (&endm);
	dif = difftime (endm,startm);
	if (dif==1) cout << "    " << dif << " second.\n";	
	else cout << "    " << dif << " seconds.\n";	

	cout << "  Building data structures...\n";
	time (&startm);
	
	//simplify bufvec2d by retaining only the first column, with POPID, then swap clear to save memory
	vector<string> PopIDvec;
	for (i=0;i<bufvec2d.size();++i)
		PopIDvec.push_back(bufvec2d[i][0]);
	vector<vector<std::string> >().swap(bufvec2d); //clear bufvec2d

	//break up vector into a 3d vector by population:  { { {pop1ind1elems},{pop1ind2elems},...}, { {pop2ind1elems},{pop2ind2elems},...} }
	vector<vector<vector<int> > > ByPop3d;
	for (i=0;i<bufvec2dint.size();++i)
	{
		//NewPopID = bufvec2d[i][0];  //access string in bufvec2d
		NewPopID = PopIDvec[i];
		IndivPerPop.push_back(NewPopID); //add the pop ID to a list to calc pop sizes later

		if (NewPopID != OldPopID) //then create a new population in ByPop3D
		{
			ByPop3d.resize(ByPop3d.size() + 1);
			
			//add the new population name to the AccessionNameList
			FullAccessionNameList.push_back(NewPopID);
		}
		
		//push row of integer elements on current line, onto last item of ByPop3d, which might be a new population as added just above
		//remember that the first three columns are 0 in bufvec2dint, so will also be 0 in ByPop3d
		ByPop3d[ByPop3d.size()-1].push_back(bufvec2dint[i]);
		
		OldPopID = NewPopID;
	}
	vector<vector<int> >().swap(bufvec2dint); //clear bufvec2dint, the integerized data is in ByPop3d

	//stop the clock
	time (&endm);
	dif = difftime (endm,startm);
	if (dif==1) cout << "    " << dif << " second.\n";	
	else cout << "    " << dif << " seconds.\n";	

/*	//print out ByPop3d
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
	
	}
*/

	cout << "  Condensing data...\n";
	time (&startm);

	//When Rarify=yes, information on allele frequencies must be retained so all alleles are
	//placed into the AllAlleleByPopList (excluding -9999 missing data).  When Rarify=no, 
	//all alleles are passaged through a set (AllAlleleByPopListSet) so that only unique 
	//alleles end up in AllAlleleByPopList.
	if (Rarify == "yes")
	{
		//resize AllAlleleByPopList
		AllAlleleByPopList.resize(ByPop3d.size());//resize number of populations
		for (i=0;i<AllAlleleByPopList.size();++i)
		{
			AllAlleleByPopList[i].resize(ColKeyToAllAlleleByPopList.size()); //resize number of loci
		}
		
		//calculate size of AllAlleles
		AllAlleles.reserve(row*AllColumnIDList.size());

		//condense alleles by locus in AllAllelesByPopList, within each population
		int AlleleIndex;
		for (i=0;i<ByPop3d.size();++i) //go thru pops
		{
			for (j=0;j<ByPop3d[i].size();++j) //go thru indivs
			{
				for (k=0;k<ColKeyToAllAlleleByPopList.size();++k) //go through each locus
				{
					for (l=0;l<ColKeyToAllAlleleByPopList[k].size();++l) //assign columns to loci
					{
						AlleleIndex = ColKeyToAllAlleleByPopList[k][l];
						int NewAllele = ByPop3d[i][j][AlleleIndex]; //get the allele in the specified column
						AllAlleles.push_back(NewAllele); //add the allele to the list of all alleles, missing data included
						if (NewAllele != -9999) //exclude missing data
							AllAlleleByPopList[i][k].push_back(NewAllele); //add the allele to the set of unique alleles at locus k, pop i	
					}
				}
			}
		}
		
 		//stop the clock on Condensing data...
		time (&endm);
		dif = difftime (endm,startm);
		if (dif==1) cout << "    " << dif << " second.\n";	
		else cout << "    " << dif << " seconds.\n";	
	}
	else if (Rarify == "no")
	{
		//resize AllAlleleByPopListSet
		AllAlleleByPopListSet.resize(ByPop3d.size());//resize number of populations
		for (i=0;i<AllAlleleByPopListSet.size();++i)
		{
			AllAlleleByPopListSet[i].resize(ColKeyToAllAlleleByPopList.size()); //resize number of loci
																		   //the index in ColKey is the locus index in AllAllelesByPopList level 2
																		   //the value of ColKey is the index of the allele in ByPop3d level 3
		}
		
		//calculate size of AllAlleles
		AllAlleles.reserve(row*AllColumnIDList.size());

		//condense alleles by locus in AllAllelesByPopListSet, within each population
		int AlleleIndex;
		for (i=0;i<ByPop3d.size();++i) //go thru pops
		{
			for (j=0;j<ByPop3d[i].size();++j) //go thru indivs
			{
				for (k=0;k<ColKeyToAllAlleleByPopList.size();++k) //go through each locus
				{
					for (l=0;l<ColKeyToAllAlleleByPopList[k].size();++l) //assign columns to loci
					{
						AlleleIndex = ColKeyToAllAlleleByPopList[k][l];
						int NewAllele = ByPop3d[i][j][AlleleIndex]; //get the allele in the specified column
						AllAlleles.push_back(NewAllele); //add the allele to the list of all alleles, missing data included
						if (NewAllele != -9999) //exclude missing data
							AllAlleleByPopListSet[i][k].insert(NewAllele); //add the allele to the set of unique alleles at locus k, pop i	
					}
				}
			}
		}
	
		//stop the clock on Condensing data...
		time (&endm);
		dif = difftime (endm,startm);
		if (dif==1) cout << "    " << dif << " second.\n";	
		else cout << "    " << dif << " seconds.\n";	
		
		//start the clock on Converting...
		time_t startd,endd;
		cout << "  Converting data structures...\n";
		time (&startd);

		//resize AllAlleleByPopList
		AllAlleleByPopList.resize(ByPop3d.size());//resize number of populations
		for (i=0;i<AllAlleleByPopList.size();++i)
		{
			AllAlleleByPopList[i].resize(ColKeyToAllAlleleByPopList.size()); //resize number of loci
		}
		//convert set to vector for further processing
		for (i=0;i<AllAlleleByPopList.size();++i)
		{
			for (j=0;j<AllAlleleByPopList[i].size();++j)
			{
				vector<int> ttvec(AllAlleleByPopListSet[i][j].begin(), AllAlleleByPopListSet[i][j].end()); //use constructor to convert set to vector	
				AllAlleleByPopList[i][j] = ttvec;
			}
		}
		vector<vector<set<int> > >().swap(AllAlleleByPopListSet); //clear variable, no longer needed
	
		//stop the clock on Converting...
		time (&endd);
		dif = difftime (endd,startd);
		if (dif==1) cout << "    " << dif << " second.\n";	
		else cout << "    " << dif << " seconds.\n";	
	}
	
	vector<vector<vector<int> > >().swap(ByPop3d); //clear ByPop3d
			
	return 0;
}

vector<unsigned int> Mysss(vector<vector<vector<int> > > AlleleByPopList)
{
	vector<unsigned int> sss;
	unsigned int s;
	for (unsigned int i=0;i<AlleleByPopList[0].size();++i)
	{
		//find the smallest non-zero sample (~= smallest pop size) for current locus, must look at all loci due to missing data
		s = std::numeric_limits<int>::max(); //set initial s to largest integer possible
		for (unsigned int j=0;j<AlleleByPopList.size();j++)
		{
			//test whether current sample size is lowest so far
			if (AlleleByPopList[j][i].size() <= s) s = AlleleByPopList[j][i].size();
			if (s == 0) break;
			//below does not allow s=0
			//if (ActiveAlleleByPopList[j][i].size() == 0) continue;
			//else if (ActiveAlleleByPopList[j][i].size() <= s) s = ActiveAlleleByPopList[j][i].size();
		}
		sss.push_back(s);  //add the smallest non-zero sample size for this locus to the list of such things.
	}
	return sss;
}

//returns maximum number of alleles possible at each locus for active and target
vector<int> MyGetMaxs(vector<vector<vector<int> > > ActiveAlleleByPopList, std::mt19937_64& rng)
{
	unsigned int i, j, k;
	vector<int> ActiveMaxAllelesList;
	set<int> NewSet;

	//determine unique alleles at each locus using a set
	for (i=0;i<ActiveAlleleByPopList[0].size();++i)
	{
		NewSet.clear();
		//traverse the locus 'column' of the 3d grid, get locus i for population j
		for (j=0;j<ActiveAlleleByPopList.size();j++)
		{
			for (k=0;k<ActiveAlleleByPopList[j][i].size();++k)
			{
				NewSet.insert(ActiveAlleleByPopList[j][i][k]);	//place alleles into set to eliminate redundancies
			}
		}
		ActiveMaxAllelesList.push_back(NewSet.size()); //the set size after adding all populations for locus i is the maximum number of alleles
	}
	return ActiveMaxAllelesList;
}

vector<int> MySubSampleData(vector<int> bb, unsigned int s, std::mt19937_64& rng)
{
	//bb is the incoming vector<int> of alleles
	//s is the non-zero minimum sample size at the locus, for rarefaction
	vector<int> aa; //rarefied vector<int> of alleles
	
	std::shuffle(bb.begin(), bb.end(), rng); //shuffle the order of the alleles so that a random sample of sss alleles can be taken for rarefaction

	for (unsigned int k=0;k<s;++k)
	{
		aa.push_back(bb[k]); //place alleles into set to eliminate redundancies, sample no more than sss[i] alleles per population
	}
	return aa;
}

//subsamples Active(Target)AlleleByPopList, reducing to smallest non-zero sample size
vector<vector<vector<int> > > MyRarifyData(vector<vector<vector<int> > > AlleleByPopList, vector<unsigned int> sv, std::mt19937_64& rng)
{
	//create receive vector, of same size as Active(Target)AlleleByPopList
	vector<vector<vector<int> > > recvec(AlleleByPopList.size()); //size only the first level
	vector<int> aa; //will hold rarefied allele vector<int>
	//unsigned long long f;
	
	//rarify each Pop x locus x allele vector independently, to make MPI_Bcast easier
	for (unsigned int i=0;i<AlleleByPopList.size();++i)
	{
		for (unsigned int j=0;j<AlleleByPopList[0].size();++j)
		{
			//f=0; //vector<int> size
			aa = MySubSampleData(AlleleByPopList[i][j], sv[j], rng);
			//f = (unsigned long long)aa.size();
			
			//printf("i=%d j=%d rank=%02d f=%d aasize=%d\n", i, j, f, aa.size());
			
			//add vector<int> onto recvec
			recvec[i].push_back(aa);				
			vector<int>().swap(aa); //clear aa
		}
	}
	//copy recvec to ActiveAlleleByPopList
	AlleleByPopList = recvec;
	vector<vector<vector<int> > >().swap(recvec); //clear recvec
	
	return AlleleByPopList;
}

bool fileExists(const char *fileName)
{
    ifstream infile(fileName);
    return infile.good();
}

vector<int> MyRemoveTargetAlleles(vector<int> AllAlleles, vector<int> AllColumnIDList, vector<int> TargetColumnIDList)
{
	unsigned int checker, j;
	unsigned int i=0;
	vector<int> AllRefAlleles;
	
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
void MyMakeRefAllelesIntoRows(vector<int> AllRefAlleles, vector<std::string> ActiveLociNameList, vector<vector<int> >& RefAllelesIntoRows)
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
void MyMakeRefAllelesByLocus(vector<vector<int> > RefAllelesIntoRows, vector<string> ActiveLociNameList, vector<std::pair<std::string, vector<int> > >& RefAllelesByLocus)
{
	vector<std::pair<std::string, vector<int> > > lola;
	std::pair<std::string, vector<int> > la; //locus name, allele list pair
	
	unsigned int i, j, k;
	int locindex;
	std::string locname, b;
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
			int al = RefAllelesIntoRows[k][i];
			if (al != -9999) lola[locindex].second.push_back(al); //ignore missing data
		}
	}
	
	//update RefAllelesByLocus
	RefAllelesByLocus = lola;

/*	//print out each locus name followed by all the alleles found within it
		vector<int> foo;
		for (i=0;i<lola.size();++i)
		{
			cout << lola[i].first << " ";
			foo = lola[i].second;
			for (j=0;j<foo.size();++j)
			{
				cout << " " << foo[j];
			}
			cout << "\n";
		}
*/
	
	vector<std::pair<std::string, vector<int> > >().swap(lola); //clear lola
}

//calculates allele frequencies for all alleles at all loci, updates vector of struct Alfreq, which contains the relational data
void MyCalculateAlleleFrequencies(vector<std::pair<std::string, vector<int> > > RefAllelesByLocus, vector<Alfreq>& AlleleFrequencies)
{
	unsigned int i, j, z;
	double freq;
	std::string b;
	vector<int> AllAlleles;
	vector<int> UniqAlleles;
	set<int> AlleleSet;

	Alfreq laf; //locusname, allelenames, frequencies
	static const struct Alfreq emptylaf; //this will be used to zero struct between loops
	
	for (i=0;i<RefAllelesByLocus.size();++i)
	{
		//empty containers
		laf = emptylaf;
		vector<int>().swap(UniqAlleles); //clear UniqAlleles
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
			int al = UniqAlleles[j];
			laf.allelenames.push_back(al); //add allele name to struct
			z = count(AllAlleles.begin(), AllAlleles.end(), al);
			freq = double(z)/double(AllAlleles.size());
			laf.frequencies.push_back(freq); //add allele frequency to struct
		}
		
		AlleleFrequencies.push_back(laf);
	}
}

//determine the size of a file on disk
std::ifstream::pos_type filesize(const char* filename)
{
	std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
	return in.tellg(); 
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
	string Rarify = "no"; //switch to use rarefaction to correct allele count for sample size differences
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

		if ( string(argv[i]) == "-r" ) 
    	{
        	Rarify = "yes";
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
        if(VarFilePath != NULL){
		    bf += VarFilePath;
		    BadFiles.push_back(bf);
        }
        else{
            BadFiles.push_back("NO VarFile");
        }
	}
 
	if (fileExists(DatFilePath) == 0) 
	{
		bf = "DatFilePath = ";
        if( DatFilePath != NULL){
		    bf += DatFilePath;
		    BadFiles.push_back(bf);
        }
        else{
            BadFiles.push_back("NO DatFile");
        }
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
	if (Rarify == "yes")
	{
		cout << "  -r invoked:\n";
		cout << "    Using rarefaction for M+ search\n";
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
	if (DoM == "yes")
	{
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
	}
	if (Kernel == "yes" && Ideal == "yes" && DoM == "no")
	{
		cout << "ERROR:  A* search cannot be performed with a kernel file.  Please correct the command line.  Quitting...\n\n";
		exit (EXIT_FAILURE);
	}

	//DETERMINE MACHINE CONFIGURATIION
	int ncpu = sysconf( _SC_NPROCESSORS_ONLN );

	//PROCESS INPUT DATA
	
	//start the clock
	time_t starti,endi;
	time (&starti);

	//read the file into a buffer using fread
	unsigned long long f = 0; //var (and later, dat) file size
	f = (unsigned long long)filesize(VarFilePath); 
	f = f + 1; //increase by 1 to accomodate \0 string terminator
	char * VarFileBuffer = (char*)malloc(f); //initialize and size the data structure to hold contents of var file

	char * d = MyBigRead(VarFilePath);
	strcpy(VarFileBuffer, d); //convert the char* to a char array, presized from prior MPI_Bcast
	strcpy(d,""); //clear char*

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

	MyProcessVarFile(VarFileBuffer, AllColumnIDList, AllLociNameList, ActiveColumnIDList, ActiveLociNameList, TargetColumnIDList, TargetLociNameList, ColKeyToAllAlleleByPopList, ReferenceOrTargetKey, PloidyList, UniqLociNamesList ); 
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

	//read the whole file into a buffer using fread
	f=0; //dat file size
	f = (unsigned long long)filesize(DatFilePath); 
	f = f + 1; //increase by 1 to accomodate \0 string terminator
	char * DatFileBuffer = (char*)malloc(f); //initialize and size the data structure to hold contents of dat file

	d = MyBigRead(DatFilePath);
	strcpy(DatFileBuffer, d); //convert the char* to a char array
	strcpy(d,""); //clear char*
	
	vector<std::string> FullAccessionNameList;
	vector<std::string> IndivPerPop;
	vector<int> AllAlleles;
	
	vector<vector<vector<int> > > AllAlleleByPopList; //structure of this 3D vector is:
	//vector<vector<vector<int> > > AllAlleleByPopList( AllAlleleByPopListSet.size(), vector<vector<int> >(UniqLociNamesList.size()) ); //structure of this 3D vector is:
	// { { {pop1,loc1 alleles},{pop1,loc2 alleles},...}, { {pop2,loc1 alleles},{pop2,loc2 alleles},...} } }
	
	//Process the dat file
	MyProcessDatFileIII(DatFileBuffer, AllColumnIDList, ColKeyToAllAlleleByPopList, AllAlleleByPopList, FullAccessionNameList, IndivPerPop, AllAlleles, Rarify);
	//AllAlleleByPopList, FullAccessionNameList, IndivPerPop, AllAlleles updated by reference

	time_t startd,endd;
	double dif = difftime (endd,startd);
	cout << "  Separating reference and target loci...\n";
	time (&startd);

	//sort AllAlleleByPopList into reference and target loci lists
	vector<vector<vector<int> > > ActiveAlleleByPopList; 
	vector<vector<vector<int> > > TargetAlleleByPopList; 
	MyReduceToRef(AllAlleleByPopList, ReferenceOrTargetKey, ActiveAlleleByPopList, TargetAlleleByPopList); //latter 2 variables updated as reference

	/*
		//Print out alleles from AllAlleleByPopList
			vector<int> si;
			for (i=0;i<AllAlleleByPopList.size() ;i++)
			{
				cout << "Population " << FullAccessionNameList[i] << "\n";
				for (j=0;j<AllAlleleByPopList[i].size();j++)
				{
					cout << "  Locus " << j << "\n    ";
					cout << "  AllAlleleByPopList[" << i << "][" << j << "].size()=" << AllAlleleByPopList[i][j].size() << "\n";
					si = AllAlleleByPopList[i][j];
					for (std::vector<int>::iterator it=si.begin(); it!=si.end(); ++it)
						cout << *it << ",";
					cout << "\n";
				}
			}
			
		//Print out FullAccessionNameList
		cout << "\n\nPopulation names\n";
		for (unsigned int i=0; i<FullAccessionNameList.size();i++) cout << FullAccessionNameList[i] << "\n";
		
		//Print out lists of unique reference alleles from ActiveAlleleByPopList
		cout << "ActiveAlleleByPopList:\n";
		for (unsigned int i=0; i<ActiveAlleleByPopList.size() ;i++)
		{
			cout << "Population " << i << "\n";
			for (unsigned int j=0;j<ActiveAlleleByPopList[i].size();j++)
			{
				cout << "Locus " << j << "\n";
				for (unsigned int k=0;k<ActiveAlleleByPopList[i][j].size();k++)
				{
					cout << ActiveAlleleByPopList[i][j][k] << ",";
				}
				cout << "\n";
			}
		
		}
		
		//Print out lists of unique target alleles from TargetAlleleByPopList
		cout << "TargetAlleleByPopList:\n";
		for (unsigned int i=0; i<TargetAlleleByPopList.size() ;i++)
		{
			cout << "Population " << i << "\n";
			for (unsigned int j=0;j<TargetAlleleByPopList[i].size();j++)
			{
				cout << "Locus " << j << "\n";
				for (unsigned int k=0;k<TargetAlleleByPopList[i][j].size();k++)
				{
					cout << TargetAlleleByPopList[i][j][k] << ",";
				}
				cout << "\n";
			}
		}
	*/	

	vector<vector<vector<int> > >().swap(AllAlleleByPopList); //clear variable, no longer needed

	time (&endd);
	dif = difftime (endd,startd);
	if (dif==1) cout << "    " << dif << " second.\n";	
	else cout << "    " << dif << " seconds.\n";	

	//CALCULATE SOME USEFUL VARIABLES
	//get total number of accessions
	unsigned int NumberOfAccessions = ActiveAlleleByPopList.size();
	
	//catch a possible error in the command line
	if (DoM == "yes")
	{
		if (MaxCoreSize > NumberOfAccessions)
		{
			cout << "ERROR:  Maximum core size of "<<MaxCoreSize<< " is greater than the number of accessions.  Please set to a lower value. Quitting...\n\n";
			exit (EXIT_FAILURE);
		}
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
	}
	
	cout << "  Calculating run specific parameters...\n";
	time (&startd);

	//1. remove target alleles from AllAlleles
	vector<int> AllRefAlleles = MyRemoveTargetAlleles(AllAlleles, AllColumnIDList, TargetColumnIDList);
	vector<int>().swap(AllAlleles); //clear variable, no longer needed
	
	//2. put AllRefAlleles into a 2d vector, samples are rows
	vector<vector<int> > RefAllelesIntoRows(IndivPerPop.size(), vector<int> ( ActiveLociNameList.size() ));
	MyMakeRefAllelesIntoRows(AllRefAlleles, ActiveLociNameList, RefAllelesIntoRows);
	vector<int>().swap(AllRefAlleles); //clear variable, no longer needed

	//3. extract alleles into vector of pairs (ignores missing data -9999)
	vector<std::pair<std::string, vector<int> > > RefAllelesByLocus; // first = locus name, second = vector of all alleles present, updated as reference below
	MyMakeRefAllelesByLocus(RefAllelesIntoRows, ActiveLociNameList, RefAllelesByLocus);
	vector<vector<int> >().swap(RefAllelesIntoRows); //clear variable, no longer needed
	
	//4. calculate allele frequencies, finally (ignores missing data -9999)
	vector<Alfreq> AlleleFrequencies;//(RefAllelesByLocus.size());  //declare vector of struct Alfreq
	MyCalculateAlleleFrequencies(RefAllelesByLocus, AlleleFrequencies);
	vector<std::pair<std::string, vector<int> > >().swap(RefAllelesByLocus); //clear variable, no longer needed
	
	time (&endd);
	dif = difftime (endd,startd);
	if (dif==1) cout << "    " << dif << " second.\n";	
	else cout << "    " << dif << " seconds.\n";	

	cout << "  Finalizing data structures...\n";
	time (&startd);

	/*
		//print out structs containing allele frequencies
		Alfreq laf;
		vector<int> anames;
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

	//seed the random number generator
	unsigned long long seed = (unsigned long long)chrono::high_resolution_clock::now().time_since_epoch().count();
	std::mt19937_64 rng(seed);    // random-number engine used (Mersenne-Twister in this case), seeded with time

	//get maximum number of alleles possible at each locus for active and target
	vector<unsigned int> sss;  //sss is a vector holding the smallest non-zero sample size for each Active locus, only relevant to Rarify=yes
	vector<unsigned int> sst;  //sst is a vector holding the smallest non-zero sample size for each Target locus, only relevant to Rarify=yes
	vector<int> ActiveMaxAllelesList, TargetMaxAllelesList;
	sss = Mysss(ActiveAlleleByPopList);
	sst = Mysss(TargetAlleleByPopList);
	ActiveMaxAllelesList = MyGetMaxs(ActiveAlleleByPopList, rng);
	TargetMaxAllelesList = MyGetMaxs(TargetAlleleByPopList, rng);

	/*
		//Print out alleles from TargetAlleleByPopList
			vector<int> si;
			for (i=0;i<TargetAlleleByPopList.size() ;i++)
			{
				cout << "PrePopulation " << FullAccessionNameList[i] << "\n";
				for (j=0;j<TargetAlleleByPopList[i].size();j++)
				{
					cout << "PreLocus " << j << ", ";
					cout << "TargetAlleleByPopList[" << i << "][" << j << "].size()=" << TargetAlleleByPopList[i][j].size() << ", ";
					si = TargetAlleleByPopList[i][j];
					for (std::vector<int>::iterator it=si.begin(); it!=si.end(); ++it)
						cout << *it << ",";
					cout << "\n";
				}
			}
	*/

	//If Rarify="yes", reduce the Active(Target)AlleleByPopList to the smallest non-zero sample size for each locus
	//This is the rarification step. It is performed by master node then passed to other procids 
	//using MPI_Bcast so that the same rarefied data set is used by all procs.
	if (Rarify == "yes")
	{
		cout << "    Rarefaction...\n";
		
		//subsample the Active(Target)MaxAllelesList
		ActiveAlleleByPopList = MyRarifyData(ActiveAlleleByPopList, sss, rng);
		TargetAlleleByPopList = MyRarifyData(TargetAlleleByPopList, sst, rng);

		//test whether all smallest sample sizes are zero for ActiveAlleleByPopList after rarefaction
		bool zerovec = std::any_of(sss.begin(), sss.end(), [](unsigned int i) { return i==0; });
		if (zerovec)
		{
			cout << "\nERROR:  Some sample sizes = 0 after rarefaction.  Remove loci with missing data for all members of an accession.  Quitting...\n\n"; 
			exit (EXIT_FAILURE); //master0 reports above, everybody quits here
		}
		
		//recalculate Active(Target)MaxAllelesList
		ActiveMaxAllelesList = MyGetMaxs(ActiveAlleleByPopList, rng);
		TargetMaxAllelesList = MyGetMaxs(TargetAlleleByPopList, rng);
	}
	
	/*
		//Print out alleles from ActiveAlleleByPopList
			vector<int> si;
			for (i=0;i<ActiveAlleleByPopList.size() ;i++)
			{
				cout << "PostPopulation " << FullAccessionNameList[i] << "\n";
				for (j=0;j<ActiveAlleleByPopList[i].size();j++)
				{
					cout << "PostLocus " << j << ", ";
					cout << "ActiveAlleleByPopList[" << i << "][" << j << "].size()=" << ActiveAlleleByPopList[i][j].size() << ", ";
					si = ActiveAlleleByPopList[i][j];
					for (std::vector<int>::iterator it=si.begin(); it!=si.end(); ++it)
						cout << *it << ",";
					cout << "\n";
				}
			}
			//print out smallest sample size
			for (unsigned int i=0;i<sss.size();++i) cout << sss[i] << ",";
			cout << "\n";
	
		//print the *MaxAllelesList
		for (unsigned int i=0;i<ActiveMaxAllelesList.size();++i)
		{
			cout << "ActiveMaxAllelesList[" << i << "]=" << ActiveMaxAllelesList[i] << "\n";
		}
		for (unsigned int i=0;i<TargetMaxAllelesList.size();++i)
		{
			cout << "TargetMaxAllelesList[" << i << "]=" << TargetMaxAllelesList[i] << "\n";
		}
	*/	
	
	//clean up a little
	vector<unsigned int>().swap(sss);
	vector<unsigned int>().swap(sst);	

	time (&endd);
	dif = difftime (endd,startd);
	if (dif==1) cout << "    " << dif << " second.\n";	
	else cout << "    " << dif << " seconds.\n";	

	//stop the clock
	time (&endi);
	dif = difftime (endi,starti);

	if (dif == 1)
		cout << "Input files processed.  Elapsed time = "<< dif << " second.\n\n";
	else
		cout << "Input files processed.  Elapsed time = "<< dif << " seconds.\n\n";
	
		cout << "Number of accessions = " << NumberOfAccessions << ", Number of reference loci = " << NumLoci 
		<< "\n  Number of target loci = "<< TargetAlleleByPopList[1].size();
	if (Kernel == "yes") cout << ", Number of kernel accessions = " << KernelAccessionList.size() << "\n\n";
	else cout << "\n\n";
	
	//PERFORM A*
	if (Ideal == "yes")
	{
		int parallelism_enabled = 1; //0=no, not 0 = yes
		if (parallelism_enabled == 0) cout << "\nBeginning serial A* search...\n\n";
		else cout << "\nBeginning parallel A* search (" << ncpu << " threads)...\n\n";
				
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
			start1
		);
		
		//stop the clock
		time (&end1);
		dif = difftime (end1,start1);
		cout << "\nA* search complete.  Elapsed time = "<< dif << " seconds.\n\n";
	}	
	
	//PERFORM M+
	if (DoM == "yes")
	{
		//compile OpenMP pragma as parallel or not?
		int parallelism_enabled = 0; //0=no, not 0 = yes
		if (parallelism_enabled == 0) cout << "\nBeginning serial M+ search...\n\n";
		else cout << "\nBeginning parallel M+ search (" << ncpu << " threads)...\n\n";

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
		rng,
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
		double dif = difftime (endm,startm);
		cout << "\nM+ search complete.  Elapsed time = "<< dif << " seconds.\n\n";
	
	}
	
	return 0;
}
