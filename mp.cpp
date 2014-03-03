#include "m+.hpp"
	
int MyCalculateDiversity(vector<vector<vector<std::string> > > AlleleList, vector<int> ActiveMaxAllelesList, std::string Standardize, double& RandomActiveDiversity, double& AltRandomActiveDiversity)
{
	/*AlleleList structure:
		  Pop1..r
			  locusarray1..n		*/
	int CoreSize = AlleleList.size();
	int NumLoci = AlleleList[0].size();
	int i, j, M;
	unsigned int k;
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
		//use set to eliminate redundancies
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
		for (k=0;k<Mlist.size();++k)
		{
			M = M + Mlist[k];
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
	for (unsigned int i=0;i<TempListStr.size();++i)
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
   		 }
    	else if( i == (percent/2)){
     		bar.replace(i,1,">");
    	}
    	else{
      		bar.replace(i,1," ");
   		 }
	}

  	std::cout<< "\r" "[" << bar << "] ";
  	std::cout.width( 3 );
  	std::cout<< percent << "% complete     thread "<<omp_get_thread_num() <<" "<< std::flush;
}

//M+	
void mp(
	int MinCoreSize,
	int MaxCoreSize,
	int SamplingFreq,
	int NumReplicates,
	char* OutFilePath,
	std::string Kernel,
	vector<int> KernelAccessionIndex,
	vector<int> AccessionNameList,
	vector<vector<vector<std::string> > > ActiveAlleleByPopList,
	vector<vector<vector<std::string> > > TargetAlleleByPopList,
	vector<int> ActiveMaxAllelesList,
	vector<int> TargetMaxAllelesList,
	vector<std::string> FullAccessionNameList,
	int parallelism_enabled
	)	
{
	
	//set up variables for monitoring progress
	int percent; //percent of analysis completed
	int progindex = 0;  //index to monitor progress, percent = 100*(progindex/l)
	//below is a stupid way to calculate the number of rows in the output file, value l (which = V1) 
	//used to monitor progress and as the maximum vector index for shared output vectors
	int l=0;
	for (int i=MinCoreSize;i<MaxCoreSize+1;i=i+SamplingFreq)
	{
	for (int j=0;j<NumReplicates;j++)
		{
			l++;
		}
	}
	
	//set up vectors to fill with results
	double V1 = (double)l; //(MaxCoreSize - MinCoreSize + 1)*NumReplicates; //number of rows in output vectors
	vector<vector<double> > Results(V1, vector<double>(9)); //will contain numerical results
	vector<vector<string> > Members(V1); //will contain core set members
	

	#pragma omp parallel if(parallelism_enabled) 
	{		
		unsigned int r; //r = core size, 
		int nr, RandAcc, b, row, bsc, plateau; //nr = controller to repeat NumReplicates times
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
		RecoveryFile << "core size	random reference diversity	optimized reference diversity	random target diversity	optimized target diversity	alt random reference diversity	alt optimized reference diversity	alt random target diversity	alt optimized target diversity	core members" << "\n";
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
				//develop random starting core set
				//clear AccessionsInCore and set size
				AccessionsInCore.clear();
				AccessionsInCore.resize(r);
				
				//add kernel accessions to core, if necessary
				if (Kernel == "yes")
				{
					for (unsigned int i=0;i<KernelAccessionIndex.size();i++)
					{
						AccessionsInCore[i] = KernelAccessionIndex[i];
					}
				}

				//clear TempList and set size					
				TempList.clear();
				TempList.resize( AccessionNameList.size() );
				
				//set list of available accessions in TempList, by erasing those already in the core
				TempList = AccessionNameList;
				//expunge the kernel accessions, so they are not available for random addition below
				//KernelAccessionIndex has been reverse sorted so you don't go outside range after automatic resize by .erase
				for (unsigned int i=0;i<KernelAccessionIndex.size();i++)
				{
					b = KernelAccessionIndex[i];
					TempList.erase(TempList.begin()+b);
				}
			
				//randomly add accessions until r accessions are in the core. if there is a kernel, include those (done above)
				//plus additional, randomly selected accessions, until you get r accessions
				//for (int i=0;i<r;i++)
				for (unsigned int i=KernelAccessionIndex.size();i<r;i++)
				{
					//choose an accession randomly from those available
					RandAcc = rand() % TempList.size();
					//add it to the list
					AccessionsInCore[i] = TempList[RandAcc];
				
					//remove it from the list of available accessions
					TempList.erase(TempList.begin()+RandAcc);
				}
		
				//assemble genotypes for random core and calculate diversity
				//1. put together initial list of active alleles
				CoreAlleles.clear();
				CoreAlleles.resize( AccessionsInCore.size() );
				for (unsigned int i=0;i<AccessionsInCore.size();i++)
				{
					b = AccessionsInCore[i];
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

				//3. calculate diversity from random selection at target loci
				AlleleList.clear();
				AlleleList.resize( AccessionsInCore.size() );
				for (unsigned int j=0;j<AccessionsInCore.size();j++)
				{
					b = AccessionsInCore[j];
					AlleleList[j] = TargetAlleleByPopList[b];
				}
				MyCalculateDiversity(AlleleList, TargetMaxAllelesList, Standardize, RandomTargetDiversity, AltRandomTargetDiversity);

	
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
					for (unsigned int i=0;i<r;i++)
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
						for (unsigned int i=0;i<r;i++)
						{
							b = AccessionsInCore[i];
							CoreAlleles[i] = ActiveAlleleByPopList[b];
						}
			
						//2. go through all possible subsets of size r-1, one at a time, noting which is best.
						//If there is a kernel, do not swap out any of those accessions (they are retained as the
						//first KernelAccessionIndex.size() items in CoreAlleles).  Accomplished by starting for loop
						//at KernelAccessionIndex.size().
						best=0;
						for (unsigned int i=KernelAccessionIndex.size();i<CoreAlleles.size();i++)
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
						for (unsigned int k=0;k<BestSubCoreRevSorted.size();k++)
						{
							bsc = BestSubCoreRevSorted[k];
							TempList.erase( TempList.begin() + bsc );
						}
				
						//shuffle the list of remaining accessions, so addition order is not predictable
						std::random_shuffle (TempList.begin(), TempList.end());
					
						//add each remaining accession consecutively, calculate diversity, test 
						//whether it is better than the prior one
						best = 0;
						for (unsigned int k=0;k<TempList.size();k++)
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
					
					} //while(true) endless loop
				}
	
				//7. Calculate diversity at target loci
				//assemble the target loci allelelist for the accessions in the best core
				AlleleList.clear();
				AlleleList.resize( AccessionsInCore.size() );
				for (unsigned int j=0;j<AccessionsInCore.size();j++)
				{
					b = AccessionsInCore[j];
					AlleleList[j] = TargetAlleleByPopList[b];
				}
		
		
				//calculate diversity at target loci based upon the optimized core selection
				MyCalculateDiversity(AlleleList, TargetMaxAllelesList, Standardize, OptimizedTargetDiversity, AltOptimizedTargetDiversity);

		
				//8. Assemble stats for optimized core and add to output vectors
				//create a list of accession names from the list of accession ID's in AccessionsInCore
				sort( AccessionsInCore.begin(), AccessionsInCore.end() );
				
				TempListStr.clear();
				TempListStr.resize(r);
				for (unsigned int i=0;i<AccessionsInCore.size();i++)
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
			} //end #pragma omp for loop
	} //end #pragma omp parallel	
	
	//set up file stream for output file
	ofstream output; 
	output.open(OutFilePath);
	output.close(); //quick open close done to clear any existing file each time program is run
	output.open(OutFilePath, ios::out | ios::app); //open file in append mode
	output << "core size	random reference diversity	optimized reference diversity	random target diversity	optimized target diversity	alt random reference diversity	alt optimized reference diversity	alt random target diversity	alt optimized target diversity	core members" << "\n";
		
	//write out results row by row
	for (int i=0;i<V1;i++)
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
		for (unsigned int j=0;j<Members[i].size();j++)
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
}