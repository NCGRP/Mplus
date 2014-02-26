#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <math.h>
#include <numeric>
#include <omp.h>
#include <set>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>
using namespace std;

/***************STRUCTS*****************/
struct Alfreq //contains allele frequencies for all alleles at one locus
{
	std::string locusname;
	vector<std::string> allelenames;
	vector<double> frequencies;
};



/***************CLASSES*****************/
class Node
{
	public:
		//variables
		Node();
		~Node();
	

		//set functions
		void Setf0(double v);
		void Setg(double s);
		void Seth0(double t);
		void Seta0(double t);
		void Seta1(double t);
		void Seta2(double t);
	
		void SetAccName(std::string p);
		void SetAlleles(vector<vector<std::string> > q);
		void SetAlleleCounts(vector<int> r);
		void SetParent(std::string u);
		void SetPopSize(int w);
		
		//get functions
		double Getf0(), Getg(), Geth0(), Geta0(), Geta1(), Geta2();
		std::string GetAccName();
		std::string GetParent();
		vector<int> GetAlleleCounts();
		vector<vector<string> > GetSetOfAlleles();
		int GetPopSize();
		
		//special functions
		/*// declare a sort functor nested within Node, old version
		struct gSort : public binary_function<Node, Node, bool>
	   	{
		  	bool operator()(Node i, Node j)
			{
				(i.Getg() < j.Getg());
			}
		};*/
		
		
		//sort functor for SortedCostNodeList::push, ranks Nodes by distance metrics
		struct fSort : public binary_function<Node, Node, bool>
	   	{
		  	bool operator()(Node i, Node j)
			{
				if (i.Getf0() < j.Getf0()) return true; 		// estimated path length is shorter
				else if (i.Getf0() == j.Getf0()) 				//path length is the same
				{
					if (i.Geta0() > j.Geta0()) return true; 	//# alleles donated is larger
					else if (i.Geta0() == j.Geta0()) 			//# alleles the same
					{
						if (i.Geta1() < j.Geta1()) return true;	//rarest allele is more rare
						else if (i.Geta1() == j.Geta1())		//rarest alleles are same frequency
						{
							if (i.Geta2() < j.Geta2()) return true;	//average allele frequency is lower = more rare alleles
							else return false;						//average allele frequency is higher = fewer rare alleles
						}
						else return false;						//rarest allele is more common
					}
					else return false;							//# alleles donated is smaller
				}
				else return false; 								//path length longer
			}
		};
		
		
	
	protected:
		double f0, g;
		double h0; //minimum accessions remaining to collect all alleles at worst collected locus
		double a0; //# alleles gained over parent path by adding this node
		double a1; //frequency of rarest allele, favor accession with rarest allele
		double a2; //average frequency of alleles, across all loci. favor accession with lowest mean freq
		std::string AccName;
		vector<vector<string> > SetOfAlleles;
		vector<int> AlleleCounts;
		std::string Parent;
		int PopSize;

};

//contains a list of nodes for the OPEN and CLOSED list, automatically sorted by cost
class SortedCostNodeList
{
	public:
		//variables
		SortedCostNodeList();
		~SortedCostNodeList();
		vector<Node> s;
		
		//functions
		void push(Node); //place a new Node on list, sort by cost (f0), then by a0, a1, a2 to break ties
		Node pop(); //return the Node with the lowest cost
		vector<Node> Gets(); //return the sorted vector of Nodes, s
		void pushs(vector<Node>); //place a vector of updated Nodes on list, sort as in ::push
	protected:

};



/***************VARIABLES*****************/


/***************SHARED FUNCTIONS IN M+.cpp*****************/
std::vector<std::string> MyFilterDuplicates (std::vector<std::string> ListToFilter);
std::vector<std::string> MyFilterDuplicatesII (std::vector<std::string> ListToFilter);

/***************SHARED FUNCTIONS IN aStar.cpp*****************/
void printOPENList(SortedCostNodeList OPENlist);
void printCLOSEDList(SortedCostNodeList CLOSEDlist);
void printAllNodes(std::vector<Node> AllNodes);

/***************FUNCTIONS SHARED BETWEEN M+.cpp and aStar.cpp*****************/
int aStar (char* IdealFilePath, vector<vector<vector<std::string> > > ActiveAllelesByPopList, std::vector<int> ActiveMaxAllelesList, std::vector<std::string> UniqLociNamesList, std::vector<int> ReferenceOrTargetKey, vector<std::string> FullAccessionNameList, vector<int> PloidyList, vector<int> PopSizes, vector<Alfreq> AlleleFrequencies, int parallelism_enabled);
