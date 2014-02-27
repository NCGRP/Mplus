#include "m+.hpp"

//reconstruct the path from the node back to the start
vector<std::string> MyReconstructPath(Node e, vector<Node> AllNodes)
{
	std::string Parent;
	std::string AccName;
	std::string f;
	vector<std::string> path; //path is going to contain the calculated path in reverse order
	
	Parent = e.GetParent();
	
	while (Parent != "start")
	{
		//get the node associated with the parent of current node in path
		for (unsigned int i=0;i<AllNodes.size();++i)
		{
			f = AllNodes[i].GetAccName();
			if (f == Parent)
			{
				path.push_back(f); //add the matching accession name to path
				Parent = AllNodes[i].GetParent(); //get the parent of node with AccName = f
				break;
			}
		}
	}

	std::reverse(path.begin(), path.end()); //reverse the path
	return path;
}

//reconstruct the path from the node back to the start, thread safe
vector<std::string> MyReconstructPathII(Node e, vector<Node> AllNodes)
{
	std::string Parent;
	std::string AccName;
	std::string f;
	vector<std::string> path; //path is going to contain the calculated path in reverse order
	path.resize(e.Getf0() + 1);
	
	Parent = e.GetParent();
	
	int p = 0;  //p indexes the position in the parent path vector
	while (Parent != "start")
	{
		//get the node associated with the parent of current node in path
		for (unsigned int i=0;i<AllNodes.size();++i)
		{
			f = AllNodes[i].GetAccName();
			if (f == Parent)
			{
				//add the matching accession name to path
				path[p] = f;
				++p; //index++ so next parent is put in the next vector position
				
				Parent = AllNodes[i].GetParent(); //get the parent of node with AccName = f
				break;
			}
		}
	}

	std::reverse(path.begin(), path.end()); //reverse the path
	return path;
}



//calculate distances, add to node
void MyCalculatef(Node& e, vector<Node> AllNodes, vector<int> goalstate, vector<int> RefPloidyList, vector<Alfreq> AlleleFrequencies)
{
	unsigned int i, j, k, l;
	double g, h0, a0, a1, a2;
	double temp;
	
	//g = number of accessions in current path, equals size of reconstructed path
	vector<std::string> ParentPath;
	ParentPath = MyReconstructPath(e, AllNodes);
	g = double(ParentPath.size()) - 1.0; //minus 1 because path length includes node_start
	e.Setg(g);
	
	 /*
	To estimate the path length remaining assuming this node is part of the final path, calculate the heuristic h.
	1. Tabulate the difference in allele count between the currentstate and goalstate
	for each locus. This equals the number of alleles remaining to be collected at each
	locus.
	2. Adjust the number of alleles needed by dividing by the ploidy for each locus.  
	This value represents the minimum number of individuals necessary to completely
	collect the locus (assumes all loci heterozygous).  Note it may be a fraction, don't worry.
	3. Find the locus that requires the most individuals to completely collect it. Save
	the number of individuals.
	4. Determine the PopSize for all accessions (nodes) not yet in the ParentPath. Sort
	these by size.
	5. Using the worst case number of individuals needed from step 3, determine how many
	accessions (from 4) are necessary to encompass that many individuals.  This is the
	minimum number of accessions needed to complete the collection of alleles, h0.  
	*/
	
	//h = minimum # accessions required to collect all remaining alleles
	//at worst collected locus, assuming Node e is part of the solution.
	vector<int> currentstate = e.GetAlleleCounts();
		
	vector<double> d(goalstate.size());  //this vector will hold the difference in allele count
									     //between current state and goal, for each locus
	//Step 1. Get # alleles needed at each locus
	for (i=0;i<d.size();++i)
	{
		d[i] = (double)goalstate[i] - (double)currentstate[i];
		//cout << goalstate[i] <<"	"<< currentstate[i] << "\n";
	}
	//cout << std::accumulate(goalstate.begin(), goalstate.end(),0) << "	" << std::accumulate(currentstate.begin(), currentstate.end(),0)<<"\n";;
	
	//Step 2. Adjust by ploidy to get minimum number of individuals needed
	vector<double> z;
	for (i=0;i<d.size();++i)
	{
		z.push_back( d[i] / RefPloidyList[i] );	
	}
	
	//Step 3. Find the locus that will require the most individuals to get all alleles
	double mostIndiv = *max_element(z.begin(), z.end());
	
	//Step 4. Get PopSize for all available nodes, sort
	//get all accession names
	vector<std::string> AllAcc;
	Node n;
	for (i=0;i<AllNodes.size();++i) 
	{
		n = AllNodes[i];
		AllAcc.push_back(n.GetAccName());
	}
	
	//remove accessions in parent path from AllAcc, put difference in vector v
	vector<std::string> v;
	std::sort(ParentPath.begin(), ParentPath.end());
	std::sort(AllAcc.begin(), AllAcc.end());
	std::set_difference( AllAcc.begin(), AllAcc.end(), ParentPath.begin(), ParentPath.end(), std::inserter(v, v.end()) );

	/*cout << "ParentPath[i]                =";
	for (i=0;i<ParentPath.size();++i) cout << ParentPath[i] << " ";
	cout << "\n";
	cout << "AllAcc[i]                    =";
	for (i=0;i<AllAcc.size();++i) cout << AllAcc[i] << " ";
	cout << "\n";
	cout << "AllAcc[i] minus ParentPath[i]=";
	for (i=0;i<v.size();++i) cout << v[i] << " ";
	cout << "\n\n";
	*/

	//get PopSizes for remaining accessions, listed in v;
	vector<int> p;
	std::string s, t;
	for (i=0;i<v.size();++i)
	{
		s = v[i]; //get one of the remaining accessions
		for (j=0;j<AllNodes.size();++j)
		{
			t = AllNodes[j].GetAccName();
			if (t == s) //you've identified the correct node, add Node to p
			{
				p.push_back(AllNodes[j].GetPopSize());
				break;	
			}
		}
	}

	//sort list of PopSizes, largest to smallest
	std::sort(p.begin(), p.end(), std::greater<double>());
	
	//count number of accessions necessary to accumulate mostIndiv individuals
	j=0;
	for (i=0;i<p.size();++i)
	{
		j = j + p[i];
		if (j >= mostIndiv) break; //i now indicates the number of accessions needed (minus 1)
	}
	
	h0 = i + 1;
	e.Seth0(h0);
	e.Setf0(g+h0);
		
	/*
	In practice, the above heuristic h0 will seldom be meaningful as it assumes complete
	heterozygosity of all individuals in an accession.  h0 will often be the same for 
	several accessions.  For this reason, to settle ties, we need to calculate ancillary metrics 
	indicating the present node's suitability to be included.  These include metrics that
	measure the value of the node in terms of the total number of new alleles it will
	contribute (a0), the rarest allele it contains that remains to be sampled (a1, because, 
	assuming h0 and a0 are equal, the accession with the rarest allele is the best to 
	include now), and a measure of the overall rarity of the alleles it can offer to the
	collection (a2).
	*/
	
	/*calculate a0
	1. Get list of all alleles in the ParentPath excluding the current node, pc.
	2. Get list of alleles for the current node, e.GetSetOfAlleles(), called c.
	3. Calculate a0 as the number of new alleles in current node relative to parent node
	*/
	
	//1.
	//traverse the ParentPath to pick up all sets of alleles
	vector<vector<std::string> > pc(RefPloidyList.size()); //will hold only unique alleles accumulated at each locus, size to number of loci
	vector<vector<std::string> > u; //holds set of alleles, temporarily
	for (i=0;i<ParentPath.size();++i)
	{
		s = ParentPath[i];
		for (j=0;j<AllNodes.size();++j)
		{
			t = AllNodes[j].GetAccName();
			if (t == s) //you've identified the correct node
			{
				//extract the set of alleles
				u = AllNodes[j].GetSetOfAlleles();
				//add them to a master 2d list of alleles
				for (k=0;k<pc.size();++k)
				{
					pc[k].insert(pc[k].end(), u[k].begin(), u[k].end()); // add the list of alleles to the current 2d list
				}
				break; //get next node in parent path
			}
		}
	}
		
	//sort pc and remove duplicate elements
	for (i=0;i<pc.size();++i)
	{
		v = pc[i];
		std::sort(v.begin(), v.end());
		v.erase(std::unique(v.begin(), v.end()), v.end());	
		pc[i] = v;
	}
		
	//2. 
	vector<vector<string> > c = e.GetSetOfAlleles();

	//3. Get unique elements in c relative to pc
	a0 = 0;
	vector<std::string> w, x;
	u.clear(); //clear old 2d string vector
	u.resize(RefPloidyList.size()); //resize it in case the clear does not destroy the top level
	for (i=0;i<pc.size();++i)
	{
		v = pc[i]; //v = alleles at locus i for ParentPath
		w = c[i];  //c = alleles at locus i for current node
		
		//use set_difference to find the elements in w that are unique relative to v
		vector<std::string>().swap(x); //clear x
		std::sort(v.begin(), v.end());
		std::sort(w.begin(), w.end());
		std::set_difference( w.begin(), w.end(), v.begin(), v.end(), std::inserter(x, x.end()) ); 
		a0 = a0 + x.size(); //size of x is the number of alleles unique to current node relative to all nodes in parent path
		
		u[i] = x; //put alleles unique to current node into a 2d vector, one item for each locus
		
		/*//print out contributions of this accession to diversity
		cout << "all,    w: ";
		for (j=0;j<w.size();++j) cout << w[j] << " ";
		cout << "\n";
		cout << "novel,  x: ";
		for (j=0;j<x.size();++j) cout << x[j] << " ";
		cout << "\n";
		cout << "parent, v: ";
		for (j=0;j<v.size();++j) cout << v[j] << " ";
		cout << "\n\n";
		*/
	}
	e.Seta0(a0);
	
	/*calculate a1, the frequency of the rarest allele added by the current node
	1. Vector u contains a list of unique alleles added by the current node, get the
	frequency of each one of these. 
	2. Find minimum frequency, this is a1
	*/
	
	//1.
	vector<double>().swap(z); //clear z
	vector<double> zz;  //zz will hold a list of all frequencies for alleles added by current node
	vector<std::string>::iterator it;
	for (i=0;i<u.size();++i)
	{
		w = AlleleFrequencies[i].allelenames; //get the list of alleles for locus i
		z = AlleleFrequencies[i].frequencies; //get the list of allele frequencies
		
		for (j=0;j<u[i].size();++j)
		{
			s = u[i][j];//get name of one of the alleles unique to this accession, at locus i
			//find its frequency by matching in w, the list of all allele names
			it = std::find(w.begin(), w.end(), s); //search for current allele in list of all alleles
			l = std::distance(w.begin(), it); //iterator to index
			temp = z[l]; //get the frequency
			zz.push_back(temp); //add it to the list of frequencies for all loci mashed together
		}
	}
	
	//2.
	if (zz.size() > 0) a1 = *std::min_element(zz.begin(), zz.end());
	else a1 = -1.0; //add a default value if no new alleles are added by current node
	e.Seta1(a1);
	
	/*calculate a2, the average frequency of alleles added by the current node
	1. Vector zz contains the list of frequencies of allunique alleles added by the 
	current node, get the mean
	*/
	if (zz.size() > 0) 	a2 = std::accumulate(zz.begin(), zz.end(), 0.0) / zz.size();
	else a2 = -1.0;
	e.Seta2(a2);
}

//counts the number of alleles at each locus, eliminating redundancies by using set
vector<int> MyCountAllelesAtEachLocus(vector<vector<std::string> > a)
{
	int b;
	unsigned int i, j;
	vector<int> AlleleCounts;
	set<std::string> AlleleSet; //non-redundant list of alleles at a locus
	vector<std::string> t;

	for (i=0;i<a.size();++i)
	{
		t = a[i]; //get all alleles at locus i
		
		//eliminate redundancy in t by putting into a set
		AlleleSet.clear(); //clear AlleleSet
		
		//filter out redundant alleles by dumping the vector into a set
		for (j=0;j<t.size();++j) AlleleSet.insert(t[j]);
		
		b = AlleleSet.size();
		AlleleCounts.push_back(b);
	}

	return AlleleCounts;
}

vector<int> MyGetUpdatedAlleleCounts(vector<std::string> ParentPath, vector<Node> AllNodes, vector<int> RefPloidyList)
{
	unsigned int i, j, k;
	std::string s, t;
	
	vector<vector<std::string> > pc(RefPloidyList.size()); //will hold only unique alleles accumulated at each locus, size to number of loci
	vector<vector<std::string> > u; //holds set of alleles, temporarily
	vector<int> w; //holds count of alleles at each locus
	
	///traverse the ParentPath to pick up all sets of alleles
	for (i=0;i<ParentPath.size();++i)
	{
		s = ParentPath[i];
		for (j=0;j<AllNodes.size();++j)
		{
			t = AllNodes[j].GetAccName();
			if (t == s) //you've identified the correct node
			{
				//extract the set of alleles
				u = AllNodes[j].GetSetOfAlleles();
				//add them to a master 2d list of alleles
				for (k=0;k<pc.size();++k)
				{
					pc[k].insert(pc[k].end(), u[k].begin(), u[k].end()); // add the list of alleles to the current 2d list
				}
				break; //get next node in parent path
			}
		}
	}
	
	//reduce to unique alleles
	w = MyCountAllelesAtEachLocus(pc);
	 
	return w;
}

void printOPENlist(SortedCostNodeList OPENlist)
{
	//print out OPENlist
	vector<Node> ess = OPENlist.Gets();
	cout << "OPENlist\nAccName	Parent	f0	g	h0	a0	a1	a2\n";
	for (unsigned int i=0;i<ess.size();++i)
	{
		cout 	<< ess[i].GetAccName() << "\t"
				<< ess[i].GetParent() << "\t" 
				<< ess[i].Getf0() << "\t" 
				<< ess[i].Getg() << "\t" 
				<< ess[i].Geth0() << "\t" 
				<< ess[i].Geta0() << "\t" 
				<< ess[i].Geta1() << "\t" 
				<< ess[i].Geta2()<<"\n";
	}
}

void printCLOSEDlist(SortedCostNodeList CLOSEDlist)
{
	//print out CLOSEDlist
	vector<Node> ess = CLOSEDlist.Gets();
	cout << "\nCLOSEDlist\nAccName	Parent	f0	g	h0	a0	a1	a2\n";
	for (unsigned int i=0;i<ess.size();++i)
	{
		cout 	<< ess[i].GetAccName() << "\t"
				<< ess[i].GetParent() << "\t" 
				<< ess[i].Getf0() << "\t" 
				<< ess[i].Getg() << "\t" 
				<< ess[i].Geth0() << "\t" 
				<< ess[i].Geta0() << "\t" 
				<< ess[i].Geta1() << "\t" 
				<< ess[i].Geta2()<<"\n";
	}
}

void printAllNodes(vector<Node> AllNodes)
{		
	//print out AllNodes
	vector<Node> ess = AllNodes;
	cout << "\nAllNodes\nAccName	Parent	f0	g	h0	a0	a1	a2\n";
	for (unsigned int i=0;i<ess.size();++i)
	{
		cout 	<< ess[i].GetAccName() << "\t"
				<< ess[i].GetParent() << "\t" 
				<< ess[i].Getf0() << "\t" 
				<< ess[i].Getg() << "\t" 
				<< ess[i].Geth0() << "\t" 
				<< ess[i].Geta0() << "\t" 
				<< ess[i].Geta1() << "\t" 
				<< ess[i].Geta2()<<"\n";
	}
}


//functor to retrieve Node that matches AccName
struct FindNodeViaAccName : public binary_function<Node, std::string, bool>
	{
		bool operator()(Node& i, const std::string& nn) const
		{
			return i.GetAccName() == nn;
		}
	};





//class Node methods

//constructor and destructor
Node::Node() {}
Node::~Node() {}

//get info from class Node
double Node::Getf0() {return f0;}
double Node::Getg() {return g;}
double Node::Geth0() {return h0;}
double Node::Geta0() {return a0;}
double Node::Geta1() {return a1;}
double Node::Geta2() {return a2;}

std::string Node::GetParent() {return Parent;}
std::string Node::GetAccName() {return AccName;}
vector<int> Node::GetAlleleCounts() {return AlleleCounts;}
vector<vector<string> > Node::GetSetOfAlleles() {return SetOfAlleles;}
int Node::GetPopSize() {return PopSize;}

//place info into class Node
void Node::SetAccName(std::string p) {AccName = p;}
void Node::SetAlleles(vector<vector<string> > q) {SetOfAlleles = q;}
void Node::SetAlleleCounts(vector<int> r) {AlleleCounts = r;}
void Node::SetParent(std::string u) {Parent = u;}

void Node::Setf0(double v) {f0 = v;}
void Node::Setg(double s) {g = s;}
void Node::Seth0(double t) {h0 = t;}
void Node::Seta0(double t) {a0 = t;}
void Node::Seta1(double t) {a1 = t;}
void Node::Seta2(double t) {a2 = t;}

void Node::SetPopSize(int w) {PopSize = w;}



//class SortedCostNodeList methods

//constructor and destructor
SortedCostNodeList::SortedCostNodeList() {}
SortedCostNodeList::~SortedCostNodeList() {}

//add a single node to the class
void SortedCostNodeList::push(Node x)
{
	s.push_back(x); //add the new node to the vector s
	std::sort(s.begin(), s.end(), Node::fSort());//sort vector ascending
	std::reverse(s.begin(), s.end()); //reverse vector so shortest is at end
}

//add a vector of nodes to the class
void SortedCostNodeList::pushs(vector<Node> x)
{
	s = x; //move input vector to s
	std::sort(s.begin(), s.end(), Node::fSort());//sort vector ascending
	std::reverse(s.begin(), s.end()); //reverse vector so shortest is at end
}

//return least cost node from the class, remove
Node SortedCostNodeList::pop()
{
	Node r = s.back(); //get the last node off the list s, already sorted so lowest cost is last
	s.pop_back(); //remove the last node
	return r;
}

//return the sorted vector of Nodes
vector<Node> SortedCostNodeList::Gets() {return s;}



//A*
int aStar (char* IdealFilePath, vector<vector<vector<std::string> > > ActiveAlleleByPopList, vector<int> ActiveMaxAllelesList, vector<std::string> UniqLociNamesList, vector<int> ReferenceOrTargetKey, vector<std::string> FullAccessionNameList, vector<int> PloidyList, vector<int> PopSizes, vector<Alfreq> AlleleFrequencies, int parallelism_enabled)
{
	//SET UP
	unsigned int i;
	
	//get total number of alleles
	int TotAlleles = 0;
	for (i=0;i<ActiveMaxAllelesList.size();++i) TotAlleles = TotAlleles + ActiveMaxAllelesList[i];
		
	//get list of unique allele names at each locus, remove target loci from UniqLociNamesList
	// and PloidyList, reference(0), target(1)
	int b;
	vector<std::string> RefLociNamesList;
	vector<int> RefPloidyList;
	for (i=0;i<UniqLociNamesList.size();++i)
	{
		b = ReferenceOrTargetKey[i];
		if (b == 0) 
		{
			RefLociNamesList.push_back(UniqLociNamesList[i]);
			RefPloidyList.push_back(PloidyList[i]);
		}
	}
	
	//PREPARE FOR A*

	//set up OPEN and CLOSED lists
	SortedCostNodeList OPENlist;
	SortedCostNodeList CLOSEDlist;
	
	//set initial state
	vector<int> startstate (ActiveMaxAllelesList.size(),0); //all loci have zero alleles for node_start

	//set goal state, goal is a vector of the count of all alleles across all reference loci
	vector<int> goalstate = ActiveMaxAllelesList;
	
	//fill in starting state for all Nodes, this is the basic description of each accession
	//do not fill in information for Parent, AlleleCounts, and distances
	vector<Node> AllNodes;
	for (i=0;i<FullAccessionNameList.size();++i)
	{
		Node node_base;
		node_base.SetAccName(FullAccessionNameList[i]);
		node_base.SetAlleles(ActiveAlleleByPopList[i]);
		node_base.SetPopSize(PopSizes[i]);
		
		AllNodes.push_back(node_base);
	}
	
	//create start node
	//set basic info
	Node node_start;
	node_start.SetAccName("start");
	
	//load a default set of alleles (just empty lists)
	vector<vector<std::string> > foo(RefLociNamesList.size());
	node_start.SetAlleles(foo); //add empty vectors for each locus
	
	node_start.SetAlleleCounts(startstate);
	node_start.SetPopSize(1);//1 used for no special reason, for other nodes it will affect calc of distances, but not for node_start
	node_start.SetParent("start");
	
	AllNodes.push_back(node_start);
	MyCalculatef(node_start, AllNodes, goalstate, RefPloidyList, AlleleFrequencies);

	OPENlist.push(node_start); //put node_start on OPENlist



	//ENTER A* ALGORITHM
	Node node_current, node_successor;
	vector<std::string> ParentPath;
	vector<Node> successorNodes;
	vector<int> currentstate;

	vector<std::string> AccNameList;
	AccNameList = FullAccessionNameList;
	AccNameList.push_back("start"); //add start state to list of node names, for use to determine valid successor nodes
	
	std::string nn;
	vector<string>::iterator it;
	vector<Node>::iterator itn;
	vector<Node> s;
	vector<std::string> v;
	vector<int> w;
	int l;
	Node tempNode;
	static const SortedCostNodeList emptyl; //this will be used to zero OPENlist between loops
	std::string Soln; //holds whether a solution was found or not

	while (true)
	{

		//exit condition 1, OPENlist is empty so there is no solution
		if (OPENlist.s.size() == 0) 
		{
			Soln = "no";
			cout << "Sorry, there is no solution for A* optimization.\n";
			break;
		}

		//get node with lowest f0 (or a0 or a1 or a2) on OPENlist, call it node_current
		node_current = OPENlist.pop();  //pop removes the best node because everything is sorted when it is pushed onto the OPENlist
	
		//add node_current to CLOSEDlist
		CLOSEDlist.push(node_current);

		//Display status of search
		currentstate = node_current.GetAlleleCounts();
		vector<std::string>().swap(ParentPath); //clear ParentPath
		ParentPath = MyReconstructPath(node_current, AllNodes); //reconstruct path
		if (node_current.GetAccName() != "start")
		{
			cout 	<< "Adding accession: " << node_current.GetAccName()
					<< "  Current core size: " << ParentPath.size() + 1
					<< "  Alleles captured: " << std::accumulate(currentstate.begin(), currentstate.end(), 0) << "/" << std::accumulate(goalstate.begin(), goalstate.end(),0) << "\n";
		}

		//exit condition 2, current state of node_current = ActiveMaxAllelesList (which is goalstate)
		if (currentstate == goalstate) 
		{
			Soln = "yes"; //node_current contains the final accession necessary, trace back its parents to get minimum core
			//trace back the ParentPath from node current to get the ideal core set
			ParentPath = MyReconstructPath(node_current, AllNodes);
			ParentPath.push_back(node_current.GetAccName());
			
			//output ordered core to terminal
			cout << "\nThe ideal core set contains "<<ParentPath.size()<<" accessions.\n";
			cout << "One ideal core = ";
			for (i=0;i<ParentPath.size();++i) 
			{
				if (i == ParentPath.size() - 1) cout << ParentPath[i] << "\n";
				else cout << ParentPath[i] << ",";
			}
			
			//write ordered core to output file
			ofstream output; //set up file stream for output file
			output.open(IdealFilePath);
			output.close(); //quick open close done to clear any existing file each time program is run
			output.open(IdealFilePath, ios::out);
			output << "Ideal core\n(";
			for (i=0;i<ParentPath.size();++i) 
			{
				if (i == ParentPath.size() - 1) output << ParentPath[i] << ")\n";
				else output << ParentPath[i] << ",";
			}
			output.close();

			//write sorted core to terminal
			cout << "Sorted = ";
			std::sort(ParentPath.begin(), ParentPath.end());
			for (i=0;i<ParentPath.size();++i) 
			{
				if (i == ParentPath.size() - 1) cout << ParentPath[i] << "\n";
				else cout << ParentPath[i] << ",";
			}
			break;
		}

		//generate each possible node_successor that can come after node_current.
		//successors include anything not in node_current's parent path, or on CLOSEDlist
		std::sort(ParentPath.begin(), ParentPath.end());
		std::sort(AccNameList.begin(), AccNameList.end());

		//elements in AccNameList but not in ParentPath are possible successors, identify with set_difference
		vector<std::string>().swap(v); //clear v, v will receive difference in two vectors
		std::set_difference( AccNameList.begin(), AccNameList.end(), ParentPath.begin(), ParentPath.end(), std::inserter(v, v.end()) );

		//elements in v but not in CLOSEDlist (or ParentPath) are successors
		s = CLOSEDlist.Gets(); //get list of Nodes in CLOSEDlist
		for (i=0;i<s.size();++i)
		{
			nn = s[i].GetAccName(); //get the name of accession i off of CLOSEDlist
			it = std::find(v.begin(), v.end(), nn);
			if (it != v.end()) //if the item is found, i.e. the iterator returned by std::find is not last (which is returned if nothing is found)
			{
				v.erase(it); //erase the common item
			}
		}

		/*
		cout << "generating successors: ParentPath=";
		for (i=0;i<ParentPath.size();++i) cout << ParentPath[i] << ",";
		cout << "\n";
		cout << "                       successors=";
		for (i=0;i<v.size();++i) cout << v[i] << ",";
		cout << "\n";
		//getchar();
		*/
		
		//create a Node for each successor, make node_current its parent
		for (i=0;i<v.size();++i)
		{
			//find the node in AllNodes, change Parent to node_current.GetAccName()
			nn = v[i];  //get the successor node's AccName
			itn = std::find_if( AllNodes.begin(), AllNodes.end(), std::bind2nd(FindNodeViaAccName(), nn) ); //find iterator in AllNodes using custom comparator
			l = std::distance(AllNodes.begin(), itn); //convert iterator to index
		
			 //set the parent to node_current
			AllNodes[l].SetParent(node_current.GetAccName());
		}
			
		//update distances for each successor node, pass to OPENlist
		OPENlist = emptyl; //clear OPENlist so you can add updated successors
		s.clear();  //clear the vector<Node>, treated as public for input by multiple threads
		s.resize(v.size()); //resize to number of possible successor nodes
		
		#pragma omp parallel if(parallelism_enabled) 
		{
			//declare private variables
			vector<std::string> ParentPath;
			vector<int> w;
			string nn;
			vector<Node>::iterator itn;
			int l;
			int vsize = v.size();  //unsigned into coerced to signed int so compiler doesn't squawk for using unsigned int for #pragma omp for iterator

			#pragma omp for
			for (int i=0;i<vsize;++i)
			{

				//find the node in AllNodes
				nn = v[i];  //get the successor node's AccName
				itn = std::find_if( AllNodes.begin(), AllNodes.end(), std::bind2nd(FindNodeViaAccName(), nn) ); //find iterator in AllNodes using custom comparator
				l = std::distance(AllNodes.begin(), itn); //convert iterator to index

				//set AlleleCounts to unique alleles at all loci in ParentPath of node_current + this successor
				ParentPath.clear();
				ParentPath.resize(node_current.Getf0() + 1);
				ParentPath = MyReconstructPathII(AllNodes[l], AllNodes);

				/*
				cout << "AlleleCounts for " << nn << " before updating: ";
				w = AllNodes[l].GetAlleleCounts();
				for (int j=0;j<w.size();++j) cout << w[j] << ",";
				cout << "\n";
				*/
		
				ParentPath.push_back(nn);//adds this successor to ParentPath

				/*
				cout << "ParentPath for successor "<<AllNodes[l].GetAccName() << ":";
				for (int j=0;j<ParentPath.size();++j) cout << ParentPath[j] << ",";
				cout << "\n";
				*/

				w = MyGetUpdatedAlleleCounts(ParentPath, AllNodes, RefPloidyList); //calculate allele counts for ParentPath + this successor
				AllNodes[l].SetAlleleCounts(w); //set new prospective allele counts

				/*
				cout << "AlleleCounts for " << nn << " after updating: ";
				w = AllNodes[l].GetAlleleCounts();
				for (int j=0;j<w.size();++j) cout << w[j] << ",";
				cout << "\n";
				*/

				//recalculate distances for node given new parent
				MyCalculatef(AllNodes[l], AllNodes, goalstate, RefPloidyList, AlleleFrequencies); //AllNodes[l] updated as reference

				//accumulate updated nodes into a single vector, to be pushed onto OPENlist later
				s[i] = AllNodes[l];

				}
			}//end pragma omp parallel if

			//add vector of updated nodes onto OPENlist via .pushs, which sorts them properly
			OPENlist.pushs(s);
		}
	return 0;
}
