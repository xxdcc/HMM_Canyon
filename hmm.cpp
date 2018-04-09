#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <cmath>
#include <iterator>
#include "float.h"
#include <iostream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h> //round() function
#include <sstream>
//#include <boost/function.hpp>
#include <boost/algorithm/string.hpp>
//#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
//#include <boost/thread.hpp>
//#include <boost/program_options.hpp>
//#include <boost/date_time.hpp>

#include <math.h> //round() function
#include <boost/math/distributions.hpp>

using namespace std;
using namespace boost;
using namespace boost::math;

#define NUM_CHROM 100      //maximum chrom numbers
#define NSTATE     2       //Number of hidden states; 0|1 = less than 0.5 | larger than 0.5 ?

#define TRANS_TABLE_ELEM_PREC 0.00001

//idea:
//Two hidden states; 0 = hypo or 1 = hyper
//Emission: for a ratio p, emission prob of (n,k) is beta_dist_pdf(p|n,k). the emmission prob from hypo state is Integrate_0^0.5[dp* pdf(p)]
//State and Transition table need be trained.

//if the probability for a state to be in hypo is greater than 95%, it is regarded as significant hypo dmc.

enum STATE_TYPE {STATE0, STATE1};

//struct for the initial differential methylation regions
typedef struct
{
	int chromIndex;
	int start;				//index of genomic coordinate of start of region
	int end;				//index of genomic coordinate of end of region
	int cxSites;			//number of cpg sites in this region.
	int *siteCountsM1;      //mC counts in L1
	int *siteCountsT1;      //tC counts in L1
//	int *siteCountsM2;	   	//mC counts in L2
//	int *siteCountsT2;	   	//tC counts in L2
	double *prob;          	//posterior probabilities of each hidden state
	int *isSignificant;    	//if the cpg is a hmc: -1 for hypo, 1 if hyper, 0 if not hypo or hyper (i.e., non-zero values are significant)
//	double *prob;          	//posterior probabilities of states
//	int *isSignificant;    	//if the cpg is a dmc: 1 if L1 small, -1 if L2 small, 0 if non-differential. ?
}DMR_REGION;

//parameters
class HMMCONFIG {
	public:
	int minDepth;			//sequencing depth for initial HMR calling
	int maxIterations;   	//maximum number of iterations for training
	int minRegionDist;     	//minimum distance between two HMRs. HMCs with distance smaller than this threshold will be merged into a region.
	double minP;           	//threshold for confidence
	int trainPct;      		// trainPct% of initial regions for training.
	int samplingSaturation; // for hmm hmr calling, precision is not that strict. 100 is good enough to set as max sequencing depth.
	int minHmcs;			//min number of hmcs to form a Hmr
	double relaxedCondition;	//relaxedCondition > 0 : allow Hmr to contain relaxedCondtion * Num_HMCs of sites at middle hmm state.
	HMMCONFIG(){
		maxIterations = 100;
		minP = 0.95;
		trainPct = 10;
		minDepth = 3;
		minRegionDist = 500;
		samplingSaturation = 100;
		minHmcs = 3;
		relaxedCondition = 0.2;
	}
};



class MultiKey {
  public:
    int  n1;
    int  k1;

    MultiKey(long N1, int K1 )
      : n1(N1), k1(K1) {}

    bool operator<(const MultiKey &right) const
    {
        if ( n1 == right.n1 ) {
        	return k1 < right.k1;
        }
        else {
            return n1 < right.n1;
        }
    }
};

map<int, map<int, STATE_TYPE> > State; //<region index, cytosine index, state number>
map<int, map<int, int> > loci; //chromIndex->startIndex => start coordinate
map<int, string> chromNameByIndex; // index -> name

DMR_REGION *regions;

//prior prability of the first bin in HMM.
double priorProb[NSTATE];

//the transition probability table
double transitionTable[NSTATE*NSTATE];

//the emmission probability table
map<int, map<MultiKey, double> > emmTable;


int regionNum;

HMMCONFIG config;

double DefRatio=0.1;

int getSiteCount(string & inf);
void getSignificantSites();
void getOptimalStates();
int buildLookupTable();
void initHMM();
void forwardBackward(double *newTransition, int pct);
double stopClimb(double *newTransition);
void outputFile(string & outputName);
void freeMem();

int string_to_int( string s){return atoi(s.c_str());}
string itos(int i)
{
    stringstream s;
    s << i;
    return s.str();
}


class cpgMeth
{
public:
	int totalC;
	int methC;
	cpgMeth()
	{
		totalC = 0;
		methC = 0;
	};
	cpgMeth(int t, int m) : totalC(t), methC(m) {};
	void out(){
		cout << totalC << "\t" << methC << endl;
	}

};

void readCompToHashSimple(string fileName, int i, int j, map <int, map <string, map<int, cpgMeth> > > & lane ){

	int chrIndex = 0;
	int startIndex = 1;
	int totalCIndex = 4;
	int methCIndex = 5;

	string line = "";
	ifstream inputf(fileName.c_str(), ios::in);

	while (inputf.good()) {
		getline(inputf, line);
		if( (line == "") || (line == "\n") ){continue;}

		vector <string> fields;
		boost::split(fields, line, boost::is_any_of("\t"));

//		if( find(line.begin(), line.begin()+1, '#') == line.begin() )
//		{
//			chrIndex = find (fields.begin(), fields.end(), "#chrom") - fields.begin();
//			startIndex = find (fields.begin(), fields.end(), "start") - fields.begin();
//			t1Index = find (fields.begin(), fields.end(), "totalC_" + itos(i-1)) - fields.begin();
//			r1Index = find (fields.begin(), fields.end(), "nominalRatio_" + itos(i-1)) - fields.begin();
//			//cout << i << "\t" << j << endl;
//			//cout << t1Index << "\t" << r1Index << endl;
//			//cout << t2Index << "\t" << r2Index << endl;
//		}
//		else{
//			string chr = fields[chrIndex];
//			int start = string_to_int(fields[startIndex]);
//			int t1 = string_to_int(fields[t1Index]);
//			double r1 = atof((fields[r1Index]).c_str());
//
//			if(t1 > config.samplingSaturation){t1=config.samplingSaturation;}
//
//			lane[i][chr][start] = cpgMeth(t1, int(t1*r1+0.5));
//		}

		if( fields[0] == "#chrom" )		//if there's header line, reset the index values; otherwise just use defalut index values
		{
			vector<string>::iterator it;

			it = find (fields.begin(), fields.end(), "totalC");
			if( it != fields.end() ){ totalCIndex = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }

			it = find (fields.begin(), fields.end(), "methC");
			if( it != fields.end() ){ methCIndex = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }

			it = find (fields.begin(), fields.end(), "#chrom");
			if( it != fields.end() ){ chrIndex = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }

			it = find (fields.begin(), fields.end(), "start");
			if( it != fields.end() ){ startIndex = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }

		} else {				//content
			string chr = fields[chrIndex];
			int start = string_to_int(fields[startIndex]);
			int totalC = string_to_int(fields[totalCIndex]);
			int methC = string_to_int(fields[methCIndex]);

			if(totalC > config.samplingSaturation){
				double r = double(methC) / totalC;
				totalC = config.samplingSaturation;
				methC = int(totalC * r + 0.5);
			}
			lane[i][chr][start] = cpgMeth(totalC, methC);
		}
	}

	inputf.close();

}


int hmm(string & inputFile, string & outputName)
{
	double newTransitionMatrix[NSTATE*NSTATE];
	double transMatrixElemPrec;
	int count=0;

//	//Initialize the configurations
//	config.maxIterations = 100;
//	config.minP = 0.95;
//	config.trainPct = 10; //25 optimal; training is about 1/25 of total initial DMRs. ~ 1 chrom ~.
//	config.minDepth = 3; //3 optimal;
//	config.minRegionDist = 500; //500 optimal
//	config.samplingSaturation = 100;
//	config.minHmcs = 3;
//	config.relaxedCondition = 0.2;

	getSiteCount(inputFile);

	printf("reading file completed. regionNum=%d\n", regionNum);

	printf("building lookup table...\n");

	//calculate the emission probability table
	buildLookupTable();

	cout << "Initializing HMM..." << endl;

	//initialize the state of the HMM
	initHMM();


	//train the transition probability table
	//todo: use supplied dmr for training; use specified transition table to calculate hidden states
	while (1)
	{
		cout << "Iteration " << count << endl;
		count++;

		if (count>config.maxIterations)
		{
			break;
		}

		forwardBackward(newTransitionMatrix, config.trainPct);

		transMatrixElemPrec = stopClimb(newTransitionMatrix);

		cout << "Iteration " << count << endl;

		if (transMatrixElemPrec<TRANS_TABLE_ELEM_PREC)
		{
			break;
		}

		memcpy(transitionTable, newTransitionMatrix, NSTATE*NSTATE*sizeof(double));
	}

	//compute the posterior probability of states for all sites
	forwardBackward(newTransitionMatrix, 100);

	//get the DMCs if prob is > 0.95 and get optimal hidden states
	getSignificantSites();
	getOptimalStates();

	outputFile(outputName);

	freeMem();

	return 0;
}


//int main()
//{
//	string inf = "xxx";
//	string ouf = "testx";
//
//	hmm(inf, ouf);
//	return 0;
//}

//read methcall file or dmc file to profile each cytosine for two samples
int getSiteCount(string & inf)
{
	int tmpStart, tmpEnd;
	int i,j,k;

	map<int, map<int, int> > selected; //chromIndex->startIndex => 1|0 : selected | not selected for consideration

	cout << "start reading" << endl;
	//string fileName = "dmc_1_v_0.txt.chr19";
	string fileName = inf;
	map <int, map <string, map<int, cpgMeth> > > lane; //lane = laneId -> chrom -> position -> (n,k)
	readCompToHashSimple(fileName, 1, 2, lane ); //read fileName into lane
	cout << "read done" << endl;

	//select sites for consideration by looping through chromIndex->startIndex

	regionNum = 0;
	i = 0;
	for( map<string, map<int, cpgMeth> >::iterator pchr = lane[1].begin(); pchr != lane[1].end(); pchr++, i++){
	//for (i=0;i<chromNum;i++){
		string chr = pchr->first;
		//i = pchr - lane[1].begin(); //error. why?
		chromNameByIndex[i] = chr;

		//read coordinates into loci by chromIndex->startIndex => start
		j = 0;
		for(map<int, cpgMeth>::iterator it = lane[1][chr].begin(); it != lane[1][chr].end(); ++it) {
			//j = it - lane[1][chr].begin(); //why compile error?
			loci[i][j] = it->first;
			j ++;
		}

		int locusSize = lane[1][chr].size();

		//set selected for each site
		//for(map<int, cpgMeth>::iterator it = lane[1][chr].begin(); it != lane[1][chr].end(); ++it) {
		for (j=0;j<locusSize;j++)
		{
			int n1 = lane[1][chr][loci[i][j]].totalC;
			int k1 = lane[1][chr][loci[i][j]].methC;
//			int n2 = lane[2][chr][loci[i][j]].totalC;
//			int k2 = lane[2][chr][loci[i][j]].methC;

			//worry about direction later //todo
			//determine the initial DMRs based on rough condition
			//(selected = 1) is selected for training, and call DMC/DMR. Neglect sites with selected = 0;
			if ( double(k1)/n1 < DefRatio  && n1 >= config.minDepth )
			{
				selected[i][j] = 1;
			}
			else
			{
				selected[i][j] = 0;
			}
			//cout << chr << "\t" << i << "\t" << j << "\t" << locus[j] << selected[i][j] << endl;
		}

		//merge close regions
		//if distance between two consective Dmcs is larger than minRegionDist, these two are still in one initial region.
		//That's because the CpG before a region is used as null.
		//Further sceening will be applied later for DMR output.
		for (j=1;j<locusSize;j++)
		{
			if ((selected[i][j-1]==1)&&(selected[i][j]==0))
			{
				tmpStart=j;
				for (tmpEnd=tmpStart+1;tmpEnd<locusSize;tmpEnd++)
				{
					if (selected[i][tmpEnd]==1)
					{
						break;
					}
				}
				if ((tmpEnd!=locusSize)&&(loci[i][tmpEnd]-loci[i][tmpStart]<=config.minRegionDist))
				{
					for (k=tmpStart;k<=tmpEnd;k++)
					{
						selected[i][k] = 1;
					}
				}
				j=tmpEnd;
			}
		}

		for (j=0;j<locusSize;j++)
		{
			if (selected[i][j]==1)
			{
				tmpStart = j;

				for (k=tmpStart+1;k<locusSize;k++)
				{
					if (selected[i][k]==0)
					{
						break;
					}
				}

				j=k;

				regionNum++;
			}
		}
	}

	//get initially identified DMRs into regions array by looping through chromIndex->startIndex
	regions = (DMR_REGION *)malloc(regionNum*sizeof(DMR_REGION));
	regionNum = 0;
	i = 0;
	for( map<string, map<int, cpgMeth> >::iterator pchr = lane[1].begin(); pchr != lane[1].end(); pchr++, i++){
	//for (i=0;i<chromNum;i++){
		string chr = pchr->first;
		int locusSize = lane[1][chr].size();

		//merge all consective sites with selected == 1
		for (j=0;j<locusSize;j++)
		{
			if (selected[i][j]==1)
			{
				tmpStart = j;

				for (k=tmpStart+1;k<locusSize;k++)
				{
					if (selected[i][k]==0)
					{
						break;
					}
				}

				j=k; //j is reset here, in addition to the for() loop

				tmpEnd = k-1;

				if (tmpStart>0)
				{
					tmpStart--;
				}

				regions[regionNum].chromIndex=i;
				regions[regionNum].start = tmpStart;
				regions[regionNum].end = tmpEnd;
				regions[regionNum].cxSites = tmpEnd-tmpStart+1;
				regions[regionNum].siteCountsM1 = (int *)malloc(regions[regionNum].cxSites*sizeof(int));
//				regions[regionNum].siteCountsM2 = (int *)malloc(regions[regionNum].cxSites*sizeof(int));
				regions[regionNum].siteCountsT1 = (int *)malloc(regions[regionNum].cxSites*sizeof(int));
//				regions[regionNum].siteCountsT2 = (int *)malloc(regions[regionNum].cxSites*sizeof(int));

				for (k=0;k<regions[regionNum].cxSites;k++)
				{
					regions[regionNum].siteCountsM1[k] = lane[1][chr][loci[i][tmpStart+k]].methC;
//					regions[regionNum].siteCountsM2[k] = lane[2][chr][loci[i][tmpStart+k]].methC;
					regions[regionNum].siteCountsT1[k] = lane[1][chr][loci[i][tmpStart+k]].totalC;
//					regions[regionNum].siteCountsT2[k] = lane[2][chr][loci[i][tmpStart+k]].totalC;
				}
				regions[regionNum].prob = (double *)malloc(regions[regionNum].cxSites*NSTATE*sizeof(double));
				regions[regionNum].isSignificant = (int *)malloc(regions[regionNum].cxSites*sizeof(int));
				regionNum++;
			}
		}
	}

	if(0){ //this is for debug

		for (i=0;i<regionNum;i++){
			printf( "%i\t%i\t%i\t%i\t%i\t%f\t%i\n",regions[i].chromIndex, regions[i].start, regions[i].end, regions[i].cxSites, *(regions[i].siteCountsM1), *(regions[i].prob), *(regions[i].isSignificant));
		}

		ofstream locFile;
		locFile.open("binsInConsideration", ios_base::out);
		//output bin file
		for (i=0;i<regionNum;i++)
		{
			string tmpChrom = chromNameByIndex[regions[i].chromIndex];

			for (j=0;j<regions[i].cxSites;j++)
			{
				int st = loci[regions[i].chromIndex][regions[i].start];
				int en = loci[regions[i].chromIndex][regions[i].end];
				locFile << tmpChrom << "\t";
				locFile << i << "\t" << j << "\t";
				locFile << selected[regions[i].chromIndex][regions[i].start+j] << "\t";
				locFile << st << "\t";
				locFile << en << "\t";
				locFile << loci[regions[i].chromIndex][regions[i].start+j] << "\t";
				locFile << regions[i].siteCountsT1[j] << "\t";
				locFile << regions[i].siteCountsM1[j] << "\t";
//				locFile << regions[i].siteCountsT2[j] << "\t";
//				locFile << regions[i].siteCountsM2[j] << "\t";


				for (k=0;k<NSTATE;k++)
				{
					//fprintf(fh, "%1.3f\t",regions[i].prob[j*NSTATE+k]);
					locFile << regions[i].prob[j*NSTATE+k] << "\t";

				}
				//fprintf(fh, "%d\n", regions[i].isSignificant[j]);
				locFile << regions[i].isSignificant[j] << endl;
			}
		}

		locFile.close();

	}

	return 0;
}

void freeMem()
{
	int i;

	for (i=0;i<regionNum;i++)
	{
		free(regions[i].siteCountsM1);
//		free(regions[i].siteCountsM2);
		free(regions[i].prob);
		free(regions[i].isSignificant);
	}

	free(regions);
}

//calculate the emission probability for each possible Multikey
int buildLookupTable()
{
	map <MultiKey, int> marked;

	for(int i=0; i<regionNum; i++)
	{
		for (int j=0;j<regions[i].cxSites;j++)
		{
			MultiKey obs(regions[i].siteCountsT1[j], regions[i].siteCountsM1[j] );
			marked[obs] = 1;
		}
	}

	//calculate the emission probability
	for (int JJ=0; JJ<config.samplingSaturation+1; JJ++){
	for (int j=0;j<=JJ;j++)
	{
		MultiKey combi(JJ, j );
		if (!marked[combi])
		{
			continue;
		}

		double a = j + 1;
		double b = JJ - j + 1;
		beta_distribution<> bdist(a, b);
		emmTable[0][combi] = cdf(bdist, DefRatio); //0.5 is the upper limit of integration. if you want to find true ratio less than 0.1, then change it to 0.1.
		emmTable[1][combi] = 1 - cdf(bdist, DefRatio);

		if(0){ //for debug
			printf("%d\t%d\t", JJ, j );
			for (int i=0;i<NSTATE;i++){
				printf( "%f\t",emmTable[i][combi]);
			}
			printf("\n");
		}
	}
	}

	//cout << "return.." << endl;
	//exit;
	return 1;
}

//initialize HMM
void initHMM()
{
	int i,j;

	priorProb[0] = 0.0;
	priorProb[1] = 1.0;

	for (i=0;i<NSTATE;i++)
	{
		for (j=0;j<NSTATE;j++)
		{
			transitionTable[i*NSTATE+j]=1.0/NSTATE;
		}
	}
}

//Forward-Backward algorithm for training transition table. step is used for selecting a subset of regions for training
void forwardBackward(double *newTransition, int pct)
{
	double *alpha, *beta, *eps;
	double sum, sum1;
	double priorSum[NSTATE];
	double transitionSum[NSTATE*NSTATE],postProbSum[NSTATE*NSTATE];
	int maxLen;
	int i,j,k,l;

	//allocate memory and initializing arrays
	memset(priorSum,0,NSTATE*sizeof(double));
	memset(transitionSum, 0, NSTATE*NSTATE*sizeof(double));
	memset(postProbSum, 0 , NSTATE*NSTATE*sizeof(double));

	maxLen = 0;
	int regNum = regionNum * pct / 100;
	for (i=0;i<regNum;i++)
	{
		if (regions[i].cxSites>maxLen)
		{
			maxLen = regions[i].cxSites;
		}
	}

	alpha = (double *)malloc(NSTATE*maxLen*sizeof(double));
	beta = (double *)malloc(NSTATE*maxLen*sizeof(double));
	eps = (double *)malloc(NSTATE*NSTATE*maxLen*sizeof(double));

	for (l=0;l<regNum;l++)
	{
		//Forward algorithm

		MultiKey obs(regions[l].siteCountsT1[0], regions[l].siteCountsM1[0]);
		sum=0;
		for (k=0;k<NSTATE;k++)
		{
			//alpha[k] = emmissionProbTable[k][ regions[l].siteCountsM1[0] ][ regions[l].siteCountsM2[0] ]*priorProb[k];
			alpha[k] = emmTable[k][ obs ]*priorProb[k];
			//cout << "region=" << l << " state=" << k << " alpha[k]=" << alpha[k] ;
			//cout << " emmTable[k][ obs ]*priorProb[k]=" << emmTable[k][ obs ] << " * " << priorProb[k] << endl;
		}

		for (i=1;i<regions[l].cxSites;i++)
		{
			sum1 = 0;

			for (j=0;j<NSTATE;j++)
			{
				sum = 0;

				for (k=0;k<NSTATE;k++)
				{
					sum+=alpha[(i-1)*NSTATE+k]*transitionTable[k*NSTATE+j];
				}

				//alpha[i*NSTATE+j] = sum*emmissionProbTable[j][ regions[l].siteCountsM1[i] ][ regions[l].siteCountsM2[i] ];
				MultiKey obs1(regions[l].siteCountsT1[i], regions[l].siteCountsM1[i]);
				alpha[i*NSTATE+j] = sum*emmTable[j][obs1];
				sum1+=alpha[i*NSTATE+j];
				//cout << "region=" << l << " state=" << i << "," << j << " alpha[i*NSTATE+j]=" << alpha[i*NSTATE+j] << endl;
			}

			//cout << "scale:" << endl;
			for (j=0;j<NSTATE;j++)
			{
				alpha[i*NSTATE+j]/=sum1;
				//cout << "region=" << l << " state=" << i << "," << j << " alpha[i*NSTATE+j]=" << alpha[i*NSTATE+j] << endl;
			}
		}

		//Backward algorithm

		for (j=0;j<NSTATE;j++)
		{
			beta[(regions[l].cxSites-1)*NSTATE+j] = 1.0/NSTATE;
		}

		for (i=regions[l].cxSites-2;i>=0;i--)
		{
			sum1 = 0;

			for (j=0;j<NSTATE;j++)
			{
				sum = 0;

				for (k=0;k<NSTATE;k++)
				{
					//sum+=beta[(i+1)*NSTATE+k]*transitionTable[j*NSTATE+k]*emmissionProbTable[k][ regions[l].siteCountsM1[i+1] ][ regions[l].siteCountsM2[i+1] ];
					MultiKey obs1(regions[l].siteCountsT1[i+1], regions[l].siteCountsM1[i+1]);
					sum+=beta[(i+1)*NSTATE+k]*transitionTable[j*NSTATE+k]*emmTable[k][obs1];
				}

				beta[i*NSTATE+j]=sum;
				sum1+=beta[i*NSTATE+j];
			}

			for (j=0;j<NSTATE;j++)
			{
				beta[i*NSTATE+j]/=sum1;
			}
		}

		//compute posterior probability

		for (i=0;i<regions[l].cxSites;i++)
		{
			sum1=0;
			for (j=0;j<NSTATE;j++)
			{
				regions[l].prob[i*NSTATE+j]=alpha[i*NSTATE+j]*beta[i*NSTATE+j];
				sum1+=regions[l].prob[i*NSTATE+j];
			}

			for (j=0;j<NSTATE;j++)
			{
				regions[l].prob[i*NSTATE+j]/=sum1;
			}
		}

		//compute eps

		for (i=0;i<regions[l].cxSites-1;i++)
		{
			sum1=0;
			for (j=0;j<NSTATE;j++)
			{
				for (k=0;k<NSTATE;k++)
				{
					//eps[i*NSTATE*NSTATE+j*NSTATE+k]=alpha[i*NSTATE+j]*beta[(i+1)*NSTATE+k]*emmissionProbTable[k][ regions[l].siteCountsM1[i+1] ][ regions[l].siteCountsM2[i+1] ]*transitionTable[j*NSTATE+k];
					MultiKey obs1(regions[l].siteCountsT1[i+1], regions[l].siteCountsM1[i+1]);
					eps[i*NSTATE*NSTATE+j*NSTATE+k]=alpha[i*NSTATE+j]*beta[(i+1)*NSTATE+k]*emmTable[k][obs1]*transitionTable[j*NSTATE+k];
					sum1+=eps[i*NSTATE*NSTATE+j*NSTATE+k];
				}
			}

			for (j=0;j<NSTATE;j++)
			{
				for (k=0;k<NSTATE;k++)
				{
					eps[i*NSTATE*NSTATE+j*NSTATE+k]/=sum1;

				}
			}
		}

		//update the cumulation for caluculating transition probability table

		for (i=0;i<regions[l].cxSites-1;i++)
		{
			for (j=0;j<NSTATE;j++)
			{
				for (k=0;k<NSTATE;k++)
				{
					transitionSum[j*NSTATE+k]+=eps[i*NSTATE*NSTATE+j*NSTATE+k];
					postProbSum[j*NSTATE+k]+=regions[l].prob[i*NSTATE+j];
				}
			}
		}

	}

	//update transition table
	for (j=0;j<NSTATE;j++)
	{
		for (k=0;k<NSTATE;k++)
		{
			newTransition[j*NSTATE+k] = transitionSum[j*NSTATE+k]/postProbSum[j*NSTATE+k];
		}
	}

	free(alpha);
	free(beta);
	free(eps);
}

//calculate the biggest difference for transition matrix elements
double stopClimb(double *newTransition)
{
	int i,j;
	double maxDiff = 0.0;

	for (i=0;i<NSTATE;i++)
	{

		for (j=0;j<NSTATE;j++)
		{
			double diff = abs(newTransition[i*NSTATE+j] - transitionTable[i*NSTATE+j]);
			if(diff > maxDiff){
				maxDiff = diff;
			}
		}
	}

	return maxDiff;
}

//output files
void outputFile(string & outputName)
{

	int i,j,k;
	string tmpChrom;
	int start,end;

	//get file names
	string locFileName = outputName + ".loc";
	string hmmFileName = outputName + ".hmm";
	string regionFileName = outputName + ".region";
	string dmrFileName = outputName + ".dmr";

	ofstream locFile;
	locFile.open(locFileName.c_str(), ios_base::out);


	//output loc file
	//loc header: chrom, start, end, cytosine start, cytosine end, totalC, methC, prob, emission[0], emission[1], State, isSig
	for (i=0;i<regionNum;i++)
	{
		tmpChrom = chromNameByIndex[regions[i].chromIndex];

		for (j=0;j<regions[i].cxSites;j++)
		{
			int st = loci[regions[i].chromIndex][regions[i].start];
			int en = loci[regions[i].chromIndex][regions[i].end];
			locFile << tmpChrom << "\t";
			locFile << st << "\t";
			locFile << en << "\t";
			locFile << loci[regions[i].chromIndex][regions[i].start+j] << "\t";
			locFile << loci[regions[i].chromIndex][regions[i].start+j] + 2 << "\t";
			locFile << regions[i].siteCountsT1[j] << "\t";
			locFile << regions[i].siteCountsM1[j] << "\t";
//			locFile << regions[i].siteCountsT2[j] << "\t";
//			locFile << regions[i].siteCountsM2[j] << "\t";


			for (k=0;k<NSTATE;k++)
			{
				locFile << regions[i].prob[j*NSTATE+k] << "\t";
			}

			MultiKey combi(regions[i].siteCountsT1[j], regions[i].siteCountsM1[j] );
			for (k=0;k<NSTATE;k++)
			{
				locFile << emmTable[k][combi] << "\t";
			}
			locFile << State[i][j] << "\t";
			locFile << regions[i].isSignificant[j] << endl;
		}
	}
	locFile.close();

	//output regions formed by DMCs only
	//region header: chrom, start, end, isSig, hmcCounts, state
	ofstream regionFile;
	regionFile.open(regionFileName.c_str(), ios_base::out);
	regionFile << "#HMR_chrom\tstart\tend\thypo.-1.hyper.1\thmcCounts" << endl;
	for (i=0;i<regionNum;i++)
	{
		tmpChrom = chromNameByIndex[regions[i].chromIndex];

		for (j=0;j<regions[i].cxSites;j++)
		{
			if (regions[i].isSignificant[j] != 0)
			{
				start = j;
			}
			else
			{
				continue;
			}

			for (k=start+1;k<regions[i].cxSites;k++)
			{
				if(loci[regions[i].chromIndex][regions[i].start + k] - loci[regions[i].chromIndex][regions[i].start + k - 1] > config.minRegionDist){
					break;
				}
				if (!(regions[i].isSignificant[k]==regions[i].isSignificant[k-1]))
				{
					break;
				}
			}

			end = k - 1;
			//buggy for MOABS DMR calling here
			int hmcCount = end - start + 1;
			if(hmcCount >= config.minHmcs){
				regionFile << tmpChrom << "\t";
				regionFile << loci[regions[i].chromIndex][regions[i].start+start] << "\t";
				regionFile << loci[regions[i].chromIndex][regions[i].start+end] + 2 << "\t";
				//regionFile << (regions[i].isSignificant[j]==1?"hypo":"hype") << endl;
				regionFile << regions[i].isSignificant[j] << "\t";
				regionFile << hmcCount << endl;
			}


			j=k-1;
		}
	}

	regionFile.close();


	//output DMRs by joining consective differential states, allowing 20% sites (around 1 or 2 cpgs) as middle states
	//dmr header: chr, start, end, isSig, dist, dmcCount, cpgCount ,State
	ofstream dmrFile;
	dmrFile.open(dmrFileName.c_str(), ios_base::out);

	for (i=0;i<regionNum;i++)
	{
		tmpChrom = chromNameByIndex[regions[i].chromIndex];

		for (j=0;j<regions[i].cxSites;j++)
		{

			if (regions[i].isSignificant[j] !=0)
			{
				start = j;
			}
			else
			{
				continue;
			}
			int dmcCount = 1;
			int cpgCount = 0;
			STATE_TYPE type = State[i][j];
			int lastDmc = j;
			for (k=start+1;k<regions[i].cxSites;k++)
			{
				if(State[i][k]==type){
					if(loci[regions[i].chromIndex][regions[i].start + k] - loci[regions[i].chromIndex][regions[i].start + lastDmc] > config.minRegionDist){
						cpgCount -= k - lastDmc - 1;
						k = lastDmc + 1;
						break;
					}
					if(regions[i].isSignificant[k] !=0){
						lastDmc = k;
						dmcCount ++;
					} else {
						cpgCount ++;
						if( config.relaxedCondition > 0 && cpgCount >int( dmcCount * config.relaxedCondition + 1 ) ){
							cpgCount -= k - lastDmc;
							k = lastDmc + 1;
							break;
						}
					}
				} else if(State[i][k]!=type){
					break;
				}

			}
			if(k== regions[i].cxSites ){
				cpgCount -= k - lastDmc - 1;
				k = lastDmc + 1;
			}

			end = k - 1;
			//dmrFile << end << endl;
			if(dmcCount >=config.minHmcs){ //3 hmc and 500bp seem to be good parameters
				int firstLoc = loci[regions[i].chromIndex][regions[i].start+start];
				int lastLoc = loci[regions[i].chromIndex][regions[i].start+end];

				int count = 0;
				double sum = 0.0;
				for(int x = start; x <= end; x++){
					count ++;
					sum += double(regions[i].siteCountsM1[x]) / regions[i].siteCountsT1[x];
				}
				if(count != (end -start + 1)){cout << "ERROR HERE" << endl;}
				double mean = sum / count;
				double cpgsPerKb = count / double(lastLoc-firstLoc) * 1000;
//				if(ratio * 100 < opts.ratioUpLimit && cpgsPerKb >= 5.0 ){
					dmrFile << tmpChrom << "\t";
					dmrFile << firstLoc << "\t";
					dmrFile << lastLoc + 2 << "\t";
					//dmrFile << (regions[i].isSignificant[j]==1?"hypo":"hype") << "\t";
					dmrFile << regions[i].isSignificant[j] << "\t";
					dmrFile << end -start + 1 << "\t" << dmcCount << "\t" << cpgCount << "\t";
					dmrFile << State[i][j] << "\t";
					dmrFile << mean << "\t";
					dmrFile << cpgsPerKb << endl;
//				}
			}

			j=k-1;
		}
	}

	dmrFile.close();


	//output the transition probability table
	ofstream hmmFile;
	hmmFile.open(hmmFileName.c_str(), ios_base::out);
	hmmFile << "transition matrix:" << endl;
	for (i=0;i<NSTATE;i++)
	{
		for (j=0;j<NSTATE;j++)
		{
			hmmFile << transitionTable[i*NSTATE+j] << "\t";
		}
		hmmFile << endl;
	}
	hmmFile.close();
}

// get DMCs if prob > minP(0.95)
void getSignificantSites()
{
	int i,j;

	for (i=0;i<regionNum;i++)
	{
		for (j=0;j<regions[i].cxSites;j++)
		{

			if (regions[i].prob[j*NSTATE]>config.minP && regions[i].siteCountsT1[j] >= config.minDepth)
			{
				regions[i].isSignificant[j] = -1;
			}
			else if (regions[i].prob[j*NSTATE+1]>config.minP && regions[i].siteCountsT1[j] >= config.minDepth)
			{
				regions[i].isSignificant[j] = 1;
			}
			else
			{
				regions[i].isSignificant[j] = 0;
			}
		}
	}
}

//get optimal state sequence
void getOptimalStates()
{
	int i,j;
	double sum;

	for (i=0;i<regionNum;i++)
	{
		for (j=0;j<regions[i].cxSites;j++)
		{
			sum = 0;

			double s0 = regions[i].prob[j*NSTATE];
			double s1 = regions[i].prob[j*NSTATE+1];
			if(s0 > s1){
				State[i][j] = STATE0;
			} else {
				State[i][j] = STATE1;
			}

			if (regions[i].prob[j*NSTATE]>config.minP && regions[i].siteCountsT1[j] >= config.minDepth)
			{
				regions[i].isSignificant[j] = -1;
			}
			else if (regions[i].prob[j*NSTATE+1]>config.minP && regions[i].siteCountsT1[j] >= config.minDepth)
			{
				regions[i].isSignificant[j] = 1;
			}
			else
			{
				regions[i].isSignificant[j] = 0;
			}
		}
	}
}
