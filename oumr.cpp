
#include <numeric>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iostream>
#include <ostream>
#include <iterator>
#include <cassert>
#include <boost/program_options.hpp>
//#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <cstring>

#include <map>
#include <cstdlib>

#include <limits>
#include <memory>
#include "hmm.h"

using namespace std;


namespace po = boost::program_options;
po::variables_map options;

class Opts //Options
{
public:
	int threads;
	string methCallFile;
	string outputFile;
	string hmmModelIn;
	string hmmModelOut;
	string outputDir;
	string webOutputDir;
	string sampleName;
	string genome;
	string reference;
	int cytosineMinScore;
	int nextBaseMinScore;

	//this option is absorbed by fullMode.
	//int reportSkippedBase;
	int qualityScoreBase;
	int trimWGBSEndRepairPE2Seq;
	int trimWGBSEndRepairPE1Seq;
	int processPEOverlapSeq;
	int trimRRBSEndRepairSeq;
	int skipRandomChrom; 	//skip random chrom?
	int requiredFlag;
	int excludedFlag;
	int minFragSize;
	int minMMFragSize;
	int reportCpX; 			//'G': output CG methy file
	int reportCHX; 			//'G': output CHG methy file
	int fullMode;			//if off, only *.G.bed and *_stat.txt are generated; if on, *.bed, *_skip.bed, and *_strand.bed are generated.
	int statsOnly;
	int verbose;
	int reportStrand;		//output *_strand.bed? //this file has been incorporated in *bam.bed file and is only for debuging purpose now.
	int reportSkip; 		//output *_skip.bed?
	int reportCombined; 	//output *.bed?
	double fdr;
	int minDepthForAllC;
	int ratioUpLimit;
	Opts(){
		threads = 1;
		outputDir = "./";
		webOutputDir = "./";
		sampleName = "mSuite";
		genome = "";
		reference = "";
		methCallFile = "";

		skipRandomChrom = 1;
		requiredFlag = 0;
		excludedFlag = 0;
		reportCpX = 'G';
		reportCHX = 'X';
		fullMode = 0;
		statsOnly = 0;
		verbose = 0;
		reportStrand = 0;
		reportSkip = 0;
		reportCombined = 0;
		fdr = 0.05;
		minDepthForAllC = 5;
		ratioUpLimit = 10;
	};
} opts;


int parse_options(int ac, char * av[]){
	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", 								"Produce help message. Common options are provided with single letter format. Parameter defaults are in brackts. Example command: hmr -m sample.G.bed; See doc for more details.")
	("methCallFile,m",						po::value<string>(), "Specify the name of file from methylation calling;")
	("outputFile,o",						po::value<string>(), "Specify the name of file to write hypomethylated regions;")
	("hmmModelOut,d",						po::value<string>(), "Specify the name of file to write HMM parameters;")
	("hmmModelIn,e",						po::value<string>(), "Specify the name of file to input HMM parameters;")
	("sampleName", 							po::value<string>(), "If two or more mappedFiles are specifed, this option generates a merged result; Ignored for one input file;")
	("outputDir",							po::value<string>(), "The name of the output directory;")
	("webOutputDir",						po::value<string>(), "The name of the web-accessible output directory for UCSC Genome Browser tracks;")
	("genome,g", 							po::value<string>(), "The UCSC Genome Browser identifier of source genome assembly; mm9 for example;")
	("reference,r", 						po::value<string>(), "Reference DNA fasta file; It's required if CHG methylation is wanted;")
	("skipRandomChrom", 					po::value<int>()->default_value(1), "Specify whether to skip random and hadrop chrom;")
	("reportCpX", 							po::value<char>()->default_value('G'), "X=G generates a file for CpG methylation; A/C/T generates file for CpA/CpC/CpT meth;")
	("reportCHX", 							po::value<char>()->default_value('X'), "X=G generates a file for CHG methylation; A/C/T generates file for CHA/CHC/CHT meth; This file is large;")
	("verbose,v", 							po::value<int>()->default_value(0), "Specify whether to enter verbose mode;")
	("fdr,f", 								po::value<double>()->default_value(0.05), "Specify the FDR cutoff;")
	("minDepthForAllC,n", 					po::value<int>()->default_value(5), "Specify the minimum depth for all loci;")
	("ratioUpLimit,u", 						po::value<int>()->default_value(10), "Specify the methylation ratio up limit for a hypomethylated region;")
	("threads,p",							po::value<int>()->default_value(1),"Number of threads on all mapped file. Suggest 1~8 on EACH input file depending RAM size and disk speed.")
	;
	//toadd
	//minDepthForSingleSampleProfiling 10//use this for various sample profiling and etc.

	po::store(po::parse_command_line(ac, av, desc), options);
	po::notify(options);

	//////////////////////////////////
	if (options.count("help") || ac == 1) {
		cout << desc << endl;
		exit(1);
	}


	//////////////////////////////////
	cout<<"Options are saved in file run.config and printed here:"<<endl;
	ofstream configFile;
	//int cur_pid = pid_t ;
	stringstream configFileName;
	configFileName << "run.config." << getpid();
	configFile.open(configFileName.str().c_str(), ios_base::app);

	for(std::map<string,po::variable_value>::iterator iter = options.begin(); iter != options.end(); ++iter)
	{
		string k =  (*iter).first;

		cout		<<k<<"=";
		configFile 	<<k<<"=";


		if( k == "threads"){
			opts.threads 	= 	options[k].as<int>();
			configFile 		<<	options[k].as<int>();
			cout 			<<	options[k].as<int>();
		}

		else if( k == "webOutputDir"){
			opts.webOutputDir 	= 	options[k].as<string>();
			configFile 			<<	options[k].as<string>();
			cout 				<<	options[k].as<string>();
		}
		else if( k == "outputDir"){
			opts.outputDir 	= 	options[k].as<string>();
			configFile 		<<	options[k].as<string>();
			cout 			<<	options[k].as<string>();
		}
		else if( k== "sampleName"){
			opts.sampleName 	= 	options[k].as<string>();
			configFile 			<<	options[k].as<string>();
			cout 				<<	options[k].as<string>();
		}
		else if( k== "methCallFile"){
					opts.methCallFile 	= 	options[k].as<string>();
					configFile 			<<	options[k].as<string>();
					cout 				<<	options[k].as<string>();
		}
		else if( k== "outputFile"){
					opts.outputFile 	= 	options[k].as<string>();
					configFile 			<<	options[k].as<string>();
					cout 				<<	options[k].as<string>();
		}
		else if( k== "hmmModelIn"){
					opts.hmmModelIn 	= 	options[k].as<string>();
					configFile 			<<	options[k].as<string>();
					cout 				<<	options[k].as<string>();
		}
		else if( k== "hmmModelOut"){
					opts.hmmModelOut 	= 	options[k].as<string>();
					configFile 			<<	options[k].as<string>();
					cout 				<<	options[k].as<string>();
		}

		else if( k == "genome"){
			opts.genome 	= 	options[k].as<string>();
			configFile 		<<	options[k].as<string>();
			cout 			<<	options[k].as<string>();
		}
		else if( k == "reference"){
			opts.reference 	= 	options[k].as<string>();
			configFile 		<<	options[k].as<string>();
			cout 			<<	options[k].as<string>();
		}

		else if( k == "skipRandomChrom"){
			opts.skipRandomChrom 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}

		else if( k == "reportCpX"){
			opts.reportCpX 	= 	options[k].as<char>();
			configFile 				<<	options[k].as<char>();
			cout 					<<	options[k].as<char>();
		}
		else if( k == "reportCHX"){
			opts.reportCHX 	= 	options[k].as<char>();
			configFile 				<<	options[k].as<char>();
			cout 					<<	options[k].as<char>();
		}
		else if( k == "fullMode"){
			opts.fullMode 			= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "statsOnly"){
			opts.statsOnly 			= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "verbose"){
			opts.verbose 			= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "fdr"){
			opts.fdr 				= 	options[k].as<double>();
			configFile 				<<	options[k].as<double>();
			cout 					<<	options[k].as<double>();
		}
		else if( k == "minDepthForAllC"){
			opts.minDepthForAllC 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "ratioUpLimit"){
			opts.ratioUpLimit 		= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else{
			cerr << "Please modify code for this new type."<<endl;
		}

		cout			<<endl;
		configFile 		<<endl;
	}

	if (options.count("hmmModelOut") == 0) {
		opts.hmmModelOut = opts.methCallFile + ".model";
	}
	if (options.count("outputFile") == 0) {
		opts.outputFile = opts.methCallFile + ".hmr";
	}

	configFile.close();
	configFile.clear();



	return 0;
};







//
//static void load_cpgs(const bool VERBOSE,
//          string cpgs_file, vector<SimpleGenomicRegion> &cpgs,
//          vector<pair<double, double> > &meth, vector<size_t> &reads) {
//  if (VERBOSE)
//    cerr << "[READING CPGS AND METH PROPS]" << endl;
//  ///////////////////////////////
//
//	int colIdForChr = 0;
//	int colIdForTotalC = 4;
//	int colIdForMethC = 5;
//	int colIdForStart = 1;
//	int colIdForStrand = 6;
//	int colIdForNext = 7;
//	cout << "if no header #chom in file, default index order is #chrom, start, end, ratio, totalC, methC, strand, nextN" << endl;
//	string chrom;
//	int start;
//	int totalC;
//	int methC;
//	char strand;
//	char next;
//
//
//	string line = "";
//	ifstream inputf(cpgs_file.c_str(), ios::in);
//	while (inputf.good()) {
//		getline(inputf, line);
//		if(line == ""){continue;}
//
//		vector <string> fields;
//		boost::split(fields, line, boost::is_any_of("\t"));
//		if( fields[0] == "#chrom" )		//if there's header line, reset the index values; otherwise just use defalut index values
//		{
//			vector<string>::iterator it;
//
//			it = find (fields.begin(), fields.end(), "totalC");
//			if( it != fields.end() ){ colIdForTotalC = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }
//
//			it = find (fields.begin(), fields.end(), "methC");
//			if( it != fields.end() ){ colIdForMethC = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }
//
//			it = find (fields.begin(), fields.end(), "#chrom");
//			if( it != fields.end() ){ colIdForChr = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }
//
//			it = find (fields.begin(), fields.end(), "start");
//			if( it != fields.end() ){ colIdForStart = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }
//
//			it = find (fields.begin(), fields.end(), "strand");
//			if( it != fields.end() ){ colIdForStrand = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }
//
//			it = find (fields.begin(), fields.end(), "next");
//			if( it != fields.end() ){ colIdForNext = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }
//
//		} else {				//content
//			chrom = fields[colIdForChr];
//			start = string_to_int(fields[colIdForStart]);
//			totalC = string_to_int(fields[colIdForTotalC]);
//			methC = string_to_int(fields[colIdForMethC]);
//			strand = fields[colIdForStrand][0];
//			next = fields[colIdForNext][0];
//			if(totalC >= opts.minDepthForAllC){
//				reads.push_back( totalC );
//				meth.push_back(std::make_pair(methC, totalC-methC));
//				stringstream ss;
//				ss << chrom << "\t" << fields[colIdForStart] << "\t" << fields[colIdForStart + 1];
//				cpgs.push_back(SimpleGenomicRegion(ss.str() ) );
//			}
//		}
//	}
//	inputf.close();
//
//
//  if (VERBOSE)
//    cerr << "TOTAL CPGS: " << cpgs.size() << endl
//         << "MEAN COVERAGE: "
//         << accumulate(reads.begin(), reads.end(), 0.0)/reads.size()
//         << endl << endl;
//}



int main(int argc, char * argv[] ) {
  try {

	parse_options(argc, argv);
	cout << "Program started" <<endl;

		string inf = opts.methCallFile;
		string ouf = opts.outputFile;

		hmm(inf, ouf);

  }
  catch (exception &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}


