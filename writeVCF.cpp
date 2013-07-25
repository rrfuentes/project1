/*
 *Last Update: July 19, 2013
 *Author: Roven Rommel B. Fuentes
 *TT-Chang Genetic Resources Center, International Rice Research Institute
 *
 *Load vcf data to HDF5.
*/

#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <map>
#include <assert.h>
#include <sys/resource.h>
#include <vector>

using namespace std;

static const char *options="f:F:p:P:n:N:i:I:c:C:r:R:";

static string datafile;
static string inputfile;
static string varName;
static string varPath;
static int col=0;
static int row=0;

/*fileformat,filter,alt,assembly,pedigree,sample,unifiedgenotyper,INFO,FORMAT*/
typedef struct{ 
    char field[20];
    char* value;  
}META_1; 

/*contigs*/
typedef struct{ 
    char* id;
    int len;
}META_2;

/*INFO&FORMAT IDs*/
typedef struct{
    char* info;
    char* format;
}INFO_FORMAT;

typedef struct{
    int chrom; /*index of META_3 or contig table*/
    int pos;
    char id[15];
    char ref;
    char* alt;
    int qual;
    char* filter;
    char* info; /*comma-separated values, ID in another table*/
    char* format; /*comma-separated values, ID in another table*/
}Misc_t;

void parseArgs(int argc, char **argv) {
    extern char *optarg;
    int c;
    while ((c = getopt(argc, argv, options)) != -1) {
        switch (c) {
        case 'f':
        case 'F': datafile = optarg; break; //file to contain data
	case 'i':
	case 'I': inputfile = optarg; break; //input data
        case 'p':
        case 'P': varPath = optarg; break; //variable path string in the file
        case 'n':
        case 'N': varName = optarg; break; //variable name string in the file
	case 'r':
	case 'R': row = atoi(optarg); break; //number of genotype rows 
        case 'c':
	case 'C': col = atoi(optarg); break; //number of samples 
	default: break;
        } // switch
    } // while
} // parseArgs

void parseHeader1(string linestream, map<string,string> &headmap1,META_1 *header){
    int count=headmap1.size();
    int first=0,len=linestream.length();
    map<string,string>::iterator ret;
    header=(META_1*)realloc(header,(count+1)*sizeof(META_1));  
    
    string key=linestream.substr(2,linestream.find_first_of("=")-2);
    strcpy(header[count].field, key.c_str());
    first = linestream.find_first_of("<"); 
    if(first!=linestream.npos){
	ret = headmap1.find(key);
	if(ret==headmap1.end()){
	    headmap1.insert(pair<string,string>(key,linestream.substr(first,len-first))); 
	}else{
	    headmap1[key].append(linestream.substr(first,len-first));
	}
    }else{
	ret = headmap1.find(key);
        first=linestream.find_first_of("=")+1; /*field w/o <...>*/
	if(ret==headmap1.end()){
	    headmap1.insert(pair<string,string>(key,linestream.substr(first,len-first))); 
	}else{
	    headmap1[key].append(linestream.substr(first,len-first));
	}
    }
}

void parseHeader2(){

}

int main(int argc, char **argv){
    const rlim_t STACK_SIZE = 1000*1024*1024; 
    struct rlimit rl;
    rl.rlim_cur = STACK_SIZE;
    int ret = setrlimit(RLIMIT_STACK,&rl);

    parseArgs(argc, argv);
    if (datafile.empty() || inputfile.empty() || varName.empty()
	|| varPath.empty()  || col==0 || row==0) {
        cerr << "Usage:\n" << *argv
                  << " -f data-file-name\n"
		  << " -i input-data-file-name\n"
		  << " -n variable-name\n"
                  << " -p variable-path\n"
		  << " -r number of SNPs\n"
                  << " -c number of samples\n"
                  << endl;
        return -1;
    }

    ifstream fp(inputfile.c_str());
    if (!fp.is_open()) { //check input file
	cout << "ERROR: Failed to open the input file \"" << inputfile.c_str() << "\"" << endl;
	cout << "REPORT: Failed to complete writing the data" << endl;
	return -1;
    }

    char* line=NULL;
    map<string,string> headmap1;
    string linestream;
    hid_t file;
    META_1* header1=NULL;
    
    /*create new file*/
    file = H5Fcreate(datafile.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT); 
    for(int i=0;getline(fp,linestream);i++){
        if(linestream.substr(0,2) != "##"){ 
	    break; /*finished reading header*/
	}else if(linestream.substr(0,6) == "#CRHOM"){
	    parseHeader2();
            break;
	}else if(linestream.substr(0,8) == "##contig"){
	    break;
	}else{
	    parseHeader1(linestream,headmap1,header1); 
	}
	
        
        
    }
    for (map<string,string>::iterator x=headmap1.begin(); x!=headmap1.end(); ++x)
        std::cout << x->first << " => " << x->second << "\n\n\n";
}
