/*
 *version: 2.0
 * Load the specified samples in memory and compare(serial/pthread) against a reference column.
 * Last Updated:July 19, 2013
 * Works on data > memorysize 
 *
 *Author: Roven Rommel B. Fuentes
 *TT-Change Genetic Resources Center, International Rice Research Institute
 *
 */

#include "queryProcessor.h"
#include "hdf5file.h"
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <pthread.h>
#include <time.h>
#include <sstream>
#include <map>
#include <sys/resource.h>

#define NUMTHREADS 100
using namespace std;

static const char *options="f:F:gs:p:P:n:N:r:R:o:O:c:C";
static string pos;
static string sample;
static string snpbound;
static string datafile;
static string varpath;
static string varname;
static string outfile;
static bool fetch=false;

void returnAlleleStates(map<int,string> &Calls){
    Calls[0]="A/A";
    Calls[1]="A/T";
    Calls[2]="A/C";
    Calls[3]="A/G";
    Calls[4]="A/N";
    Calls[5]="T/A";
    Calls[6]="T/T";
    Calls[7]="T/C";
    Calls[8]="T/G";
    Calls[9]="T/N";
    Calls[10]="C/A";
    Calls[11]="C/T";
    Calls[12]="C/C";
    Calls[13]="C/G";
    Calls[14]="C/N";
    Calls[15]="G/A";
    Calls[16]="G/T";
    Calls[17]="G/C";
    Calls[18]="G/G";
    Calls[19]="G/N";
    Calls[20]="N/A";
    Calls[21]="N/T";
    Calls[22]="N/C";
    Calls[23]="N/G";
    Calls[24]="N/N";
    Calls[25]="./.";
    Calls[26]=".";
}

void parseArgs(int argc, char **argv){
    extern char *optarg;
    int i;
    while((i=getopt(argc, argv, options)) != -1){
	switch(i){
	    case 'f':
            case 'F': datafile = optarg; break;
            case 'g': fetch = true; break;
	    case 's': sample = optarg; break;
 	    case 'p':
	    case 'P': varpath = optarg; break;
	    case 'n':
	    case 'N': varname = optarg; break;
	    case 'c':
	    case 'C': pos = optarg; break; //e.g: 'x|y:z' where x=ref sample, y=sample for comparison, y:z=range
	    case 'r':
	    case 'R': snpbound = optarg; break; //e.g: x:y, where x=start, y=end
	    case 'o':
	    case 'O': outfile = optarg; break; 
	    default: break;
	}
    }
}

int fetchData(string datafile,string sample,string snpbound,string varname, string varpath, map<int,string> Calls,string outfile){
    FILE *output;
    int temp,idx1=0,idx2=0;	
    FQ::FileFormat model = FQ::FQ_HDF5;
    QueryProcessor* fq = new QueryProcessor(datafile, model, "", 0, "",""); 
  
    if (fq->isValid() == false) {
	printf("ERROR: Failed to initiate query processor for file.\n");
 	return 1;
    }
 
    //get sample index/es
    if(sample.empty()){
	printf("ERROR: No sample number specified.\n");
        return 1;
    }else{
	
	temp = sample.find_first_of("-");
	if(temp==sample.npos){
	    idx1 = atoi(sample.c_str()); cout << idx1;
	    temp = 1;
  	}else{
	    idx1 = atoi(sample.substr(0,temp).c_str());
	    temp++;
	    idx2 = atoi(sample.substr(temp,sample.length()-temp).c_str()); cout <<idx1 << " "<<idx2;
	    temp = idx2 - idx1+1;
	    if(idx2<idx1) printf("ERROR: Invalid sample range.\n");
	}
    }

    //query
    int *data;
    string variable,param;
    vector<uint64_t> dims;
    FQ::DataType type;
    ostringstream paramtemp;
    if(!varpath.compare("/")) varpath = varpath + varname + "/FORMATfields";
    else varpath = varpath + "/" + varname + "/FORMATfields";

    if (!fq->getVariableInfo("GT", variable, dims, &type, varpath)) {
 	printf("ERROR: Failed to get the information for variable.\n");
    }else {
	if(idx2==0){ //single sample
	    data = (int*)malloc(dims[0]*sizeof(int)); //dims[0] ->length of sample
	    paramtemp << "GT[," << idx1 << "]";
	}else{
	    data = (int*)malloc(dims[0]*temp*sizeof(int));
	    paramtemp << "GT[," << idx1 << ":" << idx2+1 << "]"; //exclusive y in x-y
	}
	param = paramtemp.str();
	fq->getData(param,data); 
    }
    
    //write to file
    int len=idx2-idx1;
    output = fopen(outfile.c_str(),"w");
    for(int i=0;i<dims[0];i++){
	for(int j=0;j<temp;j++){
	    fprintf(output,"%s\t",Calls[data[i*temp+j]].c_str());
	}
	fprintf(output,"\n");
    }
    free(data);
    fclose(output);
}

int main(int argc, char **argv) {
    const rlim_t STACK_SIZE = 1000*1024*1024; 
    struct rlimit rl;
    rl.rlim_cur = STACK_SIZE;
    int ret = setrlimit(RLIMIT_STACK,&rl);

    map<int, string> Calls;
    returnAlleleStates(Calls);

    parseArgs(argc, argv);
    if(datafile.empty() || varname.empty() || varpath.empty() || outfile.empty()){
		std::cerr << "Usage:\n" << *argv
                  << " -f data-file-name"
		  << " -g flag to fetch data"
	          << " -s sample/s (x ||x-y)"
		  << " -n variable-name"
                  << " -p variable-path"
		  << " -c ref&row indices ('x|y:z' where x=ref sample,x:y as row range) "
		  << " -r snp bounds (x:y, where x=start, y=end)"
		  << " -d variable-dimension (e.g. 2:2)"
		  << " -o output-file"
                  << std::endl;
    }
    if(fetch) fetchData(datafile,sample,snpbound,varname,varpath,Calls,outfile); 	
}


