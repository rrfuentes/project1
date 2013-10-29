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
#include "parsevcf.h"
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

int getSampleNames(string datafile,string path,char **&names,int idx1,int num){
    hid_t file,dset,space,memspace,filetype,memtype;
    herr_t status;
    hsize_t dims[1] = {0},dims2[1]={num};
    hsize_t count[1],offset[1];
    size_t sdim;
    int ndims; //number of dimensions   

    path += "/samples";
    //Open file and dataset
    file = H5Fopen(datafile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT); 
    dset = H5Dopen(file,path.c_str(),H5P_DEFAULT); 
    //get datatype and its size
    filetype = H5Dget_type(dset);
    sdim = H5Tget_size(filetype);
    sdim++; /*space for null terminator*/ 
    //get dset space
    space = H5Dget_space(dset);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    if(dims[0]<(idx1+num)){
	printf("ERROR:Sample range is not within the actual size.\n"); 
	return 1;
    }

    //create memory datatype
    memtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(memtype,sdim);

    //allocate space for buffer
    names = (char**)malloc(num*sizeof(char*));
    names[0] = (char*)malloc(num*sdim*sizeof(char));
    for(int i=0;i<num;i++) names[i] = names[0] + i*sdim; 

    //Read subset
    offset[0]=idx1;
    count[0]=num;
    status = H5Sselect_hyperslab(space,H5S_SELECT_SET,offset,NULL,count,NULL);
    memspace = H5Screate_simple(ndims,dims2,NULL);
    status = H5Dread(dset,memtype,memspace,space,H5P_DEFAULT,names[0]);

    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Sclose(memspace);
    status = H5Tclose(filetype);
    status = H5Tclose(memtype);
    status = H5Fclose(file);
    return 0;
}

int getCHROM(string datafile,string path,int *&CHROM,int idx1,int num){
    hid_t file,dset,space,memspace,fieldtype;
    herr_t status;
    hsize_t dims[1] = {0},dims2[1]={num};
    hsize_t count[1],offset[1];
    int ndims; //number of dimensions   

    path += "/misc";
    //Open file and dataset
    file = H5Fopen(datafile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT); 
    dset = H5Dopen(file,path.c_str(),H5P_DEFAULT); 
    //create datatype for CHROM field
    fieldtype = H5Tcreate(H5T_COMPOUND,sizeof(int));
    status = H5Tinsert(fieldtype,"chrom",0,H5T_NATIVE_INT);
    
    //get dset space
    space = H5Dget_space(dset);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    if(dims[0]<(idx1+num)){
	printf("ERROR:Sample range is not within the actual size.\n"); 
	return 1;
    }

    //allocate space for buffer
    CHROM = (int*)malloc(dims[0]*sizeof(int));
    
    //Read subset
    offset[0]=idx1;
    count[0]=num;
    status = H5Sselect_hyperslab(space,H5S_SELECT_SET,offset,NULL,count,NULL);
    memspace = H5Screate_simple(ndims,dims2,NULL);
    status = H5Dread(dset,fieldtype,memspace,space,H5P_DEFAULT,CHROM);

    //for(int i=0;i<num;i++) cout <<CHROM[i]<<"\n";
    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Sclose(memspace);
    status = H5Tclose(fieldtype);
    status = H5Fclose(file);
    return 0;
}

int getPOS(string datafile,string path,int *&POS,int idx1,int num){
    return 0;
}

int getContigsHeader(string datafile,string path,map<int,string> &CONTIGS){
    hid_t file,dset,space,memspace,filetype,fieldtype,stringtype;
    herr_t status;
    hsize_t dims[1] = {0};
    int ndims; //number of dimensions 
    char **id;  

    path += "/meta/contigs";
    //Open file and dataset
    file = H5Fopen(datafile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT); 
    dset = H5Dopen(file,path.c_str(),H5P_DEFAULT); 
    //create datatype for CONTIGS header
    stringtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(stringtype,20);
    fieldtype = H5Tcreate(H5T_COMPOUND,sizeof(char)*20);
    status = H5Tinsert(fieldtype,"id",0,stringtype);
    
    //get dset space
    space = H5Dget_space(dset);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);

    //allocate space for buffer
    id = (char**)malloc(dims[0]*sizeof(char*));
    id[0] = (char*)malloc(dims[0]*20*sizeof(char));
    for(int i=0;i<dims[0];i++) id[i] = id[0] + i*20; 
    //read field
    status = H5Dread(dset,fieldtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,id[0]);
    //save IDs to map
    for(int i=0;i<dims[0];i++){
	CONTIGS.insert(pair<int,string>(i,string(id[i])));
    }
    //close identifiers and free resources
    free(id[0]);
    free(id);
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(fieldtype);
    status = H5Tclose(stringtype);
    status = H5Fclose(file);
    return 0;
}

int fetchData(string datafile,string sample,string snpbound,string varname, string varpath, map<int,string> Calls,string outfile){
    FILE *output;
    int temp,idx1=0,idx2=0;	
    char **samnames; 
    int *CHROM;
    int *POS;
    map<int,string> CONTIGS;
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
	    idx2 = atoi(sample.substr(temp,sample.length()-temp).c_str()); 
	    temp = idx2-idx1+1; 
	    if(idx2<=idx1){
		printf("ERROR: Invalid sample range.\n");
		return 1;
	    }
	}
    }

    if(!varpath.compare("/")) varpath = varpath + varname;
    else varpath = varpath + "/" + varname;

    //get sample name/s, CHROMs,POSs
    if(getSampleNames(datafile,varpath,samnames,idx1,temp)) return 1;
    if(getContigsHeader(datafile,varpath,CONTIGS)) return 1;
    if(getCHROM(datafile,varpath,CHROM,idx1,5)) return 1;
    if(getPOS(datafile,varpath,POS,idx1,temp)) return 1;

    //query
    int *data;
    string variable,param;
    vector<uint64_t> dims;
    FQ::DataType type;
    ostringstream paramtemp;
    
    varpath+= "/FORMATfields";
    if (!fq->getVariableInfo("GT", variable, dims, &type, varpath)) {
 	printf("ERROR: Failed to get the information for variable.\n");
	return 1;
    }else {
	output = fopen(outfile.c_str(),"w");
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
    fprintf(output,"CHROM\tPOS\t");
    for(int i=0;i<temp;i++) fprintf(output,"%s\t",samnames[i]); //print sample names
    fprintf(output,"\n");
    for(int i=0;i<dims[0];i++){
	fprintf(output,"%s\t",CONTIGS.find(CHROM[i])->second.c_str());
	for(int j=0;j<temp;j++){
	    fprintf(output,"%s\t",Calls[data[i*temp+j]].c_str()); //print GT
	}
	fprintf(output,"\n");
    }
    CONTIGS.clear();
    free(CHROM);
    free(samnames[0]);
    free(samnames);
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


