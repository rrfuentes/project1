/*
 *Last Update: Oct. 30, 2013
 *Author: Roven Rommel B. Fuentes
 *TT-Chang Genetic Resources Center, International Rice Research Institute
 *
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

static const char *options="f:F:t:gs:p:P:n:N:r:R:o:O:c:C";
static string pos;
static string sample;
static string snpbound;
static string datafile;
static string fileformat;
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
            case 't': fileformat = optarg; break;
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
    //NOTE:check dset
    
    //get datatype and its size
    filetype = H5Dget_type(dset);
    sdim = H5Tget_size(filetype);
    sdim++; /*space for null terminator*/ 
    //get dset space
    space = H5Dget_space(dset);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    if(dims[0]<num || idx1+num>dims[0]){
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
    /*Getting CHROM from VCF*/

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
    if(dims[0]<num || idx1+num>dims[0]){
	printf("ERROR:SNP range is not within the actual size.\n"); 
	return 1;
    }

    //allocate space for buffer
    CHROM = (int*)malloc(num*sizeof(int));
    
    //Read subset
    offset[0]=idx1;
    count[0]=num;
    status = H5Sselect_hyperslab(space,H5S_SELECT_SET,offset,NULL,count,NULL);
    memspace = H5Screate_simple(ndims,dims2,NULL);
    status = H5Dread(dset,fieldtype,memspace,space,H5P_DEFAULT,CHROM);

    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Sclose(memspace);
    status = H5Tclose(fieldtype);
    status = H5Fclose(file);
    return 0;
}

int getCHROM(string datafile,string path,char **&CHROM,int idx1,int num){
    hid_t file,dset,space,memspace,filetype,fieldtype,stringtype;
    herr_t status;
    hsize_t dims[1] = {0},dims2[1]={num};
    hsize_t count[1],offset[1];
    int ndims; //number of dimensions   

    path += "/misc";  
    //Open file and dataset
    file = H5Fopen(datafile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT); 
    dset = H5Dopen(file,path.c_str(),H5P_DEFAULT); 
    //create datatype for CONTIGS header
    stringtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(stringtype,20);
    fieldtype = H5Tcreate(H5T_COMPOUND,sizeof(char)*20);
    status = H5Tinsert(fieldtype,"chrom",0,stringtype);
    
    //get dset space
    space = H5Dget_space(dset);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    if(dims[0]<num || idx1+num>dims[0]){
	printf("ERROR:SNP range is not within the actual size.\n"); 
	return 1;
    }

    //allocate space for buffer
    CHROM = (char**)malloc(num*sizeof(char*));
    CHROM[0] = (char*)malloc(num*20*sizeof(char));
    for(int i=0;i<num;i++) CHROM[i] = CHROM[0] + i*20; 

    //read field
    offset[0]=idx1;
    count[0]=num;
    status = H5Sselect_hyperslab(space,H5S_SELECT_SET,offset,NULL,count,NULL);
    memspace = H5Screate_simple(ndims,dims2,NULL);
    status = H5Dread(dset,fieldtype,memspace,space,H5P_DEFAULT,CHROM[0]);
    
    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Sclose(memspace);
    status = H5Tclose(fieldtype);
    status = H5Tclose(stringtype);
    status = H5Fclose(file);
    return 0;
}

int getPOS(string datafile,string path,int *&POS,int idx1,int num){
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
    status = H5Tinsert(fieldtype,"pos",0,H5T_NATIVE_INT);
    
    //get dset space
    space = H5Dget_space(dset);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    if(dims[0]<num || idx1+num>dims[0]){
	printf("ERROR:SNP range is not within the actual size.\n"); 
	return 1;
    }

    //allocate space for buffer
    POS = (int*)malloc(num*sizeof(int));
    
    //Read subset
    offset[0]=idx1;
    count[0]=num;
    status = H5Sselect_hyperslab(space,H5S_SELECT_SET,offset,NULL,count,NULL);
    memspace = H5Screate_simple(ndims,dims2,NULL);
    status = H5Dread(dset,fieldtype,memspace,space,H5P_DEFAULT,POS);

    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Sclose(memspace);
    status = H5Tclose(fieldtype);
    status = H5Fclose(file);
    return 0;
}

int getContigsHeader(string datafile,string path,map<int,string> &CONTIGS){
    hid_t file,dset,space,filetype,fieldtype,stringtype;
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

int fetchData(string datafile,string fileformat,string sample,string snpbound,string varname, string varpath, map<int,string> Calls,string outfile){
    FILE *output;
    vector<int> col,row; //bounds
    int idx,size1=1,size2=1; //default is 1	
    char **samnames;
    char **CHROM1=NULL; 
    int *CHROM2=NULL;
    int *POS;
    map<int,string> CONTIGS;
    FQ::FileFormat model = FQ::FQ_HDF5;
    QueryProcessor* fq = new QueryProcessor(datafile, model, "", 0, "",""); 
  
    if (fq->isValid() == false) {
	printf("ERROR: Failed to initiate query processor for file.\n");
 	return 1;
    }

    if(fileformat.compare("vcf") && fileformat.compare("hapmap")){
	printf("ERROR: Invalid file format");
	return 1;
    }

    /*Get the sample bounds*/
    if(!sample.empty()){
	idx=sample.find(':',0);
	if(idx==sample.npos){
	    col.push_back(atoi(sample.substr(0,idx).c_str()));
            if(col[0]<0){
	    	printf("ERROR: Invalid sample bounds.\n");
	    	return 1;
	    }
	}else if(idx==sample.length()-1){
		printf("ERROR: Invalid sample bounds.\n");
		return 1;
	}else{
	    col.push_back(atoi(sample.substr(0,idx).c_str()));
	    col.push_back(atoi(sample.substr(idx+1,sample.length()-idx+1).c_str()));
	    if(col[0]<0 || col[1]<0 || col[0]>=col[1]){
	    	printf("ERROR: Invalid sample bounds.\n");
	    	return 1;
	    }
	    size2=col[1]-col[0]+1; //number of samples
        }
    }else{
	printf("ERROR: No sample number specified.\n");
        return 1;
    }

    /*Get the SNP bounds for subset sample*/
    if(!snpbound.empty()){
	idx=snpbound.find(':',0);
	if(idx==snpbound.npos || idx==snpbound.length()-1){
		printf("ERROR: Invalid SNP bounds.\n");
		return 1;
	}
	row.push_back(atoi(snpbound.substr(0,idx).c_str()));
	row.push_back(atoi(snpbound.substr(idx+1,snpbound.length()-idx+1).c_str()));
	if(row[0]<0 || row[1]<0 || row[0]>=row[1]){
	    printf("ERROR: Invalid SNP bounds.\n");
	    return 0;
	}
	size1=row[1]-row[0]+1; //number of samples
    }else{
	row.push_back(0);
    }

    
    int *data;
    string variable,param,temppath,dataname;
    vector<uint64_t> dims;
    FQ::DataType type;
    ostringstream paramtemp;
    
    if(!varpath.compare("/")) varpath = varpath + varname;
    else varpath = varpath + "/" + varname;

    //query GT
    if(!fileformat.compare("vcf")){
	temppath = varpath + "/FORMATfields";
  	dataname = "GT";
    }else{
 	temppath = varpath;
	dataname = "snp";
    }
   
    if (!fq->getVariableInfo(dataname.c_str(), variable, dims, &type, temppath)) {
 	printf("ERROR: Failed to get the information for variable.\n");
	return 1;
    }else {
	output = fopen(outfile.c_str(),"w");
	if(size1>1){ //set row range
	    paramtemp<< dataname <<"["<<row[0]<<":"<<row[1]+1<< ",";
	}else{
	    paramtemp<< dataname <<"[,";
            size1=dims[0];
	}

	if(size2==1){ //set sample range
	    data = (int*)malloc(size1*sizeof(int)); 
	    paramtemp << col[0] << "]";
	}else{
	    data = (int*)malloc(size1*size2*sizeof(int));
	    paramtemp << col[0] << ":" << col[1]+1 << "]"; //exclusive y in x-y
	}
	param = paramtemp.str();
	fq->getData(param,data); 
    }

    //get sample name/s, CHROMs,POSs
    if(getSampleNames(datafile,varpath,samnames,col[0],size2)) return 1;
    if(!fileformat.compare("vcf") && getContigsHeader(datafile,varpath,CONTIGS)) return 1; //only for vcf
    if(!fileformat.compare("hapmap") && getCHROM(datafile,varpath,CHROM1,row[0],size1)) return 1;
    if(!fileformat.compare("vcf") && getCHROM(datafile,varpath,CHROM2,row[0],size1)) return 1;
    if(getPOS(datafile,varpath,POS,row[0],size1)) return 1;
    
    //write to file
    if(!fileformat.compare("vcf")){ //print VCF
        fprintf(output,"CHROM\tPOS\t");
    	for(int i=0;i<size2;i++) fprintf(output,"%s\t",samnames[i]); //print sample names
    	fprintf(output,"\n");
    	for(int i=0;i<size1;i++){
	    fprintf(output,"%s\t%d\t",CONTIGS.find(CHROM2[i])->second.c_str(),POS[i]);
	    for(int j=0;j<size2;j++){
	    	fprintf(output,"%s\t",Calls[data[i*size2+j]].c_str()); //print GT
	    }
	    fprintf(output,"\n");
    	}
    }else{ //print hapmap
	fprintf(output,"CHROM\tPOS\t");
    	for(int i=0;i<size2;i++) fprintf(output,"%s\t",samnames[i]); //print sample names
    	fprintf(output,"\n");
    	for(int i=0;i<size1;i++){
	    fprintf(output,"%s\t%d\t",CHROM1[i],POS[i]);
	    for(int j=0;j<size2;j++){
	    	fprintf(output,"%s\t",Calls[data[i*size2+j]].c_str()); //print GT
	    }
	    fprintf(output,"\n");
    	}
    }
    if(CHROM1!=NULL){
	free(CHROM1[0]);
	free(CHROM1);
    }else{
    	free(CHROM2);
	CONTIGS.clear();
    }
    free(POS);
    //free(samnames[0]);
    //free(samnames);
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
    if(datafile.empty() || fileformat.empty() || varname.empty() || varpath.empty() || outfile.empty()){
		std::cerr << "Usage:\n" << *argv
                  << " -f data-file-name"
  	 	  << " -t format (vcf,hapmap)"
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
    
    if(fetch) fetchData(datafile,fileformat,sample,snpbound,varname,varpath,Calls,outfile); 	
}


