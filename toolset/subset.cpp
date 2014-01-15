/*
 *Last Update: Jan. 15, 2013
 *Author: Roven Rommel B. Fuentes
 *TT-Chang Genetic Resources Center, International Rice Research Institute
 *
 *
*/

#include "parsevcf.h"
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

using namespace std;

static const char *options="f:F:t:gs:p:P:n:N:r:R:o:O:";
static string sample;
static string snpbound;
static string datafile;
static string fileformat;
static string varpath;
static string varname;
static string outfile;
static bool fetch=false;

void returnAlleleStates(vector<string> &Calls){
    Calls.push_back("A/A");
    Calls.push_back("A/T");
    Calls.push_back("A/C");
    Calls.push_back("A/G");
    Calls.push_back("A/N");
    Calls.push_back("T/A");
    Calls.push_back("T/T");
    Calls.push_back("T/C");
    Calls.push_back("T/G");
    Calls.push_back("T/N");
    Calls.push_back("C/A");
    Calls.push_back("C/T");
    Calls.push_back("C/C");
    Calls.push_back("C/G");
    Calls.push_back("C/N");
    Calls.push_back("G/A");
    Calls.push_back("G/T");
    Calls.push_back("G/C");
    Calls.push_back("G/G");
    Calls.push_back("G/N");
    Calls.push_back("N/A");
    Calls.push_back("N/T");
    Calls.push_back("N/C");
    Calls.push_back("N/G");
    Calls.push_back("N/N");
    Calls.push_back("./.");
    Calls.push_back(".");
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
	    case 's': sample = optarg; break; //sample/s (x ||x-y) 
 	    case 'p':
	    case 'P': varpath = optarg; break;
	    case 'n':
	    case 'N': varname = optarg; break;
	    case 'r':
	    case 'R': snpbound = optarg; break; //(x-y, where x=start, y=end) or (CHROM:x-y)
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

int getContigsHeader(string datafile,string path,vector<string> &contigs){
    hid_t file,dset,space,fieldtype,stringtype;
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
	contigs.push_back(id[i]);
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

int getMeta(string datafile,string path,FILE *output){
    META_1 *data;
    hid_t file,dset,space,type;
    herr_t status;
    hsize_t dims[1] = {0};
    int ndims;
    char *temp=NULL;

    path += "/meta/meta";
    //Open file and dataset
    file = H5Fopen(datafile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT); 
    dset = H5Dopen(file,path.c_str(),H5P_DEFAULT); 

    //get dataspace and allocate memory for read buffer
    space = H5Dget_space(dset);
    type = H5Dget_type(dset);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    data = (META_1*)malloc(dims[0]*sizeof(META_1));
    
    status = H5Dread(dset,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

    for(int i=0;i<dims[0];i++){
	temp = strtok(data[i].value,";");
        if(temp==NULL){
	    fprintf(output,"##%s=%s\n",data[i].field,data[i].value);
        }else{  //multiple entries (e.g. FILTER)
	    while(temp!=NULL){
		fprintf(output,"##%s=%s\n",data[i].field,temp);
  		temp = strtok(NULL,";");
            }
	}
	
    }

    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Fclose(file);
    free(data);
    return 0;
}

int getINFOFORMAT(string datafile,string path,vector<string> &vec,FILE *output,bool flag){
    META_3 *data;
    hid_t file,dset,space,type;
    herr_t status;
    hsize_t dims[1] = {0};
    int ndims;

    if(flag) path += "/meta/format";
    else path += "/meta/info";
    //Open file and dataset
    file = H5Fopen(datafile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT); 
    dset = H5Dopen(file,path.c_str(),H5P_DEFAULT); 

    //get dataspace and allocate memory for read buffer
    space = H5Dget_space(dset);
    type = H5Dget_type(dset);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    data = (META_3*)malloc(dims[0]*sizeof(META_3));
    
    status = H5Dread(dset,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

    for(int i=0;i<dims[0];i++){
	if(flag){
	    fprintf(output,"##FORMAT=<ID=%s,",data[i].id);
	    vec.push_back(data[i].id);
        }else{
 	    fprintf(output,"##INFO=<ID=%s,",data[i].id);
            vec.push_back(data[i].id);
	}

        if(data[i].num==-1) fprintf(output,"Number=A,");
    	else if(data[i].num==-2) fprintf(output,"Number=G,");
    	else if(data[i].num==-3) fprintf(output,"Number=.,");
    	else fprintf(output,"Number=%d,",data[i].num);
        fprintf(output,"Type=%s,Description=\"%s\">\n",data[i].type,data[i].desc);
    }

    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Fclose(file);
    free(data);
    return 0;
}

int getContigs(string datafile,string path,vector<string> &vec,FILE *output){
    META_2 *data;
    hid_t file,dset,space,type;
    herr_t status;
    hsize_t dims[1] = {0};
    int ndims;

    path += "/meta/contigs";
    //Open file and dataset
    file = H5Fopen(datafile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT); 
    dset = H5Dopen(file,path.c_str(),H5P_DEFAULT); 

    //get dataspace and allocate memory for read buffer
    space = H5Dget_space(dset);
    type = H5Dget_type(dset);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    data = (META_2*)malloc(dims[0]*sizeof(META_2));
    
    status = H5Dread(dset,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

    for(int i=0;i<dims[0];i++){
	fprintf(output,"##contig=<ID=%s,length=%d>\n",data[i].id,data[i].len);
    }

    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Fclose(file);
    free(data);
    return 0;
}

int getSampleNames(string datafile,string path,FILE *output){
    hid_t file,dset,space,type;
    herr_t status;
    hsize_t dims[1] = {0};
    int ndims;
    char** samples=NULL;

    path += "/samples";
    //Open file and dataset
    file = H5Fopen(datafile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT); 
    dset = H5Dopen(file,path.c_str(),H5P_DEFAULT); 

    //get dataspace and allocate memory for read buffer
    space = H5Dget_space(dset);
    type = H5Dget_type(dset);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    
    samples = (char**)malloc(dims[0]*sizeof(char*));
    samples[0] = (char*)malloc(dims[0]*SIZE3*sizeof(char));
    for(int i=0;i<dims[0];i++){ samples[i]=samples[0]+i*SIZE3;}
	
    status = H5Dread(dset,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,samples[0]);
   
    fprintf(output,"CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for(int i=0;i<dims[0];i++){
	fprintf(output,"\t%s",samples[i]);
    }

    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Fclose(file);
    free(samples[0]);
    free(samples);
    return 0;
}

int writeHeader(string datafile,string path,vector<string> &format,vector<string> &info, vector<string> &contigs,FILE *output){
    getMeta(datafile,path,output);
    getINFOFORMAT(datafile,path,format,output,1); //FORMAT
    getINFOFORMAT(datafile,path,info,output,0); //INFO
    getContigs(datafile,path,contigs,output); //Contigs
    getSampleNames(datafile,path,output); //Sample Names
}

int getMISC(string datafile,string path,vector<string> format,vector<string> info,vector<string> contigs,int row,int size,FILE *output){
    MISC *data;
    hid_t file,dset,memspace,space,type;
    herr_t status;
    hsize_t dims[1] = {0};
    int ndims;
    hssize_t offset[1]={row};
    hsize_t count[1]={size};
    cout << "hey";
    path += "/meta/meta";
    //Open file and dataset
    file = H5Fopen(datafile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT); 
    dset = H5Dopen(file,path.c_str(),H5P_DEFAULT); 

    //get dataspace and allocate memory for read buffer
    memspace = H5Screate_simple (1, count, NULL); 
    space = H5Dget_space(dset);
    type = H5Dget_type(dset);
    //ndims = H5Sget_simple_extent_dims(space,dims,NULL);=
    status = H5Sselect_hyperslab(space, H5S_SELECT_SET,
 			(const hsize_t*)offset,NULL, count, NULL);
    data = (MISC*)malloc(size*sizeof(MISC));
    
    status = H5Dread(dset,type,memspace,space,H5P_DEFAULT,data);

    for(int i=0;i<dims[0];i++){
	    fprintf(output,"##%s\n",contigs[data[i].chrom].c_str());
    }

    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Fclose(file);
    free(data);
    return 0;
}

int fetchData(string datafile,string fileformat,string sample,string snpbound,string varname, string varpath, vector<string> Calls,string outfile){
    FILE *output;
    ibis::horometer timer1;
    timer1.start(); 
    vector<int> col,row; //bounds
    vector<string> format,info,contigs;
    int idx,size1=1,size2=1,QUERY_CHUNK=3000; //default is 1	
    char **samnames;
    char **CHROM1=NULL; //Hapmap Chroms
    int *CHROM2=NULL; //VCF Chrom indexes
    int *POS;
    char q_chrom[15]; //specific chrom
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
	idx=sample.find('-',0);
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
		idx=snpbound.find('-');
		if(idx==snpbound.npos || idx==snpbound.length()-1){
		    printf("ERROR: Invalid SNP bounds.\n");
		    return 1;
		}
	}else{
		strcpy(q_chrom,snpbound.substr(0,idx).c_str());
		idx=snpbound.find('-',idx);
	}
	cout << q_chrom;
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
	 //TO DO: check the 1st and last storage index of bounds within chromX
        // CHROM=q_chrom && pos<row[1]
        int block=0,tempsize=QUERY_CHUNK;
	size1=(size1==1)?dims[0]:size1; //assigned all rows if no range is specified
        block=size1/QUERY_CHUNK; 
	if(size1%QUERY_CHUNK) block++; //another block for remainders
        if(fetch){
	    outfile += ".txt";
	    output = fopen(outfile.c_str(),"w");
	    if(!fileformat.compare("vcf") && getContigsHeader(datafile,varpath,contigs)){ 
		printf("ERROR: Cannot fetch CONTIGS list.");
		return 1; //only for vcf
	    }
	}else{
	    if(!fileformat.compare("vcf")){ outfile += ".vcf";
            }else{ outfile += ".txt";}
	}

	if(size2==1){
	    data = (int*)malloc(tempsize*sizeof(int)); //single sample
	}else{
	    data = (int*)malloc(tempsize*size2*sizeof(int)); //multiple sample
	}

        for(int x=0;x<block;x++){ //block fetching to avoid memory overload
            //do the fetch and write here
            if(x+1==block && size1%QUERY_CHUNK>0){
		tempsize=size1%QUERY_CHUNK;
		free(data);
		if(size2==1) data = (int*)malloc(tempsize*sizeof(int)); 
		else{ data = (int*)malloc(tempsize*size2*sizeof(int));
		}
	    }
	    row[1]=row[0]+tempsize; //exclusive y in x:y
            paramtemp<< dataname <<"["<<row[0]<<":"<<row[1]<< ",";
            
	    if(size2==1) paramtemp << col[0] << "]";
	    else paramtemp << col[0] << ":" << col[1]+1 << "]"; //exclusive y in x:y 

	    param = paramtemp.str();
	    fq->getData(param,data); 

	    if(fetch){ //CHROM\tPOS\tGT
		//get sample name/s, CHROMs,POSs
               	if(getSampleNames(datafile,varpath,samnames,col[0],size2)) return 1;
    	    	if(getPOS(datafile,varpath,POS,row[0],tempsize)) return 1;
            	if(!fileformat.compare("vcf")){ //print VCF
		    if(getCHROM(datafile,varpath,CHROM2,row[0],tempsize)){
			printf("ERROR: Cannot fetch CHROM colunm.");
			return 1;
		    }
		    if(x==0){ //print in the 1st block iteration
            	    	fprintf(output,"CHROM\tPOS\t");
    	    	     	for(int i=0;i<size2;i++) 
			    fprintf(output,"%s\t",samnames[i]); //print sample names
		    }
    	    	    
    	   	    for(int i=0;i<tempsize;i++){
	    	     	fprintf(output,"\n%s\t%d\t",contigs[CHROM2[i]].c_str(),POS[i]);
	    	    	for(int j=0;j<size2;j++){
	    	    	    fprintf(output,"%s\t",Calls[data[i*size2+j]].c_str()); //print GT
	     	    	}
    	    	    }
		    free(CHROM2);
    	    	}else{ //print hapmap
		    if(getCHROM(datafile,varpath,CHROM1,row[0],tempsize)){
			printf("ERROR: Cannot fetch CHROM colunm.");
			return 1;
		    }
		    if(x==0){ //print in the 1st block iteration
	    	    	fprintf(output,"CHROM\tPOS\t");
    	    	    	for(int i=0;i<size2;i++) 
			    fprintf(output,"%s\t",samnames[i]); //print sample names
		    }
    	    	    
    	    	    for(int i=0;i<tempsize;i++){
	    	    	fprintf(output,"\n%s\t%d\t",CHROM1[i],POS[i]);
	    	    	for(int j=0;j<size2;j++)
	    	    	    fprintf(output,"%s\t",Calls[data[i*size2+j]].c_str()); //print GT
    	    	    }
		    free(CHROM1[0]);
		    free(CHROM1);
    	       	}
		free(POS);
		free(samnames[0]);
       		free(samnames);
	    }else{  //MISC\tGT
   	    	output = fopen(outfile.c_str(),"w");
	    	if(!fileformat.compare("vcf")){ 
		    if(x==0) writeHeader(datafile,varpath,format,info,contigs,output); //write Header
		    getMISC(datafile,varpath,format,info,contigs,row[0],tempsize,output);
	    	}else{

            	}
	    }
	    row[0]=row[1];
	    paramtemp.str("");
	}
   
    }

    timer1.stop();
    printf("REPORT: Successfully queried data.\n Total time elapsed:%f\n", timer1.realTime());

    if(!fetch){
       format.clear();
       info.clear();
    }
    if(!fileformat.compare("vcf")) contigs.clear();
    free(data);
    fclose(output);
    delete(fq);
}

int main(int argc, char **argv) {
    const rlim_t STACK_SIZE = 1000*1024*1024; 
    struct rlimit rl;
    rl.rlim_cur = STACK_SIZE;
    int ret = setrlimit(RLIMIT_STACK,&rl);

    vector<string> Calls;
    returnAlleleStates(Calls);

    parseArgs(argc, argv);
    if(datafile.empty() || fileformat.empty() || varname.empty() || varpath.empty() || outfile.empty()){
		std::cerr << "Usage:\n" << *argv
                  << " -f data-file-name"
  	 	  << " -t format (vcf,hapmap)"
		  << " -g flag to fetch data -- CHROM\tPOS\tGT"
	          << " -s sample/s (x ||x-y)"
		  << " -n variable-name"
                  << " -p variable-path"
		  << " -r snp bounds (x-y, where x=start, y=end) or (CHROM:x-y)" 
		  << " -o output-file"
                  << std::endl;
    }
    
    fetchData(datafile,fileformat,sample,snpbound,varname,varpath,Calls,outfile); 	
}


