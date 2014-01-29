/*
 *Last Update: Jan. 29, 2013
 *Author: Roven Rommel B. Fuentes
 *TT-Chang Genetic Resources Center, International Rice Research Institute
 *
 *
*/

#include "parsevcf.h"
#include "queryProcessor.h"
#include "hdf5file.h"
#include "hdf5.h"
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
    status = H5Tclose(type);
    status = H5Fclose(file);
    free(data);
    return 0;
}

int getFORMAT(string datafile,string path,META_3 *&format,int &formcount,FILE *output){
    hid_t file,dset,space,type;
    herr_t status;
    hsize_t dims[1] = {0};
    int ndims;

    path += "/meta/format";
    //Open file and dataset
    file = H5Fopen(datafile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT); 
    dset = H5Dopen(file,path.c_str(),H5P_DEFAULT); 

    //get dataspace and allocate memory for read buffer
    space = H5Dget_space(dset);
    type = H5Dget_type(dset);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    format = (META_3*)malloc(dims[0]*sizeof(META_3));
    formcount = dims[0];
    status = H5Dread(dset,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,format);
 
    for(int i=0;i<dims[0];i++){
	fprintf(output,"##FORMAT=<ID=%s,",format[i].id);
       
        if(format[i].num==-1) fprintf(output,"Number=A,");
    	else if(format[i].num==-2) fprintf(output,"Number=G,");
    	else if(format[i].num==-3) fprintf(output,"Number=.,");
    	else fprintf(output,"Number=%d,",format[i].num);

        fprintf(output,"Type=%s,Description=\"%s\">\n",format[i].type,format[i].desc);
        format[i].type[0]=checkType(format[i].type); //transpose to integer for faster comparison
    }

    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(type);
    status = H5Fclose(file);
    return 0;
}

int getINFO(string datafile,string path,META_3 *&info,FILE *output){
    hid_t file,dset,space,type;
    herr_t status;
    hsize_t dims[1] = {0};
    int ndims;

    path += "/meta/info";
    //Open file and dataset
    file = H5Fopen(datafile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT); 
    dset = H5Dopen(file,path.c_str(),H5P_DEFAULT); 

    //get dataspace and allocate memory for read buffer
    space = H5Dget_space(dset);
    type = H5Dget_type(dset);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    info = (META_3*)malloc(dims[0]*sizeof(META_3));

    status = H5Dread(dset,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,info);

    for(int i=0;i<dims[0];i++){
 	fprintf(output,"##INFO=<ID=%s,",info[i].id);

        if(info[i].num==-1){ fprintf(output,"Number=A,");}
    	else if(info[i].num==-2){ fprintf(output,"Number=G,");}
    	else if(info[i].num==-3){ fprintf(output,"Number=.,");}
    	else{ fprintf(output,"Number=%d,",info[i].num);}
        fprintf(output,"Type=%s,Description=\"%s\">\n",info[i].type,info[i].desc);
    }

    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(type);
    status = H5Fclose(file);
    return 0;
}

int getContigs(string datafile,string path,vector<string> &contigs,FILE *output,int vcf){
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
	if(vcf)fprintf(output,"##contig=<ID=%s,length=%d>\n",data[i].id,data[i].len);
	contigs.push_back(data[i].id);
    }

    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(type);
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
    status = H5Tclose(type);
    status = H5Fclose(file);
    free(samples[0]);
    free(samples);
    return dims[0];
}

int writeHeader(string datafile,string path,META_3 *&format,META_3 *&info,vector<string> &contigs,int &formcount,FILE *output){
    getMeta(datafile,path,output);
    getFORMAT(datafile,path,format,formcount,output); //FORMAT
    getINFO(datafile,path,info,output); //INFO
    getContigs(datafile,path,contigs,output,1); //Contigs
    return getSampleNames(datafile,path,output); //Sample Names
}

int printFORMAT(META_3 *format,int idx,int *data,FILE *output){     
    return 0;
}

int printINFO(META_3 *info,int idx,MISC *data, FILE *output){     
    int x,pos;
    unsigned int y = data[idx].info;
    char *temp=NULL;

    temp = strtok(data[idx].infoval,";");
    for (x = 0; y; x++)
    {
        pos=log2(y&-y); //get the rightmost bit and its location
	//-Brian Kernighan's Method
  	y &= y - 1; // clear the least significant bit set 
        if(info[pos].num==0){
   	    if(y==0) fprintf(output,"%s",info[pos].id);  //just for the separator(;)
	    else fprintf(output,"%s;",info[pos].id);  
	    continue; //skip flag
        }
	if(y!=0) fprintf(output,"%s=%s;",info[pos].id,temp); //just for the separator(;)
	else fprintf(output,"%s=%s",info[pos].id,temp);   
	temp = strtok(NULL,";");
    }
    return 0;
}

int readFORMATfields(string datafile,string path,int idx,META_3 *&format,int *varloc,int row,int size,int ***&intvar){
    hid_t file,dset,space,type;
    hsize_t dims[2] = {0,0};
    herr_t status;
    int ndims;
    hsize_t offset[2]={row,0};
    hsize_t count[2]={size,0};

    path += "/FORMATfields/"; 
    path += format[idx].id; 
    //Open file and dataset
    file = H5Fopen(datafile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT); 
    dset = H5Dopen(file,path.c_str(),H5P_DEFAULT); 

    //get dataspace and allocate memory for read buffer
    space = H5Dget_space(dset);
    type = H5Dget_type(dset);
    ndims = H5Sget_simple_extent_dims(space,dims,NULL);
    count[1]=dims[1]; 
    status = H5Sselect_hyperslab(space, H5S_SELECT_SET,offset,NULL,count,NULL);
    status = H5Dread(dset,type,H5S_ALL,space,H5P_DEFAULT,intvar[varloc[idx]][0]);
    printf("%d ",intvar[idx][1][0]);
    //for(int i=0;i<size-1;i++){
	//printf("%d ",intvar[idx][i][0]); //error for last field
    //}

    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(type);
    status = H5Fclose(file);
    return 0;
}

int getRow(string datafile,string path,META_3 *format,META_3 *info,vector<string> contigs,int row,int size,int *gt,FILE *output){
    MISC *data;
    hid_t file,dset,space,type,dapl;
    herr_t status;
    hsize_t offset[1]={row};
    hsize_t count[1]={size};
    
    cout << path<< " "<<offset[0] <<" " << count[0] << " " << size;
    path += "/misc";
    //Open file and dataset
    file = H5Fopen(datafile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT); 
    dset = H5Dopen(file,path.c_str(),H5P_DEFAULT); 

    //get dataspace and allocate memory for read buffer 
    space = H5Dget_space(dset);
    type = H5Dget_type(dset);
    status = H5Sselect_hyperslab(space, H5S_SELECT_SET,offset,NULL,count,NULL);
    data = (MISC*)malloc(size*sizeof(MISC));
    
    status = H5Dread(dset,type,H5S_ALL,space,H5P_DEFAULT,data);

    for(int i=0;i<size;i++){
	fprintf(output,"\n%s\t%d\t%s\t%c\t%s\t%.2f\t%s\t",contigs[data[i].chrom].c_str(),data[i].pos,data[i].id,data[i].ref,data[i].alt,data[i].qual,data[i].filter);
	printINFO(info,i,data,output);
	printFORMAT(format,i,gt,output);
    }

    //close identifiers and free resources
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(type);
    status = H5Fclose(file);
    free(data);
    return 0;
}

int fetchData(string datafile,string fileformat,string sample,string snpbound,string varname, string varpath, vector<string> Calls,string outfile){
    FILE *output;
    ibis::horometer timer1;
    timer1.start(); 
    vector<int> col,row; //bounds
    vector<string> contigs;
    META_3 *info,*format;
    int idx,size1=0,size2=1,QUERY_CHUNK=3000; //default is 1
    int samcount=0,formcount;
    char **samnames;
    char **CHROM1=NULL; //Hapmap Chroms
    int *CHROM2=NULL; //VCF Chrom indexes
    int *POS; 
    int *varloc; //location of FORMAT fields in multidimensional variables(by type)
    char q_chrom[15]; //specific chrom

    int ***intvar=NULL; //2D to accomodate fields with same datatype
    float ***floatvar=NULL;
    char ***charvar=NULL;
    char ****stringvar=NULL;

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
	    cout<<"Sample:"<<col[0];
            if(col[0]<0){ printf("ERROR: Invalid sample bounds.\n"); return 1;}
	}else if(idx==sample.length()-1){
	    printf("ERROR: Invalid sample bounds.\n"); return 1;
	}else{
	    col.push_back(atoi(sample.substr(0,idx).c_str()));
	    col.push_back(atoi(sample.substr(idx+1,sample.length()-idx+1).c_str()));
	    if(col[0]<0 || col[1]<0 || col[0]>=col[1]){
	    	printf("ERROR: Invalid sample bounds.\n"); return 1;
	    } cout<<"Sample:"<<col[0]<<"-"<<col[1];
	    size2=col[1]-col[0]+1; //number of samples
        }
    }else{
	printf("ERROR: No sample number specified.\n"); return 1;
    }

    /*Get the SNP bounds for subset sample*/
    if(!snpbound.empty()){
	int inipos=snpbound.find(':',0); 
	if(inipos==snpbound.length()-1){ //with CHROM
	    printf("ERROR: Invalid query. SNP position not specified."); return 1;
	}
	if(inipos!=snpbound.npos) strcpy(q_chrom,snpbound.substr(0,inipos++).c_str());
	else inipos=0;

	idx=snpbound.find('-',inipos);
	if(idx==snpbound.npos || idx==snpbound.length()-1){ //single SNP
	    row.push_back(atoi(snpbound.substr(inipos,snpbound.length()-inipos).c_str()));
	    size1=1; cout<<" | row:"<<row[0]<<"\n";
	    if(row[0]<0){ printf("ERROR: Invalid SNP position.\n"); return 1;}
	}else{
	    row.push_back(atoi(snpbound.substr(inipos,idx-inipos).c_str()));
	    row.push_back(atoi(snpbound.substr(idx+1,snpbound.length()-idx+1).c_str()));
	    size1=row[1]-row[0]+1; //number of rows
	    if(row[0]<0 || row[1]<0 || row[0]>=row[1]){
	    	printf("ERROR: Invalid SNP bounds.\n"); return 1;
	    } cout<<" | row:"<< row[0] <<"-"<<row[1]<<"\n";
	} 
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
	size1=(size1==0)?dims[0]:size1; //assigned all rows if no range is specified
        block=size1/QUERY_CHUNK; 
	if(size1%QUERY_CHUNK) block++; //another block for remainders
        if(fetch){
	    outfile += ".txt";
	    if(!fileformat.compare("vcf") && getContigs(datafile,varpath,contigs,output,0)){ 
		printf("ERROR: Cannot fetch CONTIGS list.");
		return 1; //only for vcf
	    }
	}else{
	    if(!fileformat.compare("vcf")) outfile += ".vcf";
            else outfile += ".txt";
	}
        output = fopen(outfile.c_str(),"w"); //file for query result

	if(size2==1){
	    data = (int*)malloc(tempsize*sizeof(int)); //single sample
	}else{
	    data = (int*)malloc(tempsize*size2*sizeof(int)); //multiple sample
	}
        int rowtemp=0;
        for(int x=0;x<block;x++){ //block fetching to avoid memory overload
            //do the fetch and write here
            if(x+1==block && size1%QUERY_CHUNK>0){
		tempsize=size1%QUERY_CHUNK;
		free(data);
		if(size2==1) data = (int*)malloc(tempsize*sizeof(int)); 
		else{ data = (int*)malloc(tempsize*size2*sizeof(int));
		}
	    }
	    rowtemp=row[0]+tempsize; //exclusive y in x:y
            paramtemp<< dataname <<"["<<row[0]<<":"<<rowtemp<< ",";
            
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
	    	if(!fileformat.compare("vcf")){ 
		    if(x==0){
			samcount=writeHeader(datafile,varpath,format,info,contigs,formcount,output); //write Header
			//allocate variables for each type of FORMAT field
			allocFieldVar(format,formcount,tempsize,samcount,intvar,floatvar,charvar,stringvar,varloc); 
		    }
		    for(int i=0;i<formcount;i++){
        		if(format[i].num==1){
    	    		    if(format[i].type[0]=='0' || !strcmp(format[i].id,"GT")){ //integer/GT
	        		readFORMATfields(datafile,varpath,i,format,varloc,row[0],tempsize,intvar);
    	    		    }else if(format[i].type[0]=='1'){ //float
             			//readFORMATfields(varloc[i],floatvar);
    	    		    }else if(format[i].type[0]=='3'){ //character
            			//readFORMATfields(varloc[i],charvar);
    	    		    }else if(format[i].type[0]=='4'){ //string
	    			//readFORMATfields(varloc[i],stringvar);
    	    		    }
        		}else{
	    		    //readFORMATfields(varloc[i],stringvar);
		    	}
    		    }
		    getRow(datafile,varpath,format,info,contigs,row[0],tempsize,data,output);	    
	    	}else{

            	}
	    }
	    row[0]=rowtemp;
	    paramtemp.str("");
	}
   
    }
    
    timer1.stop();
    printf("REPORT: Successfully queried data.\n Total time elapsed:%f\n", timer1.realTime());

    if(!fetch){
       	free(format);
       	free(info);
	freeFieldVar(intvar,floatvar,charvar,stringvar);
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


