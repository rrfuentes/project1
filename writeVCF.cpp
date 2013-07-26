/*
 *Last Update: July 25, 2013
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

#define SIZE1 20
#define SIZE2 20
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
    char field[SIZE1];
    char* value;  
}META_1; 

/*contigs*/
typedef struct{ 
    char id[SIZE2];
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

void parseHeader1(string linestream, map<string,string> &headmap){
    int count=headmap.size();
    int first=0,len=linestream.length();
    map<string,string>::iterator ret;
    
    string key=linestream.substr(2,linestream.find_first_of("=")-2);
    first = linestream.find_first_of("<"); 
    if(first!=linestream.npos){
	ret = headmap.find(key);
	if(ret==headmap.end()){
	    headmap.insert(pair<string,string>(key,linestream.substr(first,len-first))); 
	}else{
	    headmap[key].append(linestream.substr(first,len-first));
	}
    }else{
	ret = headmap.find(key);
        first=linestream.find_first_of("=")+1; /*field w/o <...>*/
	if(ret==headmap.end()){
	    headmap.insert(pair<string,string>(key,linestream.substr(first,len-first))); 
	}else{
	    headmap[key].append(linestream.substr(first,len-first));
	}
    }
}

void loadHeader1(map<string,string> &headmap,META_1 *header,hid_t file){
    /*relocate header to META_1 before loading to HDF5*/
    int x=0,count=headmap.size();
    hid_t memtype,space,dset;
    hid_t t1,t2;
    hsize_t dim[1]={count};
    herr_t status;
    string name = varName + "metainfo";
    for (map<string,string>::iterator it=headmap.begin(); it!=headmap.end();++it){
	strcpy(header[x].field, it->first.c_str());
	header[x].value = (char*)malloc(it->second.length()*sizeof(char)+1);
        strcpy(header[x].value, it->second.c_str()); 
        //cout << header[x].field << " ==> " << header[x].value << "\n\n\n";
        headmap.erase(it);
        x++;
    }
    /*Load to HDF5*/
    space = H5Screate_simple(1,dim,NULL);
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,SIZE1);
    t2 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t2,H5T_VARIABLE);
    memtype = H5Tcreate(H5T_COMPOUND,sizeof(META_1)); //create compound datatype
    H5Tinsert(memtype,"field",HOFFSET(META_1,field),t1);
    H5Tinsert(memtype,"value",HOFFSET(META_1,value),t2);
    dset = H5Dcreate (file, name.c_str(), memtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, header);
    assert(status >=0);
    /*free value pointers and variables*/
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(memtype);
    assert(status>=0);
    for(int x=0;x<count;x++){
	free(header[x].value);
    }
}

void parseHeader2(string linestream,META_2 *header,int count){
    int idx1 = 0, idx2 = 0;
    header = (META_2*)realloc(header,count*sizeof(META_2));
    idx1 = linestream.find("ID")+3;
    idx2 = linestream.find_first_of(",",idx1);
    strcpy(header[count-1].id,(linestream.substr(idx1,idx2-idx1)).c_str());
    idx1 = linestream.find_first_of("length",idx2)+7; 
    idx2 = linestream.find_first_of(">",idx1);
    header[count-1].len=atoi((linestream.substr(idx1,idx2-idx1)).c_str());
    cout << header[count-1].id << "==>" << header[count-1].len << "\n";
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
    META_2* header2=NULL;
    int contigcount=0;
    /*create new file*/
    file = H5Fcreate(datafile.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT); 
    for(int i=0;getline(fp,linestream);i++){
        if(linestream.substr(0,2) != "##"){ 
	    continue; /*finished reading header*/
	}else if(linestream.substr(0,8) == "##contig"){
	    //parseHeader2(linestream,header2,++contigcount);
        }else if(linestream.substr(0,6) == "#CHROM"){ printf("final");
	    break; 
	}else{
	    parseHeader1(linestream,headmap1); 
	}
    }
    int count=headmap1.size();
    header1=(META_1*)malloc(count*sizeof(META_1));
    loadHeader1(headmap1,header1,file); 
    free(header1);
    free(header2);
    H5Fclose(file);
    
}
