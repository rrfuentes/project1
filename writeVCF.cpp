/*
 *Last Update: July 31, 2013
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
#define SIZE3 20
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
}MISC;

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

void loadHeader1(map<string,string> &headmap,hid_t file){
    /*relocate header to META_1 before loading to HDF5*/
    int x=0,count=headmap.size(); 
    META_1 header[count];
    hid_t memtype,space,dset;
    hid_t t1,t2;
    hsize_t dim[1]={count};
    herr_t status;
    string name = varPath + varName + "_meta";
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
    status = H5Tclose(t1);
    status = H5Tclose(t2);
    assert(status>=0);
    for(int x=0;x<count;x++){
	free(header[x].value);
    }
}

void loadHeader2(META_2 *contig,hid_t file,int count,map<string,int> &contigmap){
    hid_t memtype,space,dset,t1;
    hsize_t dim[1]={count};
    herr_t status;
    string name= varPath + varName + "_contigs";
    space = H5Screate_simple(1,dim,NULL);
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,SIZE2);
    memtype = H5Tcreate(H5T_COMPOUND,sizeof(META_2));
    H5Tinsert(memtype,"id",HOFFSET(META_2,id),t1);
    H5Tinsert(memtype,"len",HOFFSET(META_2,len),H5T_NATIVE_INT);
    dset = H5Dcreate (file, name.c_str(), memtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(status >=0);
    status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, contig);
   
    assert(status >=0);
    /*free value pointers and variables*/
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(memtype);
    status = H5Tclose(t1);
    assert(status>=0);
    for(int x=0;x<count;x++){
        contigmap.insert(pair<string,int>(string(contig[x].id),x));
    }
    /*for (map<string,int>::iterator it=contigmap.begin(); it!=contigmap.end();++it){
	cout << it->first << "==>" << it->second << "\n";
    }*/
    free(contig);
}

void parseHeader2(string linestream,META_2 *&contig,int count){
    int idx1 = 0, idx2 = 0;
    contig = (META_2*)realloc(contig,count*sizeof(META_2));
    idx1 = linestream.find("ID")+3;
    idx2 = linestream.find_first_of(",",idx1);
    strcpy(contig[count-1].id,(linestream.substr(idx1,idx2-idx1)).c_str());
    idx1 = linestream.find_first_of("length",idx2)+7; 
    idx2 = linestream.find_first_of(">",idx1);
    contig[count-1].len=atoi((linestream.substr(idx1,idx2-idx1)).c_str());
    //cout << contig[count-1].id << "==>" << contig[count-1].len << "\n";
}


void loadSampleNames(string linestream, hid_t file){
    char** samples=NULL;
    int idx1=0,idx2=0,x=0,count=0;
    for(x=0;x<9;x++){ //ignore fixed table headername
        idx1 = linestream.find_first_of("\t",idx1+1); 
    } 
    for(x=0;idx1<linestream.npos;x++){ 
	samples = (char**)realloc(samples,(x+1)*sizeof(char*));
        samples[x] = new char[SIZE3];
        idx2 = linestream.find_first_of("\t",idx1+1); 
        if(idx2==linestream.npos){
	    strcpy(samples[x],(linestream.substr(idx1+1,linestream.length()-idx1)).c_str());
	}else{
	    strcpy(samples[x],(linestream.substr(idx1+1,idx2-idx1-1)).c_str());
	}
        idx1=idx2; 
        //cout << samples[x].name << "\n";
    } 
    count=x;
    
    hid_t memtype,space,dset,t1;
    hsize_t dim[1]={x};
    herr_t status;
    string name= varPath + varName + "_samples";
    space = H5Screate_simple(1,dim,NULL);
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,H5T_VARIABLE);
    dset = H5Dcreate (file, name.c_str(), t1, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(status >=0);
    status = H5Dwrite(dset, t1, H5S_ALL, H5S_ALL, H5P_DEFAULT, samples);
   
    assert(status >=0);
    //free value pointers and variables
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(t1);
    assert(status>=0);
    for(x=0;x<count;x++) free(samples[x]);
    free(samples);
}

void parseInfoFormat(string linestream,INFO_FORMAT *&infoform){}

void parseMisc(string linestream, MISC *&miscdata,INFO_FORMAT *&infoform,map<string,int> contigmap){
    int idx1=0,idx2=0,count=0;
    /*chromID*/
    idx2 = linestream.find_first_of("\t",idx1+1); 
    cout << contigmap.find(linestream.substr(idx1,idx2-idx1))->second << "\t";
    idx1=idx2+1;
    /*POS*/
    idx2 = linestream.find_first_of("\t",idx1+1); 
    cout << atoi((linestream.substr(idx1,idx2-idx1)).c_str()) << "\t";
    idx1=idx2+1;
    /*ID*/
    idx2 = linestream.find_first_of("\t",idx1+1); 
    cout << (linestream.substr(idx1,idx2-idx1)).c_str() << "\t";
    idx1=idx2+1;
    /*REF*/
    idx2 = linestream.find_first_of("\t",idx1+1); 
    cout << (linestream.substr(idx1,idx2-idx1)).c_str() << "\t";
    idx1=idx2+1;
    /*ALT*/
    idx2 = linestream.find_first_of("\t",idx1+1); 
    cout << (linestream.substr(idx1,idx2-idx1)).c_str() << "\t";
    idx1=idx2+1;
    /*QUAL*/
    idx2 = linestream.find_first_of("\t",idx1+1); 
    cout << atoi((linestream.substr(idx1,idx2-idx1)).c_str()) << "\t";
    idx1=idx2+1;
    /*FILTER*/
    idx2 = linestream.find_first_of("\t",idx1+1); 
    cout << (linestream.substr(idx1,idx2-idx1)).c_str() << "\t";
    idx1=idx2+1;
    /*INFO*/
    parseInfoFormat(linestream,infoform);
    idx2 = linestream.find_first_of("\t",idx1+1); 
    cout << (linestream.substr(idx1,idx2-idx1)).c_str() << "\n";
    idx1=idx2+1;
    /*FORMAT*/
    
}

void parseGenotypes(string linestream,int **call1,int **call2){}

int main(int argc, char **argv){
    const rlim_t STACK_SIZE = 1000*1024*1024; 
    struct rlimit rl;
    //rl.rlim_cur = STACK_SIZE;
    //int ret = setrlimit(RLIMIT_STACK,&rl);

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
    map<string,string> headmap;
    map<string,int> contigmap;
    string linestream;
    hid_t file;
    META_2* contig=NULL;
    MISC* miscdata = NULL;
    INFO_FORMAT* infoform = NULL;
    int** call1 = NULL, **call2=NULL;
    int contigcount=0;
    /*create new file*/
    file = H5Fcreate(datafile.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT); 
    for(int i=0;getline(fp,linestream);i++){
        if(linestream.substr(0,8) == "##contig"){
	    parseHeader2(linestream,contig,++contigcount);
        }else if(linestream.substr(0,2) == "##"){ 
	    parseHeader1(linestream,headmap); 
	}else if(linestream.substr(0,6) == "#CHROM"){ 
            loadSampleNames(linestream,file);
	    break; 
	}else{
            cout << "ERROR: Failed to write the input file. Unknown format of header entry. \"" << "\"" << endl;
	    cout << "REPORT: Failed to complete writing the data" << endl;
	    return 1; /*finished reading header*/
	}
    }
    loadHeader1(headmap,file); 
    loadHeader2(contig,file,contigcount,contigmap);
    for(int i=0,counter1=0,counter2=0;getline(fp,linestream),i<5;i++,counter1++,counter2++){ 
        parseMisc(linestream,miscdata,infoform,contigmap);
        parseGenotypes(linestream,call1,call2);
        //parseSubFields(linestream,);
        //break;
    }
    contigmap.clear(); 
    H5Fclose(file);
    
}
