/*
 *Last Update: Aug. 2, 2013
 *Author: Roven Rommel B. Fuentes
 *TT-Chang Genetic Resources Center, International Rice Research Institute
 *
 *Load vcf data to HDF5.
 *This version only allows SNP calls (no indels/structural variants)
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
#define SIZE4 20
#define SIZE5 15
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
    char id[SIZE4];
    int num; /*A=-1*/
    char type[9]; /*Integer, Float, Flag, Character, String*/
    char* desc;
}META_3;

typedef struct{
    int chrom; /*index of META_3 or contig table*/
    int pos;
    char id[SIZE5];
    char ref; 
    char alt[5];
    int qual;
    char* filter;
    //unsigned int info;
    //char* infoval; /*comma-separated values, ID in another table*/
    //unsigned int format; 
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
	header[x].value = (char*)malloc(it->second.length()*sizeof(char));
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

void loadHeader2(META_2 *&contig,hid_t file,int count,map<string,int> &contigmap){
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

void parseHeaderInfoFormat(string linestream,META_3 *&temp,int count){
    int idx1=0,idx2=0;
    string num;
    temp = (META_3*)realloc(temp,count*sizeof(META_3));
    //ID
    idx1 = linestream.find("ID")+3; 
    idx2 = linestream.find_first_of(",",idx1);
    strcpy(temp[count-1].id,(linestream.substr(idx1,idx2-idx1)).c_str());
    //Number
    idx1 = linestream.find_first_of("=",idx2)+1;
    idx2 = linestream.find_first_of(",",idx1);
    num = linestream.substr(idx1,idx2-idx1);
    if(!num.compare("A")) temp[count-1].num = -1;
    else if(!num.compare("G")) temp[count-1].num = -2;
    else if(!num.compare(".")) temp[count-1].num = -3;
    else temp[count-1].num = atoi(num.c_str());
    //Type
    idx1 = linestream.find_first_of("=",idx2)+1;
    idx2 = linestream.find_first_of(",",idx1);
    strcpy(temp[count-1].type,(linestream.substr(idx1,idx2-idx1)).c_str());
    //Description
    idx1 = linestream.find_first_of("=",idx2)+2;
    idx2 = linestream.length()-2;
    temp[count-1].desc = (char*)malloc((idx2-idx1+1)*sizeof(char));
    strcpy(temp[count-1].desc,(linestream.substr(idx1,idx2-idx1)).c_str());
    //cout << temp[count-1].id << "\t" << temp[count-1].num << "\t" << temp[count-1].type << "\t" << temp[count-1].desc << "\n";
}

void loadHeaderInfoFormat(bool flag,META_3 *&data,int count,map<string,int> &map,hid_t file){
    hid_t memtype,space,dset,t1,t2,t3;
    hsize_t dim[1]={count};
    herr_t status;
    string name;
    if(flag) name= varPath + varName + "_format";
    else name= varPath + varName + "_info";
    space = H5Screate_simple(1,dim,NULL);
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,SIZE4);
    t2 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t2,9);
    t3 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t3,H5T_VARIABLE);
    memtype = H5Tcreate(H5T_COMPOUND,sizeof(META_3));
    H5Tinsert(memtype,"id",HOFFSET(META_3,id),t1);
    H5Tinsert(memtype,"number",HOFFSET(META_3,num),H5T_NATIVE_INT);
    H5Tinsert(memtype,"type",HOFFSET(META_3,type),t2);
    H5Tinsert(memtype,"description",HOFFSET(META_3,desc),t3);
    dset = H5Dcreate (file, name.c_str(), memtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(status >=0);
    status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
   
    assert(status >=0);
    /*free value pointers and variables*/
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(memtype);
    status = H5Tclose(t1);
    status = H5Tclose(t2);
    status = H5Tclose(t3);
    assert(status>=0);
    for(int x=0;x<count;x++){
        map.insert(pair<string,int>(string(data[x].id),x));
        free(data[x].desc);
    }
    free(data);
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

void parseInfoFormat(string linestream,map<string,int> ids){

}

void parseMisc(string linestream, MISC *&misc,map<string,int> infomap, map<string,int> formmap, map<string,int> contigmap,int count){
    int idx1=0,idx2=0,snpidx=count-1;
    string temp;
    misc = (MISC*)realloc(misc,count*sizeof(MISC));
    //chromID
    idx2 = linestream.find_first_of("\t",idx1+1); 
    misc[snpidx].chrom = contigmap.find(linestream.substr(idx1,idx2-idx1))->second;
    idx1=idx2+1;
    //POS
    idx2 = linestream.find_first_of("\t",idx1+1); 
    misc[snpidx].pos = atoi((linestream.substr(idx1,idx2-idx1)).c_str());
    idx1=idx2+1;
    //ID
    idx2 = linestream.find_first_of("\t",idx1+1); 
    strcpy(misc[snpidx].id,(linestream.substr(idx1,idx2-idx1)).c_str());
    idx1=idx2+1;
    //REF
    idx2 = linestream.find_first_of("\t",idx1+1); 
    temp = linestream.substr(idx1,idx2-idx1);
    if(temp.length()==1) misc[snpidx].ref = temp[0];   
    else misc[snpidx].ref = 'X';
    idx1=idx2+1;
    //ALT
    idx2 = linestream.find_first_of("\t",idx1+1); 
    strcpy(misc[snpidx].alt, (linestream.substr(idx1,idx2-idx1)).c_str());
    idx1=idx2+1;
    //QUAL
    idx2 = linestream.find_first_of("\t",idx1+1); 
    misc[snpidx].qual = atoi((linestream.substr(idx1,idx2-idx1)).c_str());
    idx1=idx2+1;
    //FILTER
    idx2 = linestream.find_first_of("\t",idx1+1); 
    temp = linestream.substr(idx1,idx2-idx1);
    misc[snpidx].filter = (char*)malloc(temp.length()*sizeof(char));
    strcpy(misc[snpidx].filter,temp.c_str());
    idx1=idx2+1;
    //INFO
    temp.clear();
    //parseInfoFormat(linestream,infomap);
    idx2 = linestream.find_first_of("\t",idx1+1); 
    cout << (linestream.substr(idx1,idx2-idx1)).c_str() << "\n";
    idx1=idx2+1;
    //FORMAT*/
    //parseInfoFormat(linestream,formmap);
    //cout <<  misc[snpidx].chrom << " " << misc[snpidx].pos << " "<< misc[snpidx].id << " " << misc[snpidx].ref << " " << misc[snpidx].alt << " " << misc[snpidx].qual << " " << misc[snpidx].filter << "\n";
}

void loadMisc(MISC *&misc,int count,hid_t file){
    hid_t memtype,space,dset,t1,t2,t3;
    hsize_t dim[1]={count};
    herr_t status;
    string name = varPath + varName + "_misc";
    space = H5Screate_simple(1,dim,NULL);
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,SIZE5);
    t2 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t2,5);
    t3 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t3,H5T_VARIABLE);
    memtype = H5Tcreate(H5T_COMPOUND,sizeof(META_3));
    H5Tinsert(memtype,"chrom",HOFFSET(MISC,chrom),H5T_NATIVE_INT);
    H5Tinsert(memtype,"pos",HOFFSET(MISC,pos),H5T_NATIVE_INT);
    H5Tinsert(memtype,"id",HOFFSET(MISC,id),t1);
    H5Tinsert(memtype,"ref",HOFFSET(MISC,ref),H5T_C_S1);
    H5Tinsert(memtype,"alt",HOFFSET(MISC,alt),t2);
    H5Tinsert(memtype,"qual",HOFFSET(MISC,qual),H5T_NATIVE_INT);
    H5Tinsert(memtype,"filter",HOFFSET(MISC,filter),t3);
    dset = H5Dcreate (file, name.c_str(), memtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(status >=0);
    status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, misc);
   
    assert(status >=0);
    /*free value pointers and variables*/
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(memtype);
    status = H5Tclose(t1);
    status = H5Tclose(t2);
    status = H5Tclose(t3);
    assert(status>=0);
    for(int x=0;x<count;x++){
	free(misc[x].filter);
    }
    free(misc);
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
    map<string,int> infomap;
    map<string,int> formmap;
    string linestream;
    hid_t file;
    META_2* contig=NULL;
    MISC* misc = NULL;
    META_3* info=NULL, *format = NULL;
    int** call1 = NULL, **call2=NULL;
    int contigcount=0,infocount=0,formcount=0;
    /*create new file*/
    file = H5Fcreate(datafile.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT); 
    for(int i=0;getline(fp,linestream);i++){
        if(linestream.substr(0,2) == "##"){
            if(linestream.substr(0,8) == "##contig"){
	    	parseHeader2(linestream,contig,++contigcount);
            }else if(linestream.substr(0,6) == "##INFO"){ 
	    	parseHeaderInfoFormat(linestream,info,++infocount);
            }else if(linestream.substr(0,8) == "##FORMAT"){
                parseHeaderInfoFormat(linestream,format,++formcount);
 	    }else{
 	    	parseHeader1(linestream,headmap); 
	    }
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
    loadHeaderInfoFormat(false,info,infocount,infomap,file);
    loadHeaderInfoFormat(true,format,formcount,formmap,file);
    int i=0, counter1=0,counter2=0;
    for(i=0,counter1=0,counter2=0;getline(fp,linestream),i<5;i++,counter2++){ 
        parseMisc(linestream,misc,infomap,formmap,contigmap,++counter1);
        parseGenotypes(linestream,call1,call2);
        //parseSubFields(linestream,);
        //break;
    }
    loadMisc(misc,counter1,file);
    contigmap.clear(); 
    infomap.clear();
    formmap.clear();
    H5Fclose(file);
    
}
