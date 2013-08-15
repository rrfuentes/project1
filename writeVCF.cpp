/*
 *Last Update: Aug. 6, 2013
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

#define CHUNKSIZE1 5000
#define CHUNKSIZE2_1 10000 
#define CHUNKSIZE2_2 2 
#define SNP_CHUNK_CACHE 268435456 /*250MB*/
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
static int indelcount = 0;

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
    char alt[6];
    int qual;
    char* filter;
    unsigned int info;
    char* infoval; /*semi-colon-separated values, ID in another table*/
    unsigned int format; 
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
	default: break;
        } // switch
    } // while
} // parseArgs

int checkType(char* type){
    if(!strcmp(type,"Integer")) return 0;
    else if(!strcmp(type,"Float")) return 1;
    else if(!strcmp(type,"Flag")) return 2;
    else if(!strcmp(type,"Character")) return 3;
    else if(!strcmp(type,"String")) return 4;
}

int checkNumber(META_3 data){
    //cout << data.type << data.num << "\n";
}

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
    string name = varPath + varName + "/meta/meta";
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
    string name= varPath + varName + "/meta/contigs";
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

void loadHeaderInfoFormat(bool flag,META_3 *&data,int count,map<string,int> &map,vector<pair<int,int> > vec, hid_t file){
    hid_t memtype,space,dset,t1,t2,t3;
    hsize_t dim[1]={count};
    herr_t status;
    string name;
    if(flag) name= varPath + varName + "/meta/format";
    else name= varPath + varName + "/meta/info";
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
        vec.push_back(pair<int,int>(checkType(data[x].type),0));
        //cout << vec[x].first << " " << vec[x].second << "\n";
        free(data[x].desc);
    }
    cout << "\n";
    free(data);
}

int loadSampleNames(string linestream, hid_t file){
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
    hsize_t dim[1]={count};
    herr_t status;
    string name= varPath + varName + "/samples";
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
    for(int x=0;x<count;x++) free(samples[x]);
    free(samples);
    return count;
    
}

void parseInfo(string &linestream,int idx1, int idx2, MISC*& misc, int pos, map<string,int> ids){
    int temp1=0,temp2=idx2;
    unsigned int flag=0;
    string values(""); 
    //cout << idx2 << " "<< linestream.find_first_of("=\t",idx1+1) << " "<< linestream.substr(idx1,idx2-idx1) << "\n";
    while(idx1<idx2){
        temp1 = linestream.find_first_of("=\t;",idx1+1); //with \t or ; to include "." or flags(e.g. STR)
        flag |= 1<< (ids.find(linestream.substr(idx1,temp1-idx1))->second);
	if(temp1!=idx2 && linestream[temp1]!=';'){
            temp2 = linestream.find_first_of(";\t",temp1+1);
            values += linestream.substr(temp1+1,temp2-temp1-1)+ ";"; 
	}else{
	    temp2=idx2;
	}
        idx1=temp2+1;
    }
    
    misc[pos].info = flag;
    misc[pos].infoval = (char*)malloc((values.length()+1)*sizeof(char));
    strcpy(misc[pos].infoval,values.c_str());
}

unsigned int parseFormat(string linestream,int idx1, int idx2, map<string,int> ids){
    int temp =0, flags = 0;
    while(idx1<idx2){
        temp = linestream.find_first_of(":\t",idx1+1); 
        flags |= 1<< (ids.find(linestream.substr(idx1,temp-idx1))->second);
        idx1 = temp+1;
    }
    /*for(int x=0;x<ids.size();x++){
	cout << (flags & 1<<x) << " ";
    }
    cout << "\n";*/
    return flags;
}

int parseAlt(string linestream,int idx1, int idx2,MISC *&misc,int pos,vector<char> &callvec){
    int temp1=idx1,temp2=idx2; 
    if(idx2-idx1>1){ //check ALT
	if(idx2-idx1>5){
	    strcpy(misc[pos].alt,"X"); 
	    return 1;
	}else{
	    while(temp1<idx2){
            	temp2 = linestream.find_first_of(",\t",temp1+1); /*ALT must be single-char nucleotides*/
            	if(temp2-temp1>1){ //indels/structural variants are not supported
                    //cout <<linestream.substr(idx1,idx2-idx1) << "\n";
		    callvec.push_back('X');
		    strcpy(misc[pos].alt,"X");
	            return 1;
	    	}
		callvec.push_back(linestream[temp1]);
            	temp1 = temp2+1;
            }
	}
    }else{
	callvec.push_back(linestream[idx1]);
    }
    strcpy(misc[pos].alt,(linestream.substr(idx1,idx2-idx1)).c_str()); 
    return 0;
}

int parseMisc(string linestream, MISC *&misc,map<string,int> infomap, vector<pair<int,int> > &infovec, map<string,int> formmap, vector<pair<int,int> > &formvec, map<string,int> contigmap,vector<char> &callvec, int snpidx){
    int idx1=0,idx2=0;
    bool indel=0;
    unsigned int flag=0;
    string temp; 

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
    if(temp.length()==1){
 	misc[snpidx].ref = temp[0];  
        callvec.push_back(temp[0]); 
    }else{ 
	misc[snpidx].ref = 'X';
	callvec.push_back('X'); 
	indel=1;
    }
    idx1=idx2+1;
    //ALT
    idx2 = linestream.find_first_of("\t",idx1+1); 
    indel = (indel || parseAlt(linestream,idx1,idx2,misc,snpidx,callvec));
    idx1=idx2+1;
    //QUAL
    idx2 = linestream.find_first_of("\t",idx1+1); 
    misc[snpidx].qual = atoi((linestream.substr(idx1,idx2-idx1)).c_str());
    idx1=idx2+1;
    //FILTER
    idx2 = linestream.find_first_of("\t",idx1+1); 
    temp = linestream.substr(idx1,idx2-idx1);
    misc[snpidx].filter = (char*)malloc((temp.length()+1)*sizeof(char));
    strcpy(misc[snpidx].filter,temp.c_str()); 
    idx1=idx2+1;
    //INFO
    idx2 = linestream.find_first_of("\t",idx1+1); 
    parseInfo(linestream,idx1,idx2,misc,snpidx,infomap);
    idx1=idx2+1;
    //FORMAT
    idx2 = linestream.find_first_of("\t",idx1+1); 
    misc[snpidx].format = parseFormat(linestream,idx1,idx2,formmap);
    //cout <<  misc[snpidx].chrom << " " << misc[snpidx].pos << " "<< misc[snpidx].id << " " << misc[snpidx].ref << " " << misc[snpidx].alt << " " << misc[snpidx].qual << " " << misc[snpidx].filter << "\n";
    if(indel) indelcount++;
    
    return idx2;
}

void clearMisc(MISC *&misc,int count){
    for(int x=0;x<count;x++){
	free(misc[x].filter);
        free(misc[x].infoval);
    }
}


void setH5MiscFormat(hid_t file,hid_t &memtype,hid_t &space,hid_t &dset,hid_t &cparms,hid_t &dataprop,int chunk,string &name){
    hid_t t1,t2,t3;
    herr_t status;
    hsize_t dim[1]={chunk};
    hsize_t maxdim[1] = {H5S_UNLIMITED};
    hsize_t chkdim[1] = {chunk};
    hssize_t offset[1] = {0};

    name = varPath + varName + "/misc";
    space = H5Screate_simple(1,dim,maxdim); //create data space
    cparms = H5Pcreate(H5P_DATASET_CREATE); //create chunk 
    status = H5Pset_chunk(cparms,1,chkdim);

    //Compound MISC datatype
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,SIZE5);
    t2 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t2,5);
    t3 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t3,H5T_VARIABLE);
    memtype = H5Tcreate(H5T_COMPOUND,sizeof(MISC));
    H5Tinsert(memtype,"chrom",HOFFSET(MISC,chrom),H5T_NATIVE_INT);
    H5Tinsert(memtype,"pos",HOFFSET(MISC,pos),H5T_NATIVE_INT);
    H5Tinsert(memtype,"id",HOFFSET(MISC,id),t1);
    H5Tinsert(memtype,"ref",HOFFSET(MISC,ref),H5T_C_S1);
    H5Tinsert(memtype,"alt",HOFFSET(MISC,alt),t2);
    H5Tinsert(memtype,"qual",HOFFSET(MISC,qual),H5T_NATIVE_INT);
    H5Tinsert(memtype,"filter",HOFFSET(MISC,filter),t3);
    H5Tinsert(memtype,"infobit",HOFFSET(MISC,info),H5T_NATIVE_UINT);
    H5Tinsert(memtype,"infoval",HOFFSET(MISC,infoval),t3);
    H5Tinsert(memtype,"formatbit",HOFFSET(MISC,format),H5T_NATIVE_UINT);

    //create dataset access property list for MISC dataset
    dataprop = H5Pcreate(H5P_DATASET_ACCESS);
    status = H5Pset_chunk_cache(dataprop,H5D_CHUNK_CACHE_NSLOTS_DEFAULT,
             SNP_CHUNK_CACHE,H5D_CHUNK_CACHE_W0_DEFAULT); /*set snp chunk size*/

    dset = H5Dcreate (file, name.c_str(), memtype, space, H5P_DEFAULT, cparms, dataprop);
    assert(status >=0);
   
   /*free value pointers and variables*/
    status = H5Tclose(t1);
    status = H5Tclose(t2);
    status = H5Tclose(t3);
    assert(status>=0);
}

void setH5CallFormat(hid_t file, hid_t &space,hid_t &dset,hid_t &cparms,hid_t &dataprop,int size1,int size2,string &name){
    herr_t status;
    hsize_t dim[2]={size1,size2};
    hsize_t maxdim[2] = {H5S_UNLIMITED,size2};
    hsize_t chkdim[2] = {size1,CHUNKSIZE2_2};
    hssize_t offset[2] = {0,0};

    name = varPath + varName + "/call";
    space = H5Screate_simple(2,dim,maxdim); //create data space
    cparms = H5Pcreate(H5P_DATASET_CREATE); //create chunk 
    status = H5Pset_szip(cparms,H5_SZIP_NN_OPTION_MASK,8); //compression
    status = H5Pset_chunk(cparms,2,chkdim);
    
    //create dataset access property list for MISC dataset
    dataprop = H5Pcreate(H5P_DATASET_ACCESS);
    status = H5Pset_chunk_cache(dataprop,H5D_CHUNK_CACHE_NSLOTS_DEFAULT,
             SNP_CHUNK_CACHE*2,H5D_CHUNK_CACHE_W0_DEFAULT); //set snp chunk size

    dset = H5Dcreate (file, name.c_str(), H5T_NATIVE_INT, space, H5P_DEFAULT, cparms, dataprop);
    assert(status >=0);
}

void parseGenotypes(string linestream,int **&call,int snpidx,int offset,vector<char> &callvec,map<string,int> states){
    int idx1 = offset, idx2 = offset,x=0,counter=0; 
    int last = linestream.length(); 
    string temp;
    if(callvec[0]=='X' || callvec[1]=='X' || callvec[0]=='.' || callvec[1]=='.'){ //filter indels and structural variants
 	while(idx1<last){
            idx2 = linestream.find_first_of("\t",idx1+2); 
	    if(idx2==linestream.npos) idx2=last;
	    call[snpidx][counter] = -1;
	    idx1 = idx2+1;
            counter++;
        }
    }else{ //SNPs
	while(idx1<last){
            idx2 = linestream.find_first_of("\t",idx1+2); 
	    if(idx2==linestream.npos) idx2=last;
            //cout << linestream.substr(idx1,idx2-idx1) << "\t";

	    //0-Ref 1..N-Alt
            //get Ref
	    x=linestream[idx1]-48; //convert from char to int  
	    temp = callvec[x]; //get the call
            //get Alt 
	    x=linestream[idx1+2]-48; //convert from char to int  
	    temp += callvec[x]; //get the call

            call[snpidx][counter] = states.find(temp)->second; //save the call code
            idx1 = idx2+1;
            counter++;
    	}
    }
    
    callvec.clear();
}

void setAlleleStates(map<string,int> &Calls){
    Calls["AA"]=0;
    Calls["AT"]=1;
    Calls["AC"]=2;
    Calls["AG"]=3;
    Calls["AN"]=4;
    Calls["TA"]=5;
    Calls["TT"]=6;
    Calls["TC"]=7;
    Calls["TG"]=8;
    Calls["TN"]=9;
    Calls["CA"]=10;
    Calls["CT"]=11;
    Calls["CC"]=12;
    Calls["CG"]=13;
    Calls["CN"]=14;
    Calls["GA"]=15;
    Calls["GT"]=16;
    Calls["GC"]=17;
    Calls["GG"]=18;
    Calls["GN"]=19;
    Calls["NA"]=20;
    Calls["NT"]=21;
    Calls["NC"]=22;
    Calls["NG"]=23;
    Calls["NN"]=24;
}


int main(int argc, char **argv){
    const rlim_t STACK_SIZE = 1000*1024*1024; 
    struct rlimit rl;
    rl.rlim_cur = STACK_SIZE;
    int ret = setrlimit(RLIMIT_STACK,&rl);

    parseArgs(argc, argv);
    if (datafile.empty() || inputfile.empty() || varName.empty()
	|| varPath.empty()) {
        cerr << "Usage:\n" << *argv
                  << " -f data-file-name\n"
		  << " -i input-data-file-name\n"
		  << " -n variable-name\n"
                  << " -p variable-path\n"
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
    map<string,int> states; //allele states
    map<string,int> contigmap;
    map<string,int> infomap;
    map<string,int> formmap;
    vector<pair<int,int> > infovec;
    vector<pair<int,int> > formvec;
    vector<char> callvec;
    string linestream,groupname;
    hid_t file,group,group_sub;
    META_2* contig=NULL;
    MISC* misc = NULL;
    META_3* info=NULL, *format = NULL;
    int **call=NULL;
    int contigcount=0,infocount=0,formcount=0,samcount=0;
    
    /*create new file*/
    groupname=varPath+varName;
    file = H5Fcreate(datafile.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT); 
    group = H5Gcreate (file, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    groupname += "/meta"; 
    group_sub = H5Gcreate (file, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
            samcount = loadSampleNames(linestream,file);
	    break; 
	}else{
            cout << "ERROR: Failed to write the input file. Unknown format of header entry. \"" << "\"" << endl;
	    cout << "REPORT: Failed to complete writing the data" << endl;
	    return 1; /*finished reading header*/
	}
    }
    loadHeader1(headmap,file); 
    loadHeader2(contig,file,contigcount,contigmap);
    loadHeaderInfoFormat(false,info,infocount,infomap,infovec,file);
    loadHeaderInfoFormat(true,format,formcount,formmap,formvec,file);
    
    int i=0, counter1=0,counter2=0,lastidx;
    hid_t memtype1,space1,memspace1,dset1,cparms1,dataprop1;
    hid_t space2,memspace2,dset2,cparms2,dataprop2;
    herr_t status;
    hsize_t dim1[1]={CHUNKSIZE1}; 
    hsize_t dim2[2]={CHUNKSIZE2_1,samcount};
    hsize_t maxdim1[1] = {H5S_UNLIMITED};
    hsize_t maxdim2[2] = {H5S_UNLIMITED,samcount};
    hssize_t offset1[1]={0};
    hssize_t offset2[2]={0,0};
    hsize_t count1[1]={CHUNKSIZE1};
    hsize_t count2[2]={CHUNKSIZE2_1,samcount};
    hsize_t newsize1[1]={CHUNKSIZE1};
    hsize_t newsize2[2]={CHUNKSIZE2_1,samcount};
    string name1;
    string name2;

    //create var for misc and genotype calls
    call=(int**)malloc(CHUNKSIZE2_1*sizeof(int*));
    call[0]=(int*)malloc(CHUNKSIZE2_1*samcount*sizeof(int));
    for(int x=0;x<CHUNKSIZE2_1;x++) call[x]=call[0]+x*samcount;

    misc = (MISC*)malloc(CHUNKSIZE1*sizeof(MISC));
    if(misc==NULL){
 	cout << "ERROR: Insufficient memory. Adjust the chunksize." << endl;
        cout << "REPORT: Failed to complete writing the data" << endl;
	return 1;
    }
    setAlleleStates(states); 	

    for(i=0;fp!=NULL && i<10000;i++){ 
        getline(fp,linestream);
	if(fp!=NULL){
            lastidx = parseMisc(linestream,misc,infomap,infovec,formmap,formvec,contigmap,callvec,counter1);
            counter1++;
            parseGenotypes(linestream,call,counter2,lastidx+1,callvec,states);
            counter2++;
        }
        if(counter1==CHUNKSIZE1 || fp==NULL){
            if(i+1==CHUNKSIZE1){ //first slab
                //set compound datatype and spaces
  		setH5MiscFormat(file,memtype1,space1,dset1,cparms1,dataprop1,counter1,name1);
                memspace1 = H5Dget_space(dset1);
                //write first slab
            	status = H5Dwrite(dset1, memtype1, H5S_ALL, H5S_ALL, H5P_DEFAULT, misc); 
		clearMisc(misc,counter1); 
	    }else{ 
		offset1[0]=newsize1[0];
          	//extend data later for the 2nd to the last slab
                if(fp==NULL && counter1>0){
 		     count1[0] = counter1;
		     newsize1[0]=newsize1[0]+counter1;
		     dim1[0]=counter1;
    		     memspace1 = H5Screate_simple(1,dim1,maxdim1); 
		}else{
		     newsize1[0]=newsize1[0]+counter1; 
		}
		status = H5Dset_extent(dset1,newsize1);
	 	space1 = H5Dget_space(dset1);
                status = H5Sselect_hyperslab(space1, H5S_SELECT_SET,
 			(const hsize_t*)offset1,NULL, count1, NULL);
		status = H5Dwrite(dset1,memtype1,memspace1,space1,H5P_DEFAULT,misc);
                clearMisc(misc,counter1);  
                status = H5Sclose(space1); 
	    } 
	    cout << "MISC-Chunk" << (i+1)/CHUNKSIZE1 << "\n";
            counter1=0;
        }

        if(counter2==CHUNKSIZE2_1 || fp==NULL){
            if(i+1==CHUNKSIZE2_1){ //first slab
                //set compound datatype and spaces
  		setH5CallFormat(file,space2,dset2,cparms2,dataprop2,counter2,samcount,name2);
                memspace2 = H5Dget_space(dset2);
                //write first slab
            	status = H5Dwrite(dset2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &call[0][0]); 
	    }else{ 
		offset2[0]=newsize2[0];
          	//extend data later for the 2nd to the last slab
                if(fp==NULL && counter2>0){
 		     count2[0] = counter2;
		     newsize2[0]=newsize2[0]+counter2;
		     dim2[0]=counter2;
    		     memspace2 = H5Screate_simple(2,dim2,maxdim2); cout << "here";
		}else{
		     newsize2[0]=newsize2[0]+counter2; 
		}
		status = H5Dset_extent(dset2,newsize2);
	 	space2 = H5Dget_space(dset2);
                status = H5Sselect_hyperslab(space2, H5S_SELECT_SET,
 			(const hsize_t*)offset2,NULL, count2, NULL);
		status = H5Dwrite(dset2,H5T_NATIVE_INT,memspace2,space2,H5P_DEFAULT,&call[0][0]); 
                status = H5Sclose(space2); 
	    } 
	    cout << "Call-Chunk" << (i+1)/CHUNKSIZE2_1 << "\n";
            counter2=0;
        }
        //parseSubFields(linestream,)
       
    }

    cout << "Indel/Structural Variant Count:" << indelcount << "\n";
    /*free value pointers and variables*/
    contigmap.clear(); 
    infomap.clear();
    formmap.clear();
    infovec.clear();
    formvec.clear();
    status = H5Dclose(dset1);
    status = H5Tclose(memtype1);
    status = H5Sclose(memspace1);
    status = H5Pclose(cparms1);
    status= H5Pclose (dataprop1);
    status = H5Dclose(dset2);
    status = H5Sclose(memspace2);
    status = H5Pclose(cparms2);
    status= H5Pclose (dataprop2);
    free(call);
    free(call[0]);	
    free(misc);
    H5Fclose(file);
    
}
