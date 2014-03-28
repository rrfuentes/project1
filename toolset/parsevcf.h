/*
 *Last Update: March. 28, 2014
 *Author: Roven Rommel B. Fuentes
 *TT-Chang Genetic Resources Center, International Rice Research Institute
 *
 *Space efficient if FORMAT fields are uniform in all SNPs
 *Maximum of 32 INFO/FORMAT fields
*/

//make a preparser to compute the maximum length for string dataset

#include "hdf5.h"
#include "indexBuilder.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <queue>
#include <assert.h>
#include <sys/resource.h>
#include <vector>
#include <math.h>

#define CHUNKSIZE1 10000
#define CHUNKSIZE2_1 50000
#define CHUNKSIZE2_2 1//3//200//1
#define SNP_CHUNK_CACHE 268435456 /*250MB*/
#define SIZE1 20
#define SIZE2 20
#define SIZE3 20//50//20
#define SIZE4 20
#define SIZE5 2//25//2
#define SIZE6 2//10//2
#define SIZE7 150//70//100//80
#define SIZE8 25

using namespace std;

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
    int num; /*A=-1,G=-2,'.'=-3*/
    char type[9]; /*Integer, Float, Flag, Character, String*/
    char* desc;
}META_3;

typedef struct{
    char id[SIZE5];
    char ref; 
    char alt[6];
    float qual;
    char filter[SIZE6];
    unsigned int info;
    char infoval[SIZE7]; /*semi-colon-separated values, ID in another table*/
    unsigned int format; 
}MISC;

char checkType(char* type){
    if(!strcmp(type,"Integer")) return '0';
    else if(!strcmp(type,"Float")) return '1';
    else if(!strcmp(type,"Flag")) return '2';
    else if(!strcmp(type,"Character")) return '3';
    else if(!strcmp(type,"String")) return '4';
}

int checkNumber(META_3 data){
    //cout << data.type << data.num << "\n";
    /*if(!num.compare("A")) temp[count-1].num = -1;
    else if(!num.compare("G")) temp[count-1].num = -2;
    else if(!num.compare(".")) temp[count-1].num = -3;
    else temp[count-1].num = atoi(num.c_str());*/
}

int parseHeader1(string linestream, map<string,int> &headmap,vector<string> &meta1_val){
    int count=headmap.size();
    int first=0,len=linestream.length();
    map<string,int>::iterator ret;
    
    
    first = linestream.find_first_of("=");
    string key=linestream.substr(2,first-2); 
    if(first!=linestream.npos){ 
        first++;
	ret = headmap.find(key);
	if(ret==headmap.end()){
	    headmap.insert(pair<string,int>(key,headmap.size()));
	    meta1_val.push_back(linestream.substr(first,len-first)); 
	}else{ //multiple entries (e.g. FILTER)
	    meta1_val[headmap.find(key)->second].append(";"); 
	    meta1_val[headmap.find(key)->second].append(linestream.substr(first,len-first));
	}
    }else{ 
	cout << "ERROR: Failed to write one header row." << endl;
	exit(EXIT_FAILURE);
    }
}

void loadHeader1(hid_t file,map<string,int> &headmap,vector<string> &meta1_val,string gpath){
    /*relocate header to META_1 before loading to HDF5*/
    int x=0,count=headmap.size(); 
    META_1 header[count];
    hid_t memtype,space,dset;
    hid_t t1,t2;
    hsize_t dim[1]={count};
    herr_t status;
    string name = gpath + "/meta"; 
    map<string,int>::const_iterator it;
    for(it=headmap.begin(); it!=headmap.end();++it){
	strcpy(header[it->second].field, it->first.c_str());
	header[it->second].value = (char*)malloc((meta1_val[it->second].length()+1)*sizeof(char));
        strcpy(header[it->second].value, meta1_val[it->second].c_str());
        //cout << header[x].field << " ==> " << header[x].value << "\n\n\n";
        //x++;
    }
    /*Load to HDF5*/
    space = H5Screate_simple(1,dim,NULL);
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,SIZE1);
    t2 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t2,H5T_VARIABLE); //variable length (VLEN)
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
    headmap.clear();
    meta1_val.clear();
}

void loadHeader2(hid_t file,META_2 *&contig,int count,map<string,int> &contigmap,string gpath){
    hid_t memtype,space,dset,t1;
    int entry = 0;
    if(count==0) entry = contigmap.size();
    else entry = count; 
    hsize_t dim[1]={entry}; 
    herr_t status; 
    string name= gpath + "/contigs";
    space = H5Screate_simple(1,dim,NULL);
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,SIZE2);
    memtype = H5Tcreate(H5T_COMPOUND,sizeof(META_2));
    H5Tinsert(memtype,"id",HOFFSET(META_2,id),t1);
    H5Tinsert(memtype,"len",HOFFSET(META_2,len),H5T_NATIVE_INT);
    dset = H5Dcreate (file, name.c_str(), memtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, contig);
    assert(status >=0);
    /*free value pointers and variables*/
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(memtype);
    status = H5Tclose(t1);
    assert(status>=0);
    if(count!=0){ //Use only if contigs are declared in header
    	for(int x=0;x<count;x++){
            contigmap.insert(pair<string,int>(string(contig[x].id),x));
    	}
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
    temp[count-1].desc = (char*)malloc((idx2-idx1+2)*sizeof(char));
    strcpy(temp[count-1].desc,(linestream.substr(idx1,idx2-idx1)).c_str());
    //cout << temp[count-1].id << "\t" << temp[count-1].num << "\t" << temp[count-1].type << "\t" << temp[count-1].desc << "\n";
}

void loadHeaderInfoFormat(hid_t file,bool flag,META_3 *&field,int count,map<string,int> &map,string gpath){
    hid_t memtype,space,dset,t1,t2,t3;
    hsize_t dim[1]={count}; 
    herr_t status;
    string name; 
    if(flag) name= gpath + "/format";
    else name= gpath + "/info";
    space = H5Screate_simple(1,dim,NULL);
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,SIZE4);
    t2 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t2,9); /*Integer, Float, Flag, Character, String*/
    t3 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t3,H5T_VARIABLE); //variable length (VLEN)
    memtype = H5Tcreate(H5T_COMPOUND,sizeof(META_3));
    H5Tinsert(memtype,"id",HOFFSET(META_3,id),t1);
    H5Tinsert(memtype,"number",HOFFSET(META_3,num),H5T_NATIVE_INT);
    H5Tinsert(memtype,"type",HOFFSET(META_3,type),t2);
    H5Tinsert(memtype,"description",HOFFSET(META_3,desc),t3); 
    dset = H5Dcreate (file, name.c_str(), memtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, field);
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
        map.insert(pair<string,int>(string(field[x].id),x));
        field[x].type[0]=checkType(field[x].type);
        //field[x].type[1]='\0';
        free(field[x].desc);
    }  
    cout << "\n";
}

int loadSampleNames(hid_t file,string linestream, string inipath){
    char** samples=NULL;
    int idx1=0,idx2=0,x=0,count=0,samplesize=0;
    vector<string> temp;
    for(x=0;x<9;x++){ //ignore fixed table headername
        idx1 = linestream.find_first_of("\t",idx1+1); 
    } 
    for(x=0;idx1<linestream.npos;x++){ //use vector first since number of sample is unknown
        idx2 = linestream.find_first_of("\t",idx1+1); 
        if(idx2==linestream.npos){
	    temp.push_back(linestream.substr(idx1+1,linestream.length()-idx1));
	}else{
	    temp.push_back(linestream.substr(idx1+1,idx2-idx1-1));
	}
        idx1=idx2; 
        //cout << samples[x].name << "\n";
    } 
    samplesize = temp.size();
    samples = (char**)malloc(samplesize*sizeof(char*));
    samples[0] = (char*)malloc(samplesize*SIZE3*sizeof(char));
    
    ostringstream tempname;
    for(x=0;x<samplesize;x++){ 
        //tempname << "Sample" << x; //rename    
	samples[x]=samples[0]+x*SIZE3;
	if(SIZE3<temp[x].length()){
	    cout << "ERROR: Insufficient array space for samplenames.\n" << endl;
	    exit(EXIT_FAILURE);
	}
        strcpy(samples[x],temp[x].c_str());
	//strcpy(samples[x],tempname.str().c_str());
        //tempname.str("");
    } 
    count=x;
    
    hid_t memtype,space,dset,t1;
    hsize_t dim[1]={count};
    herr_t status;
    string name= inipath + "/samples";
    space = H5Screate_simple(1,dim,NULL);
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,SIZE3);
    dset = H5Dcreate (file, name.c_str(), t1, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, t1, H5S_ALL, H5S_ALL, H5P_DEFAULT, samples[0]);
    assert(status >=0);
    //free value pointers and variables
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(t1);
    assert(status>=0);
    free(samples[0]);
    free(samples);
    return count;
    
}

int parseInfo(string &linestream,int idx1, int idx2, MISC*& misc, int pos, map<string,int> ids){
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
    if(SIZE7<=values.length()){ 
	cout << "\nERROR:Storage size for SNP InfoVal is insufficient.\n";
	exit(EXIT_FAILURE);
    }
    strcpy(misc[pos].infoval,values.c_str());
    return 0;
}

unsigned int parseFormat(string linestream,int idx1, int idx2, map<string,int> ids,queue<int> &tokenidx){
    int temp =0, flags = 0,tokenpos=0;
    idx1+=3; //skip GT for flagging
    
    while(idx1<idx2){
        temp = linestream.find_first_of(":\t",idx1+1); 
	tokenidx.push(tokenpos++); //needed to track transfer of FORMAT field values in parseFORMATfield(by type); cannot use iterator in a loop since GT is placed in unknown order in FORMAT header.
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
            callvec.push_back('X');
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

int parseMisc(string linestream, MISC *&misc, int **&CHR_POS,map<string,int> infomap, map<string,int> formmap,META_2 *&contig,int contigcount,map<string,int> &contigmap,vector<char> &callvec, int snpidx,int &indelcount,queue<int> &tokenidx){
    int idx1=0,idx2=0;
    bool indel=0;
    unsigned int flag=0;
    string temp; 

    //chromID
    idx2 = linestream.find_first_of("\t",idx1+1); 
    temp = linestream.substr(idx1,idx2-idx1);
    if(contigcount==0){ //no contigs in header/contig will still be used to index CHROM values
	if(contigmap.empty() || contigmap.end()==contigmap.find(temp)){
	    contigmap.insert(pair<string,int>(temp,contigmap.size()));
            contig = (META_2*)realloc(contig,contigmap.size()*sizeof(META_2));
	    strcpy(contig[contigmap.size()-1].id,temp.c_str());
	    contig[contigmap.size()-1].len=0;
	    CHR_POS[snpidx][0] = contigmap.size()-1;
	}else{ 
	    CHR_POS[snpidx][0] = contigmap.find(temp)->second;
	}
    }else{
	CHR_POS[snpidx][0] = contigmap.find(temp)->second;
    }
    idx1=idx2+1;
    //POS
    idx2 = linestream.find_first_of("\t",idx1+1); 
    CHR_POS[snpidx][1] = atoi((linestream.substr(idx1,idx2-idx1)).c_str());
    idx1=idx2+1;
    //ID
    idx2 = linestream.find_first_of("\t",idx1+1); 
    if(SIZE5<idx2-idx1){ 
	cout << "\nERROR:Storage size for SNP ID is insufficient.";
	exit(EXIT_FAILURE);
    }
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
    indel = (parseAlt(linestream,idx1,idx2,misc,snpidx,callvec) || indel);
    idx1=idx2+1;
    //QUAL
    idx2 = linestream.find_first_of("\t",idx1+1); 
    misc[snpidx].qual = atof((linestream.substr(idx1,idx2-idx1)).c_str()); 
    idx1=idx2+1;
    //FILTER
    idx2 = linestream.find_first_of("\t",idx1+1); 
    temp = linestream.substr(idx1,idx2-idx1);
    strcpy(misc[snpidx].filter,temp.c_str()); 
    idx1=idx2+1;
    //INFO
    idx2 = linestream.find_first_of("\t",idx1+1); 
    int ret = parseInfo(linestream,idx1,idx2,misc,snpidx,infomap);
    if(ret) exit(EXIT_FAILURE);
    idx1=idx2+1;
    //FORMAT
    idx2 = linestream.find_first_of("\t",idx1+1); 
    misc[snpidx].format = parseFormat(linestream,idx1,idx2,formmap,tokenidx);
    //cout <<  misc[snpidx].chrom << " " << misc[snpidx].pos << " "<< misc[snpidx].id << " " << misc[snpidx].ref << " " << misc[snpidx].alt << " " << misc[snpidx].qual << " " << misc[snpidx].filter << "\n";
    if(indel) indelcount++;
    
    return idx2; //return last linestream index
}


void setH5CHRPOSFormat(hid_t file,hid_t &space,hid_t &dset,hid_t &cparms,hid_t &dataprop,int chunk,string inipath){
    herr_t status;
    hsize_t dim[2]={chunk,2};
    hsize_t maxdim[2] = {H5S_UNLIMITED,2};
    hsize_t chkdim[2] = {chunk,1};

    string name = inipath + "/pos";
    space = H5Screate_simple(2,dim,maxdim); //create data space
    cparms = H5Pcreate(H5P_DATASET_CREATE); //create chunk 
    status = H5Pset_chunk(cparms,2,chkdim);
    status = H5Pset_deflate(cparms,8); //compression

    //create dataset access property list for MISC dataset
    dataprop = H5Pcreate(H5P_DATASET_ACCESS);
    status = H5Pset_chunk_cache(dataprop,9973,
             SNP_CHUNK_CACHE*2,H5D_CHUNK_CACHE_W0_DEFAULT); /*set snp chunk size*/
    assert(status >=0);
    dset = H5Dcreate (file, name.c_str(), H5T_NATIVE_INT, space, H5P_DEFAULT, cparms, dataprop);
    //printf("%i",dset);
}

void setH5MiscFormat(hid_t file,hid_t &memtype,hid_t &space,hid_t &dset,hid_t &cparms,hid_t &dataprop,int chunk,string inipath){
    hid_t t1,t2,t3,t4;
    herr_t status;
    hsize_t dim[1]={chunk};
    hsize_t maxdim[1] = {H5S_UNLIMITED};
    hsize_t chkdim[1] = {chunk};

    string name = inipath + "/misc";
    space = H5Screate_simple(1,dim,maxdim); //create data space
    cparms = H5Pcreate(H5P_DATASET_CREATE); //create chunk 
    status = H5Pset_chunk(cparms,1,chkdim);
    status = H5Pset_deflate(cparms,8); //compression
     
    //Compound MISC datatype
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,SIZE5);
    t2 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t2,5);
    t3 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t3,SIZE6);
    t4 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t4,SIZE7);
    memtype = H5Tcreate(H5T_COMPOUND,sizeof(MISC));
    H5Tinsert(memtype,"id",HOFFSET(MISC,id),t1);
    H5Tinsert(memtype,"ref",HOFFSET(MISC,ref),H5T_C_S1);
    H5Tinsert(memtype,"alt",HOFFSET(MISC,alt),t2);
    H5Tinsert(memtype,"qual",HOFFSET(MISC,qual),H5T_NATIVE_FLOAT);
    H5Tinsert(memtype,"filter",HOFFSET(MISC,filter),t3);
    H5Tinsert(memtype,"infobit",HOFFSET(MISC,info),H5T_NATIVE_UINT); //unsigned int
    H5Tinsert(memtype,"infoval",HOFFSET(MISC,infoval),t4);
    H5Tinsert(memtype,"formatbit",HOFFSET(MISC,format),H5T_NATIVE_UINT); //unsigned int

    //create dataset access property list for MISC dataset
    dataprop = H5Pcreate(H5P_DATASET_ACCESS);
    status = H5Pset_chunk_cache(dataprop,9973,
             SNP_CHUNK_CACHE*2,H5D_CHUNK_CACHE_W0_DEFAULT); /*set snp chunk size*/
    assert(status >=0);
    dset = H5Dcreate (file, name.c_str(), memtype, space, H5P_DEFAULT, cparms, dataprop);
    //printf("%i",dset);
   
   /*free value pointers and variables*/
    status = H5Tclose(t1);
    status = H5Tclose(t2);
    status = H5Tclose(t3);
    status = H5Tclose(t4);
    assert(status>=0);
}

void setH5FORMATfield(hid_t file, hid_t *&fieldtype_a,hid_t *&space_a,hid_t *&memspace_a,hid_t *&dset_a,hid_t &cparms,hid_t &dataprop,int CHUNK1,int CHUNK2,int samplesize,string inipath,META_3 *format,int formcount){
    herr_t status;
    hsize_t dim[2]={CHUNK1,samplesize};
    hsize_t maxdim[2] = {H5S_UNLIMITED,samplesize};
    hsize_t chkdim[2] = {CHUNK1,CHUNK2};

    //max of 32 fields
    fieldtype_a = (hid_t*)malloc(32*sizeof(hid_t)); 
    space_a = (hid_t*)malloc(32*sizeof(hid_t)); 
    memspace_a = (hid_t*)malloc(32*sizeof(hid_t)); 
    dset_a = (hid_t*)malloc(32*sizeof(hid_t)); 

    //CHUNKING and CACHING
    cparms = H5Pcreate(H5P_DATASET_CREATE); //create chunk 
    status = H5Pset_chunk(cparms,2,chkdim);
    status = H5Pset_deflate(cparms,8); //compression
    
    //create dataset access property list for MISC dataset
    dataprop = H5Pcreate(H5P_DATASET_ACCESS);
    status = H5Pset_chunk_cache(dataprop,9973,
             SNP_CHUNK_CACHE*2,H5D_CHUNK_CACHE_W0_DEFAULT); //set snp chunk size
    assert(status >=0); 

    for(int i=0;i<formcount;i++){ 
	string name = inipath + "/FORMATfields/" + format[i].id; 	
        //No flag type for FORMAT fields
        space_a[i] = H5Screate_simple(2,dim,maxdim); //create data space
        if(format[i].num==1){
            //cparms and dataprop for GT is used for other fields 
    	    if(format[i].type[0]=='0' || !strcmp(format[i].id,"GT")){ //integer/GT
	        fieldtype_a[i]=H5Tcopy(H5T_NATIVE_INT);
            	dset_a[i] = H5Dcreate (file, name.c_str(), H5T_NATIVE_INT, space_a[i], H5P_DEFAULT, cparms, dataprop);
    	    }else if(format[i].type[0]=='1'){ //float
             	fieldtype_a[i]=H5Tcopy(H5T_NATIVE_FLOAT);
            	dset_a[i] = H5Dcreate (file, name.c_str(), H5T_NATIVE_FLOAT, space_a[i], H5P_DEFAULT, cparms, dataprop);
    	    }else if(format[i].type[0]=='3'){ //character
            	fieldtype_a[i]=H5Tcopy(H5T_NATIVE_CHAR);
            	dset_a[i] = H5Dcreate (file, name.c_str(), H5T_NATIVE_CHAR, space_a[i], H5P_DEFAULT, cparms, dataprop);
    	    }else if(format[i].type[0]=='4'){ //string
	    	fieldtype_a[i]=H5Tcopy(H5T_C_S1);
                H5Tset_size(fieldtype_a[i],SIZE8);
            	dset_a[i] = H5Dcreate (file, name.c_str(), fieldtype_a[i], space_a[i], H5P_DEFAULT, cparms, dataprop);
    	    }
        }else{
	    fieldtype_a[i]=H5Tcopy(H5T_C_S1);
            H5Tset_size(fieldtype_a[i],SIZE8);
            dset_a[i] = H5Dcreate (file, name.c_str(), fieldtype_a[i], space_a[i], H5P_DEFAULT, cparms, dataprop);
	}
        memspace_a[i] = H5Dget_space(dset_a[i]);
    }
}

int getEntryCount(int num,int call){
    int geno[4] = {3,6,10,15}; //possible number of genotypes for calls w/ 1,2,3,4 Alts resp.
    if(num==-1){
     	num = call-1;
    }else if(num==-2){
        num = geno[call-2]; //-1 for ref, -1 for 0-based idx
    }else if(num==-3){
	num = 2; //unknown count
    }//else default non-1 value
    return num;
}

void allocFieldVar(META_3 *format,int formcount,int CHUNK,int samcount,int ***&intvar,float ***&floatvar,char ***&charvar,char ****&stringvar,int *&varloc){
    int t1=0,t2=0,t3=0,t4=0;
    varloc = (int*)malloc(formcount*sizeof(int)); //for locating field in the multidim variable for each datatype

    for(int i=0;i<formcount;i++){
        if(format[i].num==1){
    	    if(format[i].type[0]=='0' || !strcmp(format[i].id,"GT")){ //integer/GT
	        varloc[i] = t1++;
    	    }else if(format[i].type[0]=='1'){ //float
             	varloc[i] = t2++;
    	    }else if(format[i].type[0]=='3'){ //character
            	varloc[i] = t3++;
    	    }else if(format[i].type[0]=='4'){ //string
	    	varloc[i] = t4++;
    	    }
        }else{
	    varloc[i] = t4++;
	}
    }
    if(t1>0){
	intvar = (int***)malloc(t1*sizeof(int**));
	intvar[0] = (int**)malloc(t1*CHUNK*sizeof(int*));
	intvar[0][0] = (int*)calloc(t1*CHUNK*samcount,sizeof(int));
	for (int i = 0; i < t1; i++){
            intvar[i] = intvar[0] + CHUNK*i;
            for (int j = 0; j < CHUNK; j++)
            {
            	intvar[i][j] = intvar[0][0] + CHUNK*samcount*i + samcount*j;
            }
    	}
    }
    if(t2>0){
	floatvar = (float***)malloc(t2*sizeof(float**));
	floatvar[0] = (float**)malloc(t2*CHUNK*sizeof(float*));
	floatvar[0][0] = (float*)calloc(t2*CHUNK*samcount,sizeof(float));
	for (int i = 0; i < t2; i++){
            floatvar[i] = floatvar[0] + CHUNK*i;
            for (int j = 0; j < CHUNK; j++)
            {
            	floatvar[i][j] = floatvar[0][0] + CHUNK*samcount*i + samcount*j;
            }
    	}
    }
    if(t3>0){
	charvar = (char***)malloc(t3*sizeof(char**));
	charvar[0] = (char**)malloc(t3*CHUNK*sizeof(char*));
	charvar[0][0] = (char*)calloc(t3*CHUNK*samcount,sizeof(char));
	for (int i = 0; i < t3; i++){
            charvar[i] = charvar[0] + CHUNK*i;
            for (int j = 0; j < CHUNK; j++)
            {
            	charvar[i][j] = charvar[0][0] + CHUNK*samcount*i + samcount*j;
            }
    	}
    }
    if(t4>0){
	stringvar = (char****)malloc(t4*sizeof(char***));
	stringvar[0] = (char***)malloc(t4*CHUNK*sizeof(char**));
	stringvar[0][0] = (char**)malloc(t4*CHUNK*samcount*sizeof(char*));
        stringvar[0][0][0] = (char*)calloc(t4*CHUNK*samcount*SIZE8,sizeof(char));
	for (int i = 0; i < t4; i++){
            stringvar[i] = stringvar[0] + CHUNK*i;
            for (int j = 0; j < CHUNK; j++)
            {
            	stringvar[i][j] = stringvar[0][0] + CHUNK*samcount*i + samcount*j;
		for(int k = 0; k< samcount;k++){
		    stringvar[i][j][k] = stringvar[0][0][0] + CHUNK*samcount*SIZE8*i +  samcount*j*SIZE8 + SIZE8*k;
		}
            }
    	}
    }

    //STRING
}

void freeFieldVar(int ***&intvar,float ***&floatvar,char ***&charvar,char ****&stringvar){
    if(intvar!=NULL){
    	free(intvar[0][0]);
    	free(intvar[0]);
    	free(intvar);
    }
    if(floatvar!=NULL){
    	free(floatvar[0][0]);
    	free(floatvar[0]);
    	free(floatvar);
    }
    if(charvar!=NULL){
    	free(charvar[0][0]);
    	free(charvar[0]);
    	free(charvar);
    }
    if(stringvar!=NULL){
    	free(stringvar[0][0][0]);
    	free(stringvar[0][0]);
    	free(stringvar[0]);
        free(stringvar);
    }

}

void closeIden(hid_t *&fieldtype_a,hid_t *&space_a,hid_t *&memspace_a,hid_t *&dset_a,int formcount){
    herr_t status;
    for(int i=0;i<formcount;i++){ //space_a is closed already
	status = H5Dclose(dset_a[i]);
    	status = H5Tclose(fieldtype_a[i]);
	status = H5Sclose(memspace_a[i]);
    }
    assert(status >=0);
    free(fieldtype_a);
    free(dset_a);
    free(space_a);
    free(memspace_a);
}

void parseFORMATfield(vector<vector<string> > token,queue<int> &tokenidx,unsigned int formbit,int ***&intvar,int *varloc,int fieldidx,int relidx,int samcount){
    if(formbit&(1<<fieldidx)){ //check if n-th bit/field is set
        int temp=tokenidx.front(); //index of the each FORMAT fields in a sample
	if(!token.empty()){
	    for(int x=0;x<samcount;x++){
	    	if(!token[x].empty()) 
		    intvar[varloc[fieldidx]][relidx][x] = (strcmp(token[x][temp].c_str(),"."))?atoi(token[x][temp].c_str()):-1; //-1 no value 
		else
		    intvar[varloc[fieldidx]][relidx][x] = -1; // not encoded -> ./.
            }
	}else{
	    for(int x=0;x<samcount;x++){
	    	intvar[varloc[fieldidx]][relidx][x] = -2; //indels/structural variants
            }
	}
	tokenidx.pop();
    }else{ //not a FORMAT field for current SNP
	for(int x=0;x<samcount;x++){
	    intvar[varloc[fieldidx]][relidx][x] = -3; //-3 not a field
        }
    } 
   
}

void parseFORMATfield(vector<vector<string> > token,queue<int> &tokenidx,unsigned int formbit,float ***&floatvar,int *varloc,int fieldidx,int relidx,int samcount){
    if(formbit&(1<<fieldidx)){ //check if n-th bit/field is set
        int temp=tokenidx.front(); //index of the each FORMAT fields in a sample
	if(!token.empty()){
	    for(int x=0;x<samcount;x++){
		if(!token[x].empty()) 
	    	    floatvar[varloc[fieldidx]][relidx][x] = (strcmp(token[x][temp].c_str(),"."))?atof(token[x][temp].c_str()):-1; //-1 no value 
		else
		    floatvar[varloc[fieldidx]][relidx][x] = -1;
            }
	}else{
	    for(int x=0;x<samcount;x++){
	    	floatvar[varloc[fieldidx]][relidx][x] = -2; //indels/structural variants
            }
	}
	tokenidx.pop();
    }else{ //not a FORMAT field for current SNP
	for(int x=0;x<samcount;x++){
	    floatvar[varloc[fieldidx]][relidx][x] = -3; //-3 not a field
        }
    } 
}

void parseFORMATfield(vector<vector<string> > token,queue<int> &tokenidx,unsigned int formbit,char ***&charvar,int *varloc,int fieldidx,int relidx,int samcount){
    if(formbit&(1<<fieldidx)){ //check if n-th bit/field is set
        int temp=tokenidx.front(); //index of the each FORMAT fields in a sample
	if(!token.empty()){
	    for(int x=0;x<samcount;x++){
		if(!token[x].empty()) 
	    	    charvar[varloc[fieldidx]][relidx][x] = (!token[x][temp].empty())?token[x][temp][0]:'.'; //-1 no value 
		else
		    charvar[varloc[fieldidx]][relidx][x] = '.';
            }
	}else{
	    for(int x=0;x<samcount;x++){
	    	charvar[varloc[fieldidx]][relidx][x] = '.'; //indels/structural variants
            }
	}
	tokenidx.pop();
    }else{ //not a FORMAT field for current SNP
	for(int x=0;x<samcount;x++){
	    charvar[varloc[fieldidx]][relidx][x] = 'X'; //-3 not a field
        }
    } 
}

void parseFORMATfield(vector<vector<string> > &token,queue<int> &tokenidx,unsigned int formbit,char ****&stringvar,int *varloc,int fieldidx,int relidx,int samcount){
   if(formbit&(1<<fieldidx)){ //check if n-th bit/field is set
	int temp=tokenidx.front(); //index of the each FORMAT fields in a sample
	if(!token.empty()){ 
	    for(int x=0;x<samcount;x++){
            	token[x][temp]+="\0";
		if(!token[x].empty()) 
	    	    strcpy(stringvar[varloc[fieldidx]][relidx][x],token[x][temp].c_str()); 
		else
		    strcpy(stringvar[varloc[fieldidx]][relidx][x],".\0");
		//cout << stringvar[varloc[fieldidx]][relidx][x];
		//stringvar[varloc[fieldidx]][relidx][x][token[x][temp].length()] = '\0';
            }
	}else{ //No entries
	    for(int x=0;x<samcount;x++){
	    	strcpy(stringvar[varloc[fieldidx]][relidx][x],".\0"); 
            }
	}
        tokenidx.pop(); 
    }else{ //not a FORMAT field for current SNP
	for(int x=0;x<samcount;x++){
	    	strcpy(stringvar[varloc[fieldidx]][relidx][x],".\0"); 
        }
    }
}



void parseGenotypes(hid_t file,string gpath,string linestream,int relidx,int samcount,int lastparsepos,vector<char> &callvec,map<string,int> states,META_3 *format,unsigned int formbit,map<string,int> formmap,int ***&intvar,float ***&floatvar,char ***&charvar,char ****&stringvar,int *varloc,queue<int> &tokenidx){ 
    int idx1 = lastparsepos, idx2 = lastparsepos,temp=0,counter=0; 
    int last = linestream.length(), formcount = formmap.size(),gtloc=0;  
    vector<vector<string> > token;
    string callstr;
    
    gtloc = varloc[formmap.find("GT")->second];
    if(callvec[0]=='X' || callvec[1]=='X' || callvec[0]=='.' || callvec[1]=='.'){ //REF or ALT
       //filter indels, structural variants,'.' calls
 	for(int i=0;i<samcount;i++){
            intvar[gtloc][relidx][i] = 25; //Indels/Incomplete Ref/Alt
        }
    }else{ //SNPs
	while(idx1<last){ 
             //get GT where 0-Ref 1..N-Alt
            if(linestream[idx1]=='.' || linestream[idx1+2]=='.'){
                intvar[gtloc][relidx][counter] = 26; //No entry
	    	if(linestream[idx1+3]==':') //if ./.:x:y:z
		    idx1+=2;
		else //if ./.
		    idx1++; 
		//idx1=linestream.find_first_of(":",idx1)-1; //some with '.' (NOT ./.) may still have FORMAT values
		
	    }else{
                //heterozygous
		temp=linestream[idx1]-48; //convert from char to int ('1'/'1' -> 1/1)  
	    	callstr = callvec[temp]; //get the call
            	idx1+=2;  
	    	temp=linestream[idx1]-48; //convert from char to int ('1'/'1' -> 1/1) 
	    	callstr += callvec[temp]; //get the call
                intvar[gtloc][relidx][counter] = states.find(callstr)->second; //save the call code
	    }
	  
            idx2 = linestream.find_first_of("\t",idx1); 
	    if(idx2==linestream.npos) idx2=last; //check if last sample
	    if(idx2>idx1) idx1+=2; 
            //cout << linestream.substr(idx1,idx2-idx1) << "\t";
 
            //PARSE FORMAT fields
            token.push_back( vector<string>() ); //token of strings
            while(idx1<idx2){ 
            	//other info for each sample aside from GT 
                temp = linestream.find_first_of(":\t",idx1+1);
		if(temp==linestream.npos) temp = last;
                token[counter].push_back(linestream.substr(idx1,temp-idx1));
                //cout << token[counter][token[counter].size()-1] << "|";
                idx1=temp+1;
	    }
	    //cout << "|";
            idx1 = idx2+1;
            counter++; 
    	}

        
    } 

    //store FORMAT fields
    for(int i=0;i<formcount;i++){
    	if(format[i].num==1){
            if(!strcmp(format[i].id,"GT")) continue; //already parsed
    	    if(format[i].type[0]=='0'){ //integer
	    	parseFORMATfield(token,tokenidx,formbit,intvar,varloc,i,relidx,samcount);
    	    }else if(format[i].type[0]=='1'){ //float
             	parseFORMATfield(token,tokenidx,formbit,floatvar,varloc,i,relidx,samcount);
    	    }else if(format[i].type[0]=='3'){ //character
            	parseFORMATfield(token,tokenidx,formbit,charvar,varloc,i,relidx,samcount);
    	    }else if(format[i].type[0]=='4'){ //string
	    	parseFORMATfield(token,tokenidx,formbit,stringvar,varloc,i,relidx,samcount);
    	    }
     	}else{
	     	parseFORMATfield(token,tokenidx,formbit,stringvar,varloc,i,relidx,samcount);
	}
    }
    	
    if(!token.empty()){	//free
    	for(int x=0;x<samcount;x++) 
	   token[x].clear();
    }
    callvec.clear(); 
    token.clear();
    
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

int buildPosIndex(string datafile,string inipath){
    string varpath,indexfile,varname;
    FQ::FileFormat model = FQ::FQ_HDF5; 
    indexfile = datafile+".index";
    char binning[]= "precision=2";
    ibis::gParameters().add(FQ_REPORT_STATISTIC, "false");
    //ibis::gParameters().add("fileManager.maxBytes", "2GB");
    
    IndexBuilder* indexBuilder = new IndexBuilder(datafile, model, indexfile, 0,"",""); 
    if (!indexBuilder->isValid()) {
	cout << "ERROR: Failed to initialize the IndexBuilder object for file.\n";
	delete(indexBuilder);
	return 1;
    }
    varname = "/pos[,0]";
    indexBuilder->buildIndexes(binning, inipath, varname.c_str());
    varname = "/pos[,1]";
    indexBuilder->buildIndexes(binning, inipath, varname.c_str());
    delete(indexBuilder);
}


int buildFieldIndex(string datafile,META_3 *format,int formcount,string inipath){
    string varpath,indexfile,varname;
    FQ::FileFormat model = FQ::FQ_HDF5; 
    indexfile = datafile+".index";
    char binning[]= "precision=2";
    ibis::gParameters().add(FQ_REPORT_STATISTIC, "false");
    //ibis::gParameters().add("fileManager.maxBytes", "2GB");
    
    IndexBuilder* indexBuilder = new IndexBuilder(datafile, model, indexfile, 0,"",""); 
    if (!indexBuilder->isValid()) {
	cout << "ERROR: Failed to initialize the IndexBuilder object for file.\n";
	delete(indexBuilder);
	return 1;
    }
    delete(indexBuilder);
    varpath = inipath + "/FORMATfields";
    for(int i=0;i<formcount;i++){
        if(format[i].num==1){//only single-value field can be indexed
            if(!strcmp(format[i].id,"GT") || format[i].type[0]=='3' || format[i].type[0]=='4') 	
		continue; //no indexing for GT/char/string
	    varname = format[i].id;
            IndexBuilder* indexBuilder = new IndexBuilder(datafile, model, indexfile,0,"",""); 
	    indexBuilder->buildIndexes(binning, varpath, varname.c_str());
	    delete(indexBuilder);
     	}
    }
}

int writeVCF(string inputfile,string varPath, string varName, string datafile){
    ibis::horometer timer1,timer2;
    timer1.start(); 
    const rlim_t STACK_SIZE = 1000*1024*1024; 
    struct rlimit rl;
    rl.rlim_cur = STACK_SIZE;
    int ret = setrlimit(RLIMIT_STACK,&rl);

    ifstream fp(inputfile.c_str());
    if (!fp.is_open()) { //check input file
	cout << "ERROR: Failed to open the input file \"" << inputfile.c_str() << "\"" << endl;
	cout << "REPORT: Failed to complete writing the data" << endl;
	return 1;
    }
    
    char* line=NULL;
    map<string,int> headmap;
    map<string,int> states; //allele states
    map<string,int> contigmap;
    map<string,int> infomap;
    map<string,int> formmap;
    vector<char> callvec; //vec for REF & ALT calls
    vector<string> meta1_val;
    string linestream,inipath,gpath1,gpath2;
    hid_t file,group,group_sub1,group_sub2;
    META_2* contig=NULL;
    MISC* misc = NULL;
    META_3* info=NULL, *format = NULL;
    int contigcount=0,infocount=0,formcount=0,samcount=0, indelcount = 0;
    int **CHR_POS=NULL;
    
    /*create new file*/
    inipath=varPath+varName;
    hid_t fapl;
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_libver_bounds(fapl, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    file = H5Fcreate(datafile.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,fapl); 
    group = H5Gcreate (file, inipath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    gpath1 = inipath + "/meta"; 
    group_sub1 = H5Gcreate (file, gpath1.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    gpath2 = inipath + "/FORMATfields"; 
    group_sub2 = H5Gcreate (file, gpath2.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i=0;getline(fp,linestream);i++){
        if(linestream.substr(0,2) == "##"){
            if(linestream.substr(0,8) == "##contig"){
	    	parseHeader2(linestream,contig,++contigcount);
            }else if(linestream.substr(0,6) == "##INFO"){ 
	    	parseHeaderInfoFormat(linestream,info,++infocount);
            }else if(linestream.substr(0,8) == "##FORMAT"){
                parseHeaderInfoFormat(linestream,format,++formcount);
 	    }else{
 	    	parseHeader1(linestream,headmap,meta1_val); 
	    }
	}else if(linestream.substr(0,6) == "#CHROM"){ 
	    loadHeader1(file,headmap,meta1_val,gpath1); 
	    if(contigcount){ loadHeader2(file,contig,contigcount,contigmap,gpath1);}
	    loadHeaderInfoFormat(file,false,info,infocount,infomap,gpath1); //INFO fields 
    	    loadHeaderInfoFormat(file,true,format,formcount,formmap,gpath1); //FORMAT fields 
            samcount = loadSampleNames(file,linestream,inipath);
	    break; 
	}else{
            cout << "ERROR: Failed to write the input file. Unknown format of header entry. \"" << "\"" << endl;
	    cout << "REPORT: Failed to complete writing the data" << endl;
	    return 1; /*finished reading header*/
	}
    } 
    
    int counter1=0,counter2=0,lastparsepos; 
    int *varloc; //location of FORMAT fields in multidimensional variables(by type) below
    queue<int> tokenidx; //index of each FORMAT field values in the tokenized data
			 //this is necessary since FORMAT fields/SNP may not be in order 
       //variables to contain parsed fields(single-value)
    int ***intvar=NULL; //2D to accomodate fields with same datatype
    float ***floatvar=NULL;
    char ***charvar=NULL;
      //variables to contain parsed fields(multi-value/string type)
    char **** stringvar=NULL; 
    
    hid_t memtype1,space1,memspace1_1,memspace1_2,dset1_1,dset1_2,cparms1,dataprop1; //variables for MISC
    hid_t *fieldtype_a, *space_a, *memspace_a,*dset_a,cparms2,dataprop2; //variables for FORMAT field
    herr_t status;
    hsize_t dim1_1[1]={CHUNKSIZE1}; 
    hsize_t dim1_2[2]={CHUNKSIZE1,2}; 
    hsize_t dim2[2]={CHUNKSIZE2_1,samcount};
    hsize_t maxdim1_1[1] = {H5S_UNLIMITED};
    hsize_t maxdim1_2[2] = {H5S_UNLIMITED,2};
    hsize_t maxdim2[2] = {H5S_UNLIMITED,samcount};
    hssize_t offset1_1[1]={0};
    hssize_t offset1_2[2]={0,0};
    hssize_t offset2[2]={0,0};
    hsize_t count1_1[1]={CHUNKSIZE1};
    hsize_t count1_2[2]={CHUNKSIZE1,2};
    hsize_t count2[2]={CHUNKSIZE2_1,samcount};
    hsize_t newsize1_1[1]={CHUNKSIZE1};
    hsize_t newsize1_2[2]={CHUNKSIZE1,2};
    hsize_t newsize2[2]={CHUNKSIZE2_1,samcount};

    //allocate space for FORMAT field container/var
    allocFieldVar(format,formcount,CHUNKSIZE2_1,samcount,intvar,floatvar,charvar,stringvar,varloc);

    misc = (MISC*)malloc(CHUNKSIZE1*sizeof(MISC));
    //separate Chrom/Contig index and Position to enable indexing using FastQuery
    CHR_POS = (int**)malloc(CHUNKSIZE1*sizeof(int*));
    CHR_POS[0] = (int*)calloc(CHUNKSIZE1*2,sizeof(int));
    for(int i=0;i<CHUNKSIZE1;i++){ CHR_POS[i]=CHR_POS[0]+i*2;} 

    if(misc==NULL){
 	cout << "ERROR: Insufficient memory. Adjust the chunksize." << endl;
        cout << "REPORT: Failed to complete writing the data" << endl;
	return 1;
    }
    setAlleleStates(states); 	
    
    for(int x=0;fp!=NULL;x++){ 
        getline(fp,linestream); 
	if(fp!=NULL){ 
            lastparsepos = parseMisc(linestream,misc,CHR_POS,infomap,formmap,contig, contigcount,contigmap,callvec,counter1, indelcount,tokenidx); 
            counter1++; 
            parseGenotypes(file,gpath2,linestream,counter2,samcount,lastparsepos+1,callvec,states,format,misc[counter1-1].format,formmap,intvar,floatvar,charvar, stringvar,varloc,tokenidx);
            counter2++; 
        } 
        if(counter1==CHUNKSIZE1 || (fp==NULL && counter1>0)){
            if(x+1==CHUNKSIZE1){ //first slab
                //set compound datatype and spaces
		setH5MiscFormat(file,memtype1,space1,dset1_1,cparms1,dataprop1,CHUNKSIZE1,inipath);
		setH5CHRPOSFormat(file,space1,dset1_2,cparms1,dataprop1,CHUNKSIZE1,inipath);
                memspace1_1 = H5Dget_space(dset1_1);
                memspace1_2 = H5Dget_space(dset1_2); 
                //write first slab
            	status = H5Dwrite(dset1_1, memtype1, H5S_ALL, H5S_ALL, H5P_DEFAULT, misc);  
		status = H5Dwrite(dset1_2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &CHR_POS[0][0]); 
	    }else{  
		offset1_1[0]=newsize1_1[0];
                offset1_2[0]=newsize1_2[0];
          	//extend data later for the next slab
                if(fp==NULL && counter1>0){
 		     count1_1[0]= count1_2[0] = counter1;
		     newsize1_1[0]=newsize1_2[0]=newsize1_1[0]+counter1;
		     dim1_1[0]=dim1_2[0]=counter1;
    		     memspace1_1 = H5Screate_simple(1,dim1_1,maxdim1_1); 
		     memspace1_2 = H5Screate_simple(2,dim1_2,maxdim1_2); 
		}else{
		     newsize1_1[0]=newsize1_2[0]=newsize1_1[0]+counter1; 
		}
		status = H5Dset_extent(dset1_1,newsize1_1);
                status = H5Dset_extent(dset1_2,newsize1_2);
	 	space1 = H5Dget_space(dset1_1);
                status = H5Sselect_hyperslab(space1, H5S_SELECT_SET,
 			(const hsize_t*)offset1_1,NULL, count1_1, NULL);
		status = H5Dwrite(dset1_1,memtype1,memspace1_1,space1,H5P_DEFAULT,misc); 
		status = H5Sclose(space1); 
		space1 = H5Dget_space(dset1_2);
                status = H5Sselect_hyperslab(space1, H5S_SELECT_SET,
 			(const hsize_t*)offset1_2,NULL, count1_2, NULL);
		status = H5Dwrite(dset1_2,H5T_NATIVE_INT,memspace1_2,space1,H5P_DEFAULT,&CHR_POS[0][0]); 
                status = H5Sclose(space1); 
	    } 
	    cout << "MISC-Chunk" << (x+1)/CHUNKSIZE1 <<"\n";
            counter1=0;
        }
	
        if(counter2==CHUNKSIZE2_1 || (fp==NULL && counter2>0)){
            if(x+1==CHUNKSIZE2_1){ //first slab
                //set GT and FORMAT fields datatype and spaces
                setH5FORMATfield(file,fieldtype_a,space_a,memspace_a, dset_a,cparms2,dataprop2,CHUNKSIZE2_1,CHUNKSIZE2_2,samcount,inipath,format,formcount);
             
		//write first slab for other FORMAT fields
		for(int i=0;i<formcount;i++){ 
  		    if(format[i].num==1){
    	    	    	if(format[i].type[0]=='0' || !strcmp(format[i].id,"GT")){ //integer/GT
	        	    status = H5Dwrite(dset_a[i],fieldtype_a[i], H5S_ALL, H5S_ALL, H5P_DEFAULT, &intvar[varloc[i]][0][0]); 
    	    	    	}else if(format[i].type[0]=='1'){ //float
             		    status = H5Dwrite(dset_a[i],fieldtype_a[i], H5S_ALL, H5S_ALL, H5P_DEFAULT, &floatvar[varloc[i]][0][0]);
    	    	    	}else if(format[i].type[0]=='3'){ //character
            		    status = H5Dwrite(dset_a[i],fieldtype_a[i], H5S_ALL, H5S_ALL, H5P_DEFAULT, &charvar[varloc[i]][0][0]);
    	    	    	}else if(format[i].type[0]=='4'){ //string
	    		    status = H5Dwrite(dset_a[i],fieldtype_a[i], H5S_ALL, H5S_ALL, H5P_DEFAULT, stringvar[varloc[i]][0][0]); 
    	    	    	}
        	    }else{
	    		status = H5Dwrite(dset_a[i],fieldtype_a[i], H5S_ALL, H5S_ALL, H5P_DEFAULT, stringvar[varloc[i]][0][0]); 
		    }
	        }
	    }else{ 
		offset2[0]=newsize2[0];
          	//extend data later for the 2nd to the last slab
                if(fp==NULL && counter2>0){
 		    count2[0] = counter2;
		    newsize2[0]=newsize2[0]+counter2;
		    dim2[0]=counter2;
		    for(int i=0;i<formcount;i++){ 
  			memspace_a[i] = H5Screate_simple(2,dim2,maxdim2); 
		    }
		}else{
		     newsize2[0]=newsize2[0]+counter2; 
		}
		  
                for(int i=0;i<formcount;i++){ 
		    status = H5Dset_extent(dset_a[i],newsize2);
		    space_a[i] = H5Dget_space(dset_a[i]);
		    status = H5Sselect_hyperslab(space_a[i], H5S_SELECT_SET,
 			(const hsize_t*)offset2,NULL, count2, NULL);
		    if(format[i].num==1){
    	    	    	if(format[i].type[0]=='0' || !strcmp(format[i].id,"GT")){ //integer/GT
	        	    status = H5Dwrite(dset_a[i],fieldtype_a[i], memspace_a[i], space_a[i], H5P_DEFAULT, &intvar[varloc[i]][0][0]);  
    	    	    	}else if(format[i].type[0]=='1'){ //float
             		    status = H5Dwrite(dset_a[i],fieldtype_a[i], memspace_a[i], space_a[i], H5P_DEFAULT, &floatvar[varloc[i]][0][0]);  
    	    	    	}else if(format[i].type[0]=='3'){ //character
            		    status = H5Dwrite(dset_a[i],fieldtype_a[i], memspace_a[i], space_a[i], H5P_DEFAULT, &charvar[varloc[i]][0][0]);  
    	    	    	}else if(format[i].type[0]=='4'){ //string
	    		    status = H5Dwrite(dset_a[i],fieldtype_a[i], memspace_a[i], space_a[i], H5P_DEFAULT, stringvar[varloc[i]][0][0]);
    	    	    	}
        	    }else{
	    		status = H5Dwrite(dset_a[i],fieldtype_a[i], memspace_a[i], space_a[i], H5P_DEFAULT, stringvar[varloc[i]][0][0]);  
		    }
		    status = H5Sclose(space_a[i]); 
		}
	    } 
	    cout << "\nCall-Chunk" << (x+1)/CHUNKSIZE2_1 <<"\n";
            counter2=0;
        }
        
    }
    int row = newsize2[0]; //total SNPs
    //load CHROM IDs
    if(contigcount==0){ loadHeader2(file,contig,contigcount,contigmap,gpath1);}

    cout << "Indel/Structural Variant Count:" << indelcount << "\n";
    //free value pointers and variables
    contigmap.clear(); 
    infomap.clear();
    formmap.clear();
    states.clear();
    status = H5Dclose(dset1_1);
    status = H5Dclose(dset1_2);
    status = H5Tclose(memtype1);
    status = H5Sclose(memspace1_1);
    status = H5Sclose(memspace1_2);
    status = H5Pclose(cparms1);
    status= H5Pclose (dataprop1);
    status = H5Pclose(cparms2);
    status= H5Pclose (dataprop2); 
    status = H5Gclose(group);
    status = H5Gclose(group_sub1);
    status = H5Gclose(group_sub2);
    freeFieldVar(intvar,floatvar,charvar,stringvar);
    closeIden(fieldtype_a,space_a,memspace_a,dset_a,formcount);
    H5Fclose(file); //close file

    timer1.stop();
    printf("REPORT: Successfully completed loading data.\n Total time elapsed:%f\n", timer1.realTime());
    timer2.start();
    free(misc);
    free(CHR_POS[0]);
    //build index
    buildPosIndex(datafile,inipath);
    buildFieldIndex(datafile,format,formcount,inipath);
    timer2.stop();
    printf("REPORT: Successfully completed indexing data.\n Total time elapsed:%f\n", timer2.realTime());
    free(info);
    free(format);
    free(varloc);
        
    return 0;
}
