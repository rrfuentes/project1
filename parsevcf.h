/*
 *Last Update: Aug. 16, 2013
 *Author: Roven Rommel B. Fuentes
 *TT-Chang Genetic Resources Center, International Rice Research Institute
 *
*/

#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <assert.h>
#include <sys/resource.h>
#include <vector>
#include <math.h>

#define CHUNKSIZE1 5
#define CHUNKSIZE2_1 5 
#define CHUNKSIZE2_2 2 
#define SNP_CHUNK_CACHE 268435456 /*250MB*/
#define SIZE1 20
#define SIZE2 20
#define SIZE3 20
#define SIZE4 20
#define SIZE5 15

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

void loadHeader1(map<string,string> &headmap,hid_t file,string gpath){
    /*relocate header to META_1 before loading to HDF5*/
    int x=0,count=headmap.size(); 
    META_1 header[count];
    hid_t memtype,space,dset;
    hid_t t1,t2;
    hsize_t dim[1]={count};
    herr_t status;
    string name = gpath + "/meta"; 
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

void loadHeader2(META_2 *&contig,hid_t file,int count,map<string,int> &contigmap,string gpath){
    hid_t memtype,space,dset,t1;
    hsize_t dim[1]={count};
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

void loadHeaderInfoFormat(bool flag,META_3 *&field,int count,map<string,int> &map, hid_t file,string gpath){
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
    H5Tset_size(t2,9);
    t3 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t3,H5T_VARIABLE);
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
        field[x].type[1]='\0';
        free(field[x].desc);
    }
    cout << "\n";
}

int loadSampleNames(string linestream, hid_t file,string inipath){
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
    string name= inipath + "/samples";
    space = H5Screate_simple(1,dim,NULL);
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,H5T_VARIABLE);
    dset = H5Dcreate (file, name.c_str(), t1, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
    idx1+=3; //skip GT for flag
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

int parseMisc(string linestream, MISC *&misc,map<string,int> infomap, map<string,int> formmap, map<string,int> contigmap,vector<char> &callvec, int snpidx,int &indelcount){
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


void setH5MiscFormat(hid_t file,hid_t &memtype,hid_t &space,hid_t &dset,hid_t &cparms,hid_t &dataprop,int chunk,string inipath){
    hid_t t1,t2,t3;
    herr_t status;
    hsize_t dim[1]={chunk};
    hsize_t maxdim[1] = {H5S_UNLIMITED};
    hsize_t chkdim[1] = {chunk};
    hssize_t offset[1] = {0};

    string name = inipath + "/misc";
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
    assert(status >=0);
    dset = H5Dcreate (file, name.c_str(), memtype, space, H5P_DEFAULT, cparms, dataprop);
   
   /*free value pointers and variables*/
    status = H5Tclose(t1);
    status = H5Tclose(t2);
    status = H5Tclose(t3);
    assert(status>=0);
}

void setH5CallFormat(hid_t file, hid_t &space,hid_t &dset,hid_t &cparms,hid_t &dataprop,int size1,int size2,string inipath){
    herr_t status;
    hsize_t dim[2]={size1,size2};
    hsize_t maxdim[2] = {H5S_UNLIMITED,size2};
    hsize_t chkdim[2] = {size1,CHUNKSIZE2_2};
    hssize_t offset[2] = {0,0};

    string name = inipath + "/call";
    space = H5Screate_simple(2,dim,maxdim); //create data space
    cparms = H5Pcreate(H5P_DATASET_CREATE); //create chunk 
    status = H5Pset_szip(cparms,H5_SZIP_NN_OPTION_MASK,8); //compression
    status = H5Pset_chunk(cparms,2,chkdim);
    
    //create dataset access property list for MISC dataset
    dataprop = H5Pcreate(H5P_DATASET_ACCESS);
    status = H5Pset_chunk_cache(dataprop,H5D_CHUNK_CACHE_NSLOTS_DEFAULT,
             SNP_CHUNK_CACHE*2,H5D_CHUNK_CACHE_W0_DEFAULT); //set snp chunk size
    assert(status >=0);
    dset = H5Dcreate (file, name.c_str(), H5T_NATIVE_INT, space, H5P_DEFAULT, cparms, dataprop);
    
}

int getBitPos(int b){
    return log2(b&-b);
}

int getEntryCount(int num,int call){
    int geno[4] = {3,6,10,15}; //possible number of genotypes for calls w/ 1,2,3,4 Alts resp.
    if(num==-1){
     	num = call-1;
    }else if(num==-2){
        num = geno[call-2]; //-1 for ref, -1 for 0-based idx
    }else if(num==-3){
	num = 1; //unknown count
    }//else default non-1 value
    return num;
}

void setType(META_3 *format,unsigned int formbit,hid_t *&fieldtype,hid_t &memtype,vector<pair<int,int> > &bitpos,int call){
    int x,bitcount,size=0;
    hid_t t1;
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,H5T_VARIABLE);
    hid_t *temptype = (hid_t*)malloc(32*sizeof(hid_t));
    fieldtype = (hid_t*)malloc(32*sizeof(hid_t)); //max of 32 fields

    for(x=0; formbit; formbit&=formbit-1,x++){
	bitpos.push_back(make_pair(getBitPos(formbit&(-formbit)),0)); //get the rightmost bit
	bitpos[x].second = getEntryCount(format[bitpos[x].first].num,call); //save the entry count
	if(bitpos[x].second){
	    if(format[bitpos[x].first].type[0]=='0'){ //integer
            	temptype[x] = H5Tcopy(H5T_NATIVE_INT);
    	    }else if(format[bitpos[x].first].type[0]=='1'){ //float
             	temptype[x] = H5Tcopy(H5T_NATIVE_FLOAT);
    	    }else if(format[bitpos[x].first].type[0]=='2'){ //flag
	    	continue;
    	    }else if(format[bitpos[x].first].type[0]=='3'){ //character
            	temptype[x] = H5Tcopy(H5T_NATIVE_CHAR);
    	    }else if(format[bitpos[x].first].type[0]=='4'){ //string
	    	temptype[x] = H5Tcopy(t1);
	    }
	}else{
	    hsize_t dim[1]={bitpos[x].second};
            if(format[bitpos[x].first].type[0]=='0'){ //integer
            	temptype[x] = H5Tarray_create(H5T_NATIVE_INT,1,dim);
    	    }else if(format[bitpos[x].first].type[0]=='1'){ //float
             	temptype[x] = H5Tarray_create(H5T_NATIVE_FLOAT,1,dim);
    	    }else if(format[bitpos[x].first].type[0]=='2'){ //flag
	    	continue;
    	    }else if(format[bitpos[x].first].type[0]=='3'){ //character
            	temptype[x] = H5Tarray_create(H5T_NATIVE_CHAR,1,dim);
    	    }else if(format[bitpos[x].first].type[0]=='4'){ //string
	    	temptype[x] = H5Tarray_create(t1,1,dim);;
	    }
        }
        size+=H5Tget_size(temptype[x]);
        cout << format[bitpos[x].first].id << "-";
    } 

    bitcount=x;
    memtype = H5Tcreate(H5T_COMPOUND,size);
    int offset=0; cout << "here";
    for(x=0;x<bitcount;x++){
        H5Tinsert(memtype,format[bitpos[x].first].id,offset,temptype[x]); //for whole compound type
	offset+=H5Tget_size(temptype[x]);
        fieldtype[x] = H5Tcreate(H5T_COMPOUND,sizeof(temptype[x])); //single field for write-by-field
        H5Tinsert(fieldtype[x],format[bitpos[x].first].id,0,temptype[x]); //insert field with same name
        H5Tclose(temptype[x]);
    }
}

int allocVar(vector<vector<string> > token,int idx,int **&intvar,int num,int samcount){
    int count = 0;   
    int geno[4] = {3,6,10,15}; //possible number of genotypes for calls w/ 1,2,3,4 Alts resp.
    if(num==1){
    	intvar = (int**)malloc(sizeof(int*));
	intvar[0] = (int*)malloc(samcount*sizeof(int));
        for(int x=0;x<samcount;x++){
	    intvar[0][x] = (token[x].size())?atoi(token[x][idx].c_str()):0;
	    cout << intvar[0][x] << "+";
	}
        return 1; //return num of values in the field
    }else{
	intvar = (int**)malloc(samcount*sizeof(int*));
        intvar[0] = (int*)malloc(samcount*count*sizeof(int));
	for(int x=0;x<samcount;x++){
	    intvar[x] = intvar[0]+x*count;
	    int off1=0,off2=0;
	    for(int y=0;y<count;y++){ 
	        if(token[x].size()){ 
                    off2 = (token[x][idx]).find_first_of(",",off1);
		    intvar[x][y] = atoi((token[x][idx].substr(off1,off2-off1)).c_str());
                    cout << intvar[x][y] <<"*";
		    off1=off2+1;
		}else{
		    intvar[x][y] = 0;
		}
	    }
	}
	return count;
    }
}

void loadFormatFields(hid_t dset,hid_t fieldtype,int **intvar,int call,int samcount,int retcount){
    herr_t status;
    if(call==1){
    	status = H5Dwrite(dset, fieldtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, intvar[0]);
    }else{
	status = H5Dwrite(dset, fieldtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &intvar[0]);
    }
    assert(status >=0);
}

void parseGenotypes(hid_t file,string gpath,string linestream,int **&call,int genidx,int snpidx,int samcount,int offset,vector<char> &callvec,map<string,int> states,META_3 *format,map<string,int> formmap,unsigned int formbit){
    int idx1 = offset, idx2 = offset,temp=0,counter=0; 
    int last = linestream.length(); 
    //variables to contain parsed fields(single/multi-value)
    int **intvar; //3D to accomodate fields using same datatype
    float **floatvar;
    char **charvar;
    vector<vector<string> > stringvar;
    vector<pair<int,int> > bitpos;
    vector<vector<string> > token;
    string callstr;
    herr_t status;
    hid_t memtype,space,dset,xfer_id;
    hid_t* fieldtype;
    hsize_t dim[1]={samcount};
    ostringstream num;
    num << genidx;
    string name = gpath + "/" + num.str();
    space = H5Screate_simple(1,dim,NULL);

    if(callvec[0]=='X' || callvec[1]=='X' || callvec[0]=='.' || callvec[1]=='.'){ //filter indels, structural variants,'.' calls
 	while(idx1<last){
            idx2 = linestream.find_first_of("\t",idx1+2); 
	    if(idx2==linestream.npos) idx2=last;
	    call[snpidx][counter] = 25;
	    idx1 = idx2+1;
            counter++;
        }
    }else{ //SNPs
	while(idx1<last){
             //get Ref where 0-Ref 1..N-Alt
            if(linestream[idx1]=='.' || linestream[idx1+2]=='.'){
		call[snpidx][counter] = 25;
		idx1+=2;
	    }else{
		temp=linestream[idx1]-48; //convert from char to int  
	    	callstr = callvec[temp]; //get the call
            	//get Alt 
            	idx1+=2;  
	    	temp=linestream[idx1]-48; //convert from char to int  
	    	callstr += callvec[temp]; //get the call
             	call[snpidx][counter] = states.find(callstr)->second; //save the call code
	    }
	    
	    token.push_back( vector<string>() );
            idx2 = linestream.find_first_of("\t",idx1); 
	    if(idx2==linestream.npos) idx2=last;
	    if(idx2>idx1) idx1+=2;
            //cout << linestream.substr(idx1,idx2-idx1) << "\t";
            int fieldctr=0;
            while(idx1<idx2){
            	//other info for each sample aside of GT
                temp = linestream.find_first_of(":\t",idx1+1); //skip : after GT 
                token[counter].push_back(linestream.substr(idx1,temp-idx1));
                //cout << token[counter][token[counter].size()-1] << "-";
                idx1=temp+1;
	    }
	    //cout << "|";
            idx1 = idx2+1;
            counter++;
    	}
        cout <<"\n";
    }

    //set datatype for each SNP
    setType(format,formbit,fieldtype,memtype,bitpos,callvec.size()); 
    //create dataset
    dset = H5Dcreate (file,name.c_str(), memtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    int typecount = bitpos.size();
    int retcount=0;
    
    for(int x=0;x<typecount; x++){
	if(format[bitpos[x].first].type[0]=='0'){ //integer
	    retcount = allocVar(token,x,intvar,bitpos[x].second,samcount);  
	    loadFormatFields(dset,fieldtype[x],intvar,bitpos[x].second,samcount,retcount);
    	}else if(format[bitpos[x].first].type[0]=='1'){ //float
	    	    
    	}else if(format[bitpos[x].first].type[0]=='2'){ //flag
	    	    
    	}else if(format[bitpos[x].first].type[0]=='3'){ //character
	    	    
    	}else if(format[bitpos[x].first].type[0]=='4'){ //string
	    	   
	}
	H5Tclose(fieldtype[x]);
    }
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(memtype);
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
