/*
 *Last update: July 19, 2013
 *Author: Roven Rommel B. Fuentes
 *TT-Change Genetic Resources Center, International Rice Research Institute
 *
 * Load hapmap data to HDF5.
 *load sample names
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

#define CHUNKSIZE1 300 //500 //10
#define CHUNKSIZE2_1 1200 //500 //10
#define CHUNKSIZE2_2 100 //2 //3
#define FIXCOL 11
#define SNP_CHUNK_CACHE 1048576000 /*536870912 /*500MB*/
#define SIZE1 15
#define SIZE2 3
#define SIZE3 5
#define SIZE4 8
#define SIZE5 20
#define SIZE6 20
#define SIZE7 20
#define SIZE8 20
#define SIZE9 3
#define SIZE10 30

using namespace std;

static char *fileModel = 0;

typedef struct{
    char rs[SIZE1];
    char snpallele[SIZE2];
    char chrom[SIZE3];
    int pos;
    char strand;
    char build[SIZE4];
    char center[SIZE5];
    char prot[SIZE6];
    char assay[SIZE7];
    char panel[SIZE8];
    char cq[SIZE9];
}MISC;

int loadSampleNames(string linestream, hid_t file,string headervar){
    char** samples=NULL;
    int idx1=0,idx2=0,x=0,count=0;
    for(x=0;x<11;x++){ //ignore fixed table headername
        idx1 = linestream.find_first_of("\t",idx1+1); 
    } 
    for(x=0;idx1<linestream.npos;x++){ 
	samples = (char**)realloc(samples,(x+1)*sizeof(char*));
        samples[x] = (char*)malloc(SIZE10*sizeof(char));
        idx2 = linestream.find_first_of("\t",idx1+1); 
        if(idx2==linestream.npos){
	    strcpy(samples[x],(linestream.substr(idx1+1,linestream.length()-idx1)).c_str());
	}else{
	    strcpy(samples[x],(linestream.substr(idx1+1,idx2-idx1-1)).c_str());
	}
        idx1=idx2; 
    } 
    count=x;
    
    hid_t memtype,space,dset,t1;
    hsize_t dim[1]={count};
    herr_t status;
    space = H5Screate_simple(1,dim,NULL);
    t1 = H5Tcopy(H5T_C_S1);
    H5Tset_size(t1,H5T_VARIABLE);
    dset = H5Dcreate (file, headervar.c_str(), t1, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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

int writeHapmap(string inputfile,string varPath, string varName, string datafile){
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

    map<string, int> Calls;

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
    
    /*create new file*/
    string miscvar, snpvar,inipath,headervar,linestream;
    hid_t file,group;
    file = H5Fcreate(datafile.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

    //create group
    inipath = varPath + varName;
    group = H5Gcreate (file, inipath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    miscvar = inipath + "/misc";
    snpvar = inipath + "/snp";
    headervar = inipath + "/samples"; 

    /*create another dataset for genoytpe field names*/
    getline(fp,linestream);
    int samcount = loadSampleNames(linestream,file,headervar);

    char *tok=NULL;
    int total=0;
    hid_t memtype1, memtype2, memspace1,memspace2, space1, space2, dset1, dset2, dataprop; /*handles*/
    hid_t t1,t2,t3,t4,t5,t6,t7,t8,t9;
    hid_t cparms1,cparms2;
    herr_t status;
    hsize_t dim1[1]={CHUNKSIZE1}; /*dataspace dimension for misc data*/
    hsize_t dim2[2]={CHUNKSIZE2_1,samcount}; /*dataspace dimension for genotype data*/
    hsize_t maxdim1[1]={H5S_UNLIMITED};
    hsize_t maxdim2[2]={H5S_UNLIMITED,samcount};
    hsize_t chkdim1[1]={CHUNKSIZE1}; /*misc chunk size*/
    hsize_t chkdim2[2]={CHUNKSIZE2_1,CHUNKSIZE2_2}; /*genotype chunk size*/
    hssize_t offset1[1]={0}; /*offset for slabs*/
    hssize_t offset2[2]={0,0};
    hsize_t newsize1[1]={CHUNKSIZE1}; /*iterate offset for slabs*/
    hsize_t newsize2[2]={CHUNKSIZE2_1,samcount};
    hsize_t count1[1]={CHUNKSIZE1}; /*param fo hyperslab*/
    hsize_t count2[2]={CHUNKSIZE2_1,samcount};
    
    if(CHUNKSIZE1<=1 && CHUNKSIZE2_1<=1){
	cout << "ERROR: Invalid chunk size." << endl;
	cout << "REPORT: Failed to complete writing the data" << endl;
	return 1;
    }

    /*create data space*/
    space1 = H5Screate_simple(1,dim1,maxdim1); 
    space2 = H5Screate_simple(2,dim2,maxdim2);
    /*create chunk*/
    cparms1 = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_chunk(cparms1,1,chkdim1);
    cparms2 = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_chunk(cparms2,2,chkdim2);

    t1 = H5Tcopy(H5T_C_S1);
    t2 = H5Tcopy(H5T_C_S1);
    t3 = H5Tcopy(H5T_C_S1);
    t4 = H5Tcopy(H5T_C_S1);
    t5 = H5Tcopy(H5T_C_S1);
    t6 = H5Tcopy(H5T_C_S1);
    t7 = H5Tcopy(H5T_C_S1);
    t8 = H5Tcopy(H5T_C_S1);
    t9 = H5Tcopy(H5T_C_S1);
    /*create compound datatype for fixed columns*/
    memtype1 = H5Tcreate(H5T_COMPOUND,sizeof(MISC));
    status = H5Tset_size(t1,SIZE1);
    status = H5Tinsert(memtype1,"rs#",HOFFSET(MISC,rs),t1);
     assert (status >= 0);
    status = H5Tset_size(t2,SIZE2);
    status = H5Tinsert(memtype1,"snpallele",HOFFSET(MISC,snpallele),t2);
    assert (status >= 0);
    status = H5Tset_size(t3,SIZE3);
    status = H5Tinsert(memtype1,"chrom",HOFFSET(MISC,chrom),t3);
    assert (status >= 0);
    status = H5Tinsert(memtype1,"pos",HOFFSET(MISC,pos),H5T_NATIVE_INT);
    assert (status >= 0);
    status = H5Tinsert(memtype1,"strand",HOFFSET(MISC,strand),H5T_NATIVE_CHAR);
    assert (status >= 0);
    status = H5Tset_size(t4,SIZE4);
    status = H5Tinsert(memtype1,"build",HOFFSET(MISC,build),t4);
    assert (status >= 0);
    status = H5Tset_size(t5,SIZE5);
    status = H5Tinsert(memtype1,"center",HOFFSET(MISC,center),t5);
    assert (status >= 0);
    status = H5Tset_size(t6,SIZE6);
    status = H5Tinsert(memtype1,"prot",HOFFSET(MISC,prot),t6);
    assert (status >= 0);
    status = H5Tset_size(t7,SIZE7);
    status = H5Tinsert(memtype1,"assay",HOFFSET(MISC,assay),t7);
    assert (status >= 0);
    status = H5Tset_size(t8,SIZE8);
    status = H5Tinsert(memtype1,"panel",HOFFSET(MISC,panel),t8);
    assert (status >= 0);
    status = H5Tset_size(t9,SIZE9);
    status = H5Tinsert(memtype1,"cq",HOFFSET(MISC,cq),t9);
    assert (status >= 0);
    /*create dataset*/
    dset1 = H5Dcreate (file, miscvar.c_str(), memtype1, space1, H5P_DEFAULT, cparms1, H5P_DEFAULT);  
    /*create memory for slab writes*/  
    memspace1 = H5Dget_space(dset1);

    /*Create dataset access property list for SNP dataset*/
    dataprop = H5Pcreate(H5P_DATASET_ACCESS);
    status = H5Pset_chunk_cache(dataprop,H5D_CHUNK_CACHE_NSLOTS_DEFAULT,
             SNP_CHUNK_CACHE,H5D_CHUNK_CACHE_W0_DEFAULT); /*set snp chunk size*/
    dset2 = H5Dcreate(file,snpvar.c_str(),H5T_NATIVE_INT,space2,H5P_DEFAULT,cparms2,dataprop);
    /*create memory for slab writes*/  
    memspace2 = H5Dget_space(dset2);
    
    MISC wdata[CHUNKSIZE1]; /*Write buffer*/
    int snpdata[CHUNKSIZE2_1][samcount];
    printf("sizeof 1 misc data row: %lu ; total: %lu\n",sizeof(wdata)/CHUNKSIZE1, sizeof(wdata)); 
    printf("sizeof 1 genotype data row: %lu; total: %lu\n",sizeof(snpdata)/CHUNKSIZE2_1, sizeof(snpdata)); 
    printf("linestream size:%lu\n",linestream.capacity());

    /*subetting data thru slab to avoid overloading memory */
    for(int i=0,counter1=0,counter2=0;fp!=NULL;i++){ 
	getline(fp,linestream);
	if(fp!=NULL){
            tok=NULL;
	    tok=strtok(const_cast<char*>(linestream.c_str()),"\t\n");
            /*Load miscellaneous data*/
            for (uint64_t j=0; j<FIXCOL; j++) {  
	    	if (tok==NULL){
		    cout << "ERROR: Incomplete entries on miscellaneous data of sample " << i << endl;
		    cout << "REPORT: Failed to complete writing the data" << endl;
			return 1;
	    	}else{
           	    switch(j){
		     	case 0: strcpy(wdata[counter1].rs,tok); break;
		    	case 1: strcpy(wdata[counter1].snpallele,tok); break; 
		    	case 2: strcpy(wdata[counter1].chrom,tok); break;
  		    	case 3: wdata[counter1].pos=atoi(tok); break;
                    	case 4: wdata[counter1].strand=tok[0]; break;
		    	case 5: strcpy(wdata[counter1].build,tok); break;
		    	case 6: strcpy(wdata[counter1].center,tok); break;
		    	case 7: strcpy(wdata[counter1].prot,tok); break;
		    	case 8: strcpy(wdata[counter1].assay,tok); break;
		    	case 9: strcpy(wdata[counter1].panel,tok); break;
		    	case 10: strcpy(wdata[counter1].cq,tok); break;
		    	default: break;
	   	    } /*end switch*/
           	    tok=strtok(NULL,"\t\n");
            	}
	    }
            /*Load genotype data*/
            for (uint64_t j=0; j<samcount; j++) {  
	        /*write snps*/
            	if(tok==NULL){
		    cout << "ERROR: Incomplete entries on genotype data of sample " << i << endl;
		    cout << "REPORT: Failed to complete writing the data" << endl;
		    return 1;
	    	}else{
	    	    snpdata[counter2][j]=Calls[tok];
	    	    //cout << snpdata[counter2][j] << " ";
            	    tok=strtok(NULL,"\t\n");
	    	}
	    }
	    counter1++; counter2++;
	}
        /*WRITING MISC DATA*/
        if(counter1==CHUNKSIZE1 || fp==NULL){ /*saving of slab*/
            if(i+1==CHUNKSIZE1){ /*first slab*/
		 status = H5Dwrite(dset1, memtype1, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
            }else{
		offset1[0]=newsize1[0];
		/*extend data later for the 2nd to the last slab*/
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
		status = H5Dwrite(dset1,memtype1,memspace1,space1,H5P_DEFAULT,wdata);
                status = H5Sclose(space1);
            }
 	    counter1=0; /*increment by 1 before next iteration*/
        }
	/*WRITING GENOTYPE DATA*/
        if(counter2==CHUNKSIZE2_1 || fp==NULL){
	    if(i+1==CHUNKSIZE2_1){ /*first slab*/
		status = H5Dwrite(dset2,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,snpdata);
	    }else{
		offset2[0]=newsize2[0];
		/*extend data later for the 2nd to the last slab*/
                if(fp==NULL && counter2>0){
 		     count2[0] = counter2;
		     newsize2[0]=newsize2[0]+counter2;
		     dim2[0]=counter2;
    		     memspace2 = H5Screate_simple(2,dim2,maxdim2); 
		}else{
		     newsize2[0]=newsize2[0]+counter2;
		}
		status = H5Dset_extent(dset2,newsize2);
	 	space2 = H5Dget_space(dset2);
                status = H5Sselect_hyperslab(space2, H5S_SELECT_SET,
 			(const hsize_t*)offset2,NULL, count2, NULL);
	    	status = H5Dwrite(dset2,H5T_NATIVE_INT,memspace2,space2,H5P_DEFAULT,snpdata);
	    	status = H5Sclose(space2);
	    }
	    counter2=0;
	}
 	//cout << "\n";
    }/*end subset loading*//*end subset loading*/
 
    
    status = H5Dclose(dset1);
    status = H5Dclose(dset2);
    status = H5Sclose(memspace1);
    status = H5Sclose(memspace2);
    status = H5Tclose(memtype1);
    status = H5Fclose(file);
    status= H5Pclose (cparms1);
    status= H5Pclose (cparms2);
    status= H5Pclose (dataprop);

    return 0;
}

