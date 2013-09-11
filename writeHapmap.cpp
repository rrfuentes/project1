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
#define size1 15
#define size2 3
#define size3 5
#define size4 8
#define size5 20
#define size6 20
#define size7 20
#define size8 20
#define size9 3
#define size10 3

using namespace std;

static const char *options="f:F:p:P:n:N:i:I:c:C:r:R:";

static char *fileModel = 0;
static string datafile;
static string inputfile;
static string varName;
static string varPath;
static int LENGTH=0;
static int snp=0;



typedef struct{
    char rs[size1];
    char snpallele[size2];
    char chrom[size3];
    int pos;
    char strand;
    char build[size4];
    char center[size5];
    char prot[size6];
    char assay[size7];
    char panel[size8];
    char cq[size9];
}column_t;


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
	case 'R': snp = atoi(optarg); break; //number of genotype columns 
        case 'c':
	case 'C': LENGTH = atoi(optarg); break; //number of samples 
	default: break;
        } // switch
    } // while
} // parseArgs

int main(int argc, char **argv){
    const rlim_t STACK_SIZE = 1000*1024*1024; 
    struct rlimit rl;
    rl.rlim_cur = STACK_SIZE;
    int ret = setrlimit(RLIMIT_STACK,&rl);

    parseArgs(argc, argv);
    if (datafile.empty() || inputfile.empty() || varName.empty()
	|| varPath.empty()  || LENGTH==0 || snp==0) {
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
  
    char *tok=NULL;
    int total=0, temp=0;
    string miscvar, snpvar, headervar;
    string linestream;
    hid_t file, memtype1, memtype2, memspace1,memspace2, space1, space2, dset1, dset2, dataprop; /*handles*/
    hid_t t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
    hid_t cparms1,cparms2;
    herr_t status;
    hsize_t dim1[1]={CHUNKSIZE1}; /*dataspace dimension for misc data*/
    hsize_t dim2[2]={CHUNKSIZE2_1,LENGTH}; /*dataspace dimension for genotype data*/
    hsize_t maxdim1[1]={H5S_UNLIMITED};
    hsize_t maxdim2[2]={H5S_UNLIMITED,LENGTH};
    hsize_t chkdim1[1]={CHUNKSIZE1}; /*misc chunk size*/
    hsize_t chkdim2[2]={CHUNKSIZE2_1,CHUNKSIZE2_2}; /*genotype chunk size*/
    hssize_t offset1[1]={0}; /*offset for slabs*/
    hssize_t offset2[2]={0,0};
    hsize_t newsize1[1]={CHUNKSIZE1}; /*iterate offset for slabs*/
    hsize_t newsize2[2]={CHUNKSIZE2_1,LENGTH};
    hsize_t count1[1]={CHUNKSIZE1}; /*param fo hyperslab*/
    hsize_t count2[2]={CHUNKSIZE2_1,LENGTH};
    
    if(CHUNKSIZE1<=1 && CHUNKSIZE2_1<=1){
	cout << "ERROR: Invalid chunk size." << endl;
	cout << "REPORT: Failed to complete writing the data" << endl;
	return -1;
    }

    /*create new file*/
    file = H5Fcreate(datafile.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT); 

    miscvar = varName + "_misc";
    snpvar = varName + "_snp";

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
    t10 = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(t10,size10);
    /*create compound datatype for fixed columns*/
    memtype1 = H5Tcreate(H5T_COMPOUND,sizeof(column_t));
    status = H5Tset_size(t1,size1);
    status = H5Tinsert(memtype1,"rs#",HOFFSET(column_t,rs),t1);
     assert (status >= 0);
    status = H5Tset_size(t2,size2);
    status = H5Tinsert(memtype1,"snpallele",HOFFSET(column_t,snpallele),t2);
    assert (status >= 0);
    status = H5Tset_size(t3,size3);
    status = H5Tinsert(memtype1,"chrom",HOFFSET(column_t,chrom),t3);
    assert (status >= 0);
    status = H5Tinsert(memtype1,"pos",HOFFSET(column_t,pos),H5T_NATIVE_INT);
    assert (status >= 0);
    status = H5Tinsert(memtype1,"strand",HOFFSET(column_t,strand),H5T_NATIVE_CHAR);
    assert (status >= 0);
    status = H5Tset_size(t4,size4);
    status = H5Tinsert(memtype1,"build",HOFFSET(column_t,build),t4);
    assert (status >= 0);
    status = H5Tset_size(t5,size5);
    status = H5Tinsert(memtype1,"center",HOFFSET(column_t,center),t5);
    assert (status >= 0);
    status = H5Tset_size(t6,size6);
    status = H5Tinsert(memtype1,"prot",HOFFSET(column_t,prot),t6);
    assert (status >= 0);
    status = H5Tset_size(t7,size7);
    status = H5Tinsert(memtype1,"assay",HOFFSET(column_t,assay),t7);
    assert (status >= 0);
    status = H5Tset_size(t8,size8);
    status = H5Tinsert(memtype1,"panel",HOFFSET(column_t,panel),t8);
    assert (status >= 0);
    status = H5Tset_size(t9,size9);
    status = H5Tinsert(memtype1,"cq",HOFFSET(column_t,cq),t9);
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
  
    /*create another dataset for genoytpe field names*/
    temp=sizeof(column_t);
    //getline(fp,linestream);
    //tok=strtok(const_cast<char*>(linestream.c_str()),"\t\n");
    /*for (uint64_t j=0; j<snp+FIXCOL; j++){
        if(j>=FIXCOL){
            if (tok==NULL){
	    	cout << "ERROR: Specified number of SNPs didn't match the data." << endl;
	    	cout << "REPORT: Failed to complete writing the data" << endl;
	    	return -1;
	    }else{
    		status = H5Tinsert(memtype1,tok,temp,t10);
                assert (status >= 0);
		temp+=sizeof(t10);
            	tok=strtok(NULL,"\t\n");
            }
	}
    }*/
     
    column_t wdata[CHUNKSIZE1]; /*Write buffer*/
    int snpdata[CHUNKSIZE2_1][LENGTH];
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
			return -1;
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
            for (uint64_t j=0; j<LENGTH; j++) {  
	        /*write snps*/
            	if(tok==NULL){
		    cout << "ERROR: Incomplete entries on genotype data of sample " << i << endl;
		    cout << "REPORT: Failed to complete writing the data" << endl;
		    return -1;
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
