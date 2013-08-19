/*
 *Last Update: Aug. 16, 2013
 *Author: Roven Rommel B. Fuentes
 *TT-Chang Genetic Resources Center, International Rice Research Institute
 *
 *Load vcf data to HDF5.
 *This version only allows SNP calls (no indels/structural variants)
*/

#include "parsevcf.h"

static const char *options="f:F:p:P:n:N:i:I:c:C:r:R:";

static string datafile;
static string inputfile;
static string varName;
static string varPath;

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
    vector<char> callvec;
    string linestream,inipath,gpath1,gpath2;
    hid_t file,group,group_sub1,group_sub2;
    META_2* contig=NULL;
    MISC* misc = NULL;
    META_3* info=NULL, *format = NULL;
    int **call=NULL;
    int contigcount=0,infocount=0,formcount=0,samcount=0, indelcount = 0;;
    
    /*create new file*/
    inipath=varPath+varName;
    file = H5Fcreate(datafile.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT); 
    group = H5Gcreate (file, inipath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    gpath1 = inipath + "/meta"; 
    group_sub1 = H5Gcreate (file, gpath1.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    gpath2 = inipath + "/genfields"; 
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
 	    	parseHeader1(linestream,headmap); 
	    }
	}else if(linestream.substr(0,6) == "#CHROM"){ 
            samcount = loadSampleNames(linestream,file,inipath);
	    break; 
	}else{
            cout << "ERROR: Failed to write the input file. Unknown format of header entry. \"" << "\"" << endl;
	    cout << "REPORT: Failed to complete writing the data" << endl;
	    return 1; /*finished reading header*/
	}
    }
    loadHeader1(headmap,file,gpath1); 
    loadHeader2(contig,file,contigcount,contigmap,gpath1);
    loadHeaderInfoFormat(false,info,infocount,infomap,file,gpath1);
    loadHeaderInfoFormat(true,format,formcount,formmap,file,gpath1);
    
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

    for(i=0;fp!=NULL && i<10;i++){ 
        getline(fp,linestream);
	if(fp!=NULL){
            lastidx = parseMisc(linestream,misc,infomap,formmap,contigmap,callvec,counter1,indelcount);
            counter1++;
            parseGenotypes(linestream,call,counter2,lastidx+1,callvec,states,format,formmap,misc[counter1-1].format);
            counter2++;
        }
        if(counter1==CHUNKSIZE1 || fp==NULL){
            if(i+1==CHUNKSIZE1){ //first slab
                //set compound datatype and spaces
  		setH5MiscFormat(file,memtype1,space1,dset1,cparms1,dataprop1,counter1,inipath);
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
  		setH5CallFormat(file,space2,dset2,cparms2,dataprop2,counter2,samcount,inipath);
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
    status = H5Dclose(dset1);
    status = H5Tclose(memtype1);
    status = H5Sclose(memspace1);
    status = H5Pclose(cparms1);
    status= H5Pclose (dataprop1);
    status = H5Dclose(dset2);
    status = H5Sclose(memspace2);
    status = H5Pclose(cparms2);
    status= H5Pclose (dataprop2); 
    status = H5Gclose(group);
    status = H5Gclose(group_sub1);
    status = H5Gclose(group_sub2);
    free(info);
    free(format);
    free(call[0]);
    free(call);	
    free(misc);
    H5Fclose(file);
    
}
