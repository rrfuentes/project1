/*
 *version: 2.0
 * Load the specified row/s in memory and compare(serial/pthread) against a reference row.
 * Last Updated:July 17, 2013
 * Works on data > memorysize 
 *
 *Author: Roven Rommel B. Fuentes
 *TT-Change Genetic Resources Center, International Rice Research Institute
 *
 */

#include "queryProcessor.h"
#include "hdf5file.h"
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <pthread.h>
#include <time.h>
#include <sstream>
#include <map>
#include <sys/resource.h>

#define NUMTHREADS 2
using namespace std;

static const char *options="f:F:p:P:n:N:r:R:o:O:c:C";
static string pos;
static string snpbound;
static string datafile;
static string varPath;
static string varName;
static string outfile;
int mpi_size, mpi_rank;



struct threadData {
    int* data;
    int* refdata;
    int* consensus;
    int snpcount;
    int blocksize;
    int tid;
    int tcount;
};

void parseArgs(int argc, char **argv){
    extern char *optarg;
    int i;
    while((i=getopt(argc, argv, options)) != -1){
	switch(i){
	    case 'f':
            case 'F': datafile = optarg; break;
	    case 'p':
	    case 'P': varPath = optarg; break;
	    case 'n':
	    case 'N': varName = optarg; break;
	    case 'r':
	    case 'R': pos = optarg; break; //e.g: 'x|y:z' where x=ref sample, y=sample for comparison, y:z=range
	    case 'c':
	    case 'C': snpbound = optarg; break; //e.g: x:y, where x=start, y=end
	    case 'o':
	    case 'O': outfile = optarg; break; 
	    default: break;
	}
    }


}

void *compareSample(void* thread){
    threadData *t = (threadData*) thread;
    int snpcount = t->snpcount;
    int blocksize = t->blocksize;
    int tid = t->tid;
    int start=0, end=0,idxrel=0,temp=0;
    /*VERTICAL CHOPPING of data*/
    int skip=snpcount/t->tcount; /*rows/snps per chunk*/
    if(snpcount%t->tcount>0 && tid==t->tcount-1) /*Last chunk is larger*/
		skip += (snpcount%t->tcount);

    for(int x=0;x<skip;x++){
	start=(tid*skip*blocksize)+(blocksize*x); 
	end=start+blocksize; 
	idxrel=(tid*skip)+x; /*snp index*/			
	for(int i=start; i<end;i++){ /*compare samples*/
	    t->consensus[idxrel]=t->consensus[idxrel]||(t->data[i]!=t->refdata[idxrel]);
	    t->data[i]=(t->data[i]!=t->refdata[idxrel])?t->data[i]:25;
	}
    }
    return 0;
}

int main(int argc, char **argv) {
    FILE *output;
    ibis::horometer timer1,timer2;
    timer1.start(); 
    timer2.start();
    int CHOP=1;
    
    const rlim_t STACK_SIZE = 1000*1024*1024; 
    struct rlimit rl;
    rl.rlim_cur = STACK_SIZE;
    int ret = setrlimit(RLIMIT_STACK,&rl);

    map<int, string> Calls;
    Calls[0]="AA";
    Calls[1]="AT";
    Calls[2]="AC";
    Calls[3]="AG";
    Calls[4]="AN";
    Calls[5]="TA";
    Calls[6]="TT";
    Calls[7]="TC";
    Calls[8]="TG";
    Calls[9]="TN";
    Calls[10]="CA";
    Calls[11]="CT";
    Calls[12]="CC";
    Calls[13]="CG";
    Calls[14]="CN";
    Calls[15]="GA";
    Calls[16]="GT";
    Calls[17]="GC";
    Calls[18]="GG";
    Calls[19]="GN";
    Calls[20]="NA";
    Calls[21]="NT";
    Calls[22]="NC";
    Calls[23]="NG";
    Calls[24]="NN";
    Calls[25]="--";

    parseArgs(argc, argv);
    if(datafile.empty() || pos.empty() || varName.empty() || varPath.empty() || outfile.empty()){
		std::cerr << "Usage:\n" << *argv
                  << " -f data-file-name"
		  << " -n variable-name"
                  << " -p variable-path"
		  << " -r ref&row indices ('x|y:z' where x=ref sample,x:y as row range) "
		  << " -r snp bounds (x:y, where x=start, y=end)"
		  << " -d variable-dimension (e.g. 2:2)"
		  << " -o output-file"
                  << std::endl;
    }

#ifndef FQ_NOMPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

    FQ::FileFormat model = FQ::FQ_HDF5;
    bool berr = true;
    QueryProcessor* queryProcessor = new QueryProcessor(datafile, model, "", 0, "",""); 
	
    if (queryProcessor->isValid() == false) {
	printf("ERROR: Failed to initiate query processor for file.\n");
 	berr = false;
    }

    string variable;
    vector<uint64_t> dims;
    FQ::DataType type;
    if (! queryProcessor->getVariableInfo(varName, variable, dims, &type, varPath)) {
 	printf("ERROR: Failed to get the information for variable\n");
	berr = false;
    } else {
	if(dims.size()!=2){ /*dims is derived from the data*/
	    printf("ERROR: The data has an invalid dimension. SNP data should be in 2D matrix only.\n");
	    berr=false;
	}
	    	
	string str,param; 
	vector<uint64_t> sample;
	vector<uint64_t> row; //snps
	int prePos = 0, idx = 0,blocksize=1;
	FastQuery* fq = new FastQuery(datafile, model, "", 0, "",""); 
	int *refdata=NULL, *data=NULL, *consensus=NULL;
	ostringstream paramtemp,ref;

	/*Get the index for reference sample and subrows*/
	idx = pos.find('|',prePos);
	if(idx!=pos.npos && idx!=pos.length()-1){
	    str = pos.substr(prePos,idx - prePos); 
	    sample.push_back(atoi(str.c_str())); /*get the ref. sample index*/
	    prePos=idx+1;
	    idx=pos.find(':',prePos);
	    if(idx!=pos.npos){ /*indicates multiple comparison*/
		if(idx==pos.length()-1){ /*string ends with ':'*/
		    printf("ERROR: Incomplete indices specified for sample comparison.\n");
		    return 0;
		}
		str = pos.substr(prePos,idx - prePos); 
		sample.push_back(atoi(str.c_str())); /*get the start index*/
		prePos=idx+1;
		str = pos.substr(prePos,pos.length() - prePos); 
		sample.push_back(atoi(str.c_str())); /*get the end index*/
		if(sample[2]<0 || sample[2]>=dims[1] || sample[1]>=sample[2]){
		    printf("ERROR: Indices out of bounds/invalid range.\n");
		    return 0;
		}
	    }else{ 
		str = pos.substr(prePos,idx - prePos); 
		sample.push_back(atoi(str.c_str())); /*get the index of another sample*/	
	    }

	    if(sample[0]<0 || sample[0]>=dims[1] || sample[1]<0 || sample[1]>=dims[1]){
		printf("ERROR: Indices out of bounds.\n");
		return 0;
	    }
	}else{
	    printf("ERROR: Invalid indices specified for sample comparison.\n");
	    return 0;
	}

	/*Get the SNP bounds for subset sample*/
	if(!snpbound.empty()){
	    idx=snpbound.find(':',0);
	    if(idx==snpbound.npos || idx==snpbound.length()-1){
		printf("ERROR: Invalid SNP bounds.");
		return 1;
	    }
	    str = snpbound.substr(0,idx);
	    row.push_back(atoi(str.c_str()));
	    str= snpbound.substr(idx+1,snpbound.length()-idx+1);
	    row.push_back(atoi(str.c_str()));
	    if(row[0]<0 || row[1]>=dims[0] || row[1]<0 || row[0]>=row[1]){
	    	printf("ERROR: Invalid SNP bounds.\n");
	        return 0;
	    }
	    dims[0]=row[1]-row[0]+1; /*dims is now the SNP bounds for subregion*/
	}

	if(sample.size()==3){
	    blocksize=sample[2]-sample[1]+1; /*block of data*/
	}
	
        if((dims[0]+(dims[0]*blocksize)+dims[0])*sizeof(int)>(1000*1024*1024)){ 
	    //printf("Error: Insufficient memory to handle huge block.\nREPORT: Failed to complete comparing data.\n");
	    //return 1;
	    printf("\nData is greater than the available/alloted memory space.\n");
            CHOP=10;
     	}

        printf("Running with %d thread/s.\n",NUMTHREADS);
        int rowchunk=dims[0]/CHOP, offset; /*CHOP is 1 if data fits in memory*/
	float com_time=0;
	ostringstream outtext;
	pthread_t threads[NUMTHREADS];
	threadData *thread_data = (threadData*)malloc(NUMTHREADS*sizeof(threadData));
        refdata=(int*)malloc(rowchunk*sizeof(int)); 
	data=(int*)malloc((rowchunk*blocksize)*sizeof(int));
	consensus=(int*)calloc(rowchunk,sizeof(int));
        output=fopen(outfile.c_str(),"w");
        
        /*Print output header*/
	if(sample.size()==3) 
	    fprintf(output,"Reference:%lu\nBlock samples:%lu-%lu\n",sample[0],sample[1],sample[2]);
	else 
	    fprintf(output,"Reference:%lu\nSample:%lu\n",sample[0],sample[1]);
	fprintf(output,"SNPIdx\tRef\t");
	for(int i=sample[1];i<sample[1]+blocksize;i++){ 
	    fprintf(output,"%d\t",i);
	}
	fprintf(output,"\n");

	if(!snpbound.empty()){  /*set start position if bounded*/
	    offset=(int)row[0];
	}
	int rem=0;
        if(dims[0]%CHOP!=0){
            rem = dims[0]%CHOP;
	    CHOP++; /*another chunk for the remainder*/
	}

        for(int h=0;h<CHOP;h++){
	    if(h+1==CHOP && rem!=0){
		rowchunk=rem;
		free(data);
		free(refdata);
		free(consensus);
		data=(int*)malloc((rowchunk*blocksize)*sizeof(int));
		refdata=(int*)malloc(rowchunk*sizeof(int));
		consensus=(int*)calloc(rowchunk,sizeof(int));
		printf("Remainder chunk.%d\n",rowchunk);
	    }
            
	    if(!snpbound.empty()){
		ref << variable << "[" << offset <<":"<< offset+rowchunk << "," <<  sample[0] << "]";
		if(sample.size()==3)
		    paramtemp << variable << "[" << offset <<":"<< offset+rowchunk << "," << sample[1] << ":" << sample[2]+1  <<"]"; 
		else
		    paramtemp << variable << "[" << offset <<":"<< offset+rowchunk << "," << sample[1] << "]";
	    }else{
		ref << variable << "[:," << sample[0] << "]";
		if(sample.size()==3)
		    paramtemp << variable << "[:," << sample[1] << ":" << sample[2]+1 << "]";
		else
		    paramtemp << variable << "[:," << sample[1] << "]";
	    }
	
	    param = ref.str(); 
	    fq->getData(param,refdata); 
	    param = paramtemp.str();  
	    cout << h << "Ref:"<< ref.str() <<"\tParam:"<<paramtemp.str()<<"\n\n";
	    fq->getData(param,data); /*param=var[:,0:2]*/

            /*printf("Reference Row:\n");
	    for(int i=0; i<rowchunk;i++) printf("%d ",refdata[i]);
	    printf("\n\n");

	    printf("Comparison Row Block:\n");
	    for(int x=0;x<(blocksize*rowchunk);x++){
	    	printf("%d ",data[x]); 
	    	if((x+1)%blocksize==0) printf("\n\n");
	    }*/

            /*THREADING of the comparison*/	
	    for(int i=0;i<NUMTHREADS; i++){
	    	thread_data[i].refdata = refdata;
	    	thread_data[i].data = data;
	    	thread_data[i].consensus=consensus; 
	    	thread_data[i].snpcount = rowchunk;
	    	thread_data[i].blocksize = blocksize;
	    	thread_data[i].tid = i;
	    	thread_data[i].tcount = NUMTHREADS;
	    	pthread_create(&threads[i],NULL,compareSample, (void*) &thread_data[i]);
	    }

	    for(int i=0;i<NUMTHREADS;i++){
	    	pthread_join(threads[i],NULL);
	    }

            /*RESULT printing*/
            timer1.stop();
	    for(int i=0;i<rowchunk;i++){
	    	if(consensus[i]==1){
		    outtext << i+offset << "\t" << Calls[refdata[i]] << "\t";
		    for(int x=0;x<blocksize;x++){
		    	outtext << Calls[data[i*blocksize+x]] << "\t"; 
		    }
		    fprintf(output,"%s\n",outtext.str().c_str());
		    outtext.str("");
	    	}
			
	    }
 	    offset+=rowchunk;
	    paramtemp.str("");
	    ref.str("");
            outtext.str("");
	    timer1.resume();
        }
	timer1.stop();
	printf("Comparison Time:%f\n", timer1.realTime());
	
	
	free(data);
	free(refdata);
	free(consensus);
	free(thread_data);
	fclose(output);
    }	
	delete(queryProcessor);
#ifndef FQ_NOMPI
    MPI_Finalize();
#endif
    timer2.stop();
    if (berr) {
    	printf("REPORT: Successfully completed comparing data.\n Total time elapsed:%f\n", timer2.realTime());
    	return 0;	
    } else {
    	printf("REPORT: Failed to complete comparing data.\n");
    	return -1;	
    }
	
}


