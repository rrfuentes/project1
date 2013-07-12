/**
   A simple program to test the capability of building and storing
   indexes to a file.

   binning option is specified in the form of "<binning ... />".  
   NOTE that on most systems, the binning option needs to be quoted 
   because it involved characters that have special meaning to most shells.
*/

#include "indexBuilder.h"
#include <iostream>	
#define SET 14

static std::string datafile;
static std::string indexfile;
static char binning[] = "precision=2";
static bool forcerebuild = false;
static std::string varPath;
static std::string varName;
static std::string logfile;
static int mpi_len = 0; //1600; //1000000;
static int mpi_dim = 0;

void parseArgs(int argc, char **argv) {
    static const char *options="f:F:p:P:n:N:i:I:v:V:l:L:g:G:rR";
    extern char *optarg;
    int c;
    while ((c = getopt(argc, argv, options)) != -1) {
        switch (c) {
	case 'f':
	case 'F': datafile = optarg; break;
	case 'i':
	case 'I': indexfile = optarg; break;
	case 'g':
	case 'G': logfile = optarg; break;
	case 'p':
	case 'P': varPath = optarg; break;
	case 'n':
	case 'N': varName = optarg; break;
	case 'r':
	case 'R': forcerebuild = true; break;
	case 'l':
	case 'L': mpi_len = atoi(optarg);break; 
	default: break;
        } // switch
    } // while
} // parseArgs


int main(int argc, char** argv) {
    ibis::horometer timer;
    bool berr = true;
    int ret = 0;
    timer.start();
    
    parseArgs(argc, argv);
    if (datafile.empty() || indexfile.empty() || varName.empty() || varPath.empty()) {
        std::cerr << "Usage:\n" << *argv 
                  << " -f data-file-name[hdf5]" 
		  << "-i index-file-name"
		   "[-g log-file]" 
	    	   "[-n variable-name]"
		   "[-r forcerebuild]"
		   "[-p variable-path]"
	    	   "[-b '<binning nbins=1000 />' (default unbinned)]"
	    	   "[-l mpi_subarray_size(default=100000)]\n"
	      << std::endl;
        return -1;
    }

#ifndef FQ_NOMPI 
    MPI_Init(&argc, &argv);
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

    if (! logfile.empty()){
#ifndef FQ_NOMPI
	std::ostringstream oss;
	oss << logfile << "-" << mpi_rank << ".log";
	logfile = oss.str();
#endif
        std::cout << *argv << " is to redirect log messages to \""
		  << logfile << "\" ..." << std::endl;
	ibis::util::setLogFileName(logfile.c_str());
    }
    //strcpy(binning,"precision=2");
    ibis::gParameters().add(FQ_REPORT_STATISTIC, "false");
    ibis::gParameters().add("fileManager.maxBytes", "2GB");
    if (forcerebuild) {
	ibis::gParameters().add("FastQuery.forceIndexRebuild", "true");
    } else {
	ibis::gParameters().add("FastQuery.forceIndexRebuild", "false");
    }
    
    FQ::FileFormat model = FQ::FQ_HDF5;

    if (!logfile.empty()) {
	ibis::util::logger lg;
	lg() << *argv << " data file \"" << datafile << "\"";
	if (! indexfile.empty())
	    lg() << "\tindexfile \"" << indexfile.c_str() << "\"";
	if (binning !=0)
	    lg() << "\tbinning option \"" << binning << "\"";
	if (! varPath.empty())
	    lg() << "\tvariable path \"" << varPath << "\"";
	if (! varName.empty())
	    lg() << "\tvariable name \"" << varName;
    }

    ibis::util::timer totTimer(*argv, 1);
    
    IndexBuilder* indexBuilder = new IndexBuilder(datafile, model, indexfile, 0,"",""); 
    if (!indexBuilder->isValid()) {
	LOGGER(ibis::gVerbose >= 0)
	    << "ERROR: Failed to initialize the IndexBuilder object for file \"" 
	    << datafile << "\"";
	delete(indexBuilder);
#ifndef FQ_NOMPI
        MPI_Finalize();
#endif
	return -2;
    }

    std::string variable;
    std::vector<uint64_t> dims;
    FQ::DataType type;
    char buffer[150];
    std::string vartemp;
    int beg=0, end=0;

    if (! indexBuilder->getVariableInfo(varName, variable, dims, &type, varPath)) {
	    	printf("ERROR: Failed to get the information for variable\n");
	    	berr = false;
                delete(indexBuilder);
    }else {
 	if(dims.size()!=2){
	    printf("ERROR: The data has an invalid dimension. SNP data should be in 2D matrix only.\n");
	    delete(indexBuilder);
#ifndef FQ_NOMPI
    MPI_Finalize();
#endif
	    return -1;
	}

        delete(indexBuilder);
	//building bitmap indices for rows
	if(mpi_len>0){ 
#ifndef FQ_NOMPI
        mpi_len=mpi_len/mpi_size;
#endif
	     
	     for (int i =0; i<dims[1]; i++) {
		IndexBuilder* indexBuilder = new IndexBuilder(datafile, model, indexfile, 0,"",""); 
		sprintf(buffer,"%s[,%d]",varName.c_str(),i);
		vartemp = buffer;
		ret += indexBuilder->buildIndexes(binning, varPath, vartemp.c_str(), mpi_dim, mpi_len);
		delete(indexBuilder);
	    }
	}else{
	    for (int i =0; i<dims[1]; i++) { //index for all rows/SNPs
		IndexBuilder* indexBuilder = new IndexBuilder(datafile, model, indexfile, 0,"",""); 
		sprintf(buffer,"%s[,%d]",varName.c_str(),i);
		vartemp = buffer;
		//std::cout << vartemp << " ";
		ret += indexBuilder->buildIndexes(binning, varPath, vartemp.c_str()); 
		delete(indexBuilder);
	    }
	}
    }
    timer.stop();
    printf("Bitmap Indexing Time:%f\n", timer.realTime());

#ifndef FQ_NOMPI
    MPI_Finalize();
#endif
    if (berr) {
    		printf("REPORT: Successfully completed indexing data.\n");
    		return 0;	
    } else {
    		printf("REPORT: Failed to complete indexing data.\n");
    		return -1;	
    }
    return ret;
} // main
	
