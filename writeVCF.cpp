/*
 *Last Update: Aug. 23, 2013
 *Author: Roven Rommel B. Fuentes
 *TT-Chang Genetic Resources Center, International Rice Research Institute
 *
 *
 *
 *Descrption: Load vcf data to HDF5.
 *This version only allows diploid SNP calls (no indels/structural variants)
*/

//Update: loaded FORMAT fields; include string and field w/ Number='.'

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
    parseArgs(argc, argv);
    if (datafile.empty() || inputfile.empty() || varName.empty()
	|| varPath.empty()) {
        cerr << "Usage:\n" << *argv
                  << " -f data-file-name\n"
		  << " -i input-data-file-name\n"
		  << " -n variable-name\n"
                  << " -p variable-path\n"
                  << endl;
        return 1;
    }
    
    return writeVCF(inputfile,varPath,varName,datafile);
}
