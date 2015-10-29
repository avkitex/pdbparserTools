#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstdlib>


#include "../common/residue/pdbqtFileStringsRW.h"
#include "../common/residue/pdbqtFileStringsRW.cpp"

#include "../common/residue/commonFuncs.h"
#include "../common/residue/commonFuncs.cpp"

#include "elementsRepo/elementsRepo.h"
#include "elementsRepo/elementsRepo.cpp"

#define VERSION "0.1"


using namespace std;

void help(){
    cout << "Mass counter v " << VERSION << "\n";
    cout << "Usage:\n";
    cout << "massCounter multipdbqtfile\n";
    return;
}



int main(int argc, char ** argv)
{
    string inputFile;
    if (argc != 2){
        cerr << "There must be pdbqt file specified\n";
        help();
        exit(0);
    }
    else{
        inputFile = argv[1];
    }
    elementsRepo elements();
    pdbqtFileStrings file;
    multipdbqtFileStringsReader multipdbqtfile(inputFile, 0);



}

