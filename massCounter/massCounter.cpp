#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstdlib>


#include "../common/pdbqtFileStringsRW.h"

#include "../common/commonFuncs.h"

#include "elementsRepo/elementsRepo.h"

#define VERSION "1.0"


using namespace std;

void help(){
    cout << "Mass counter v " << VERSION << "\n";
    cout << "Usage:\n";
    cout << "massCounter multipdbqtfile\n";
    return;
}

string getElementPdbqtStrings(string s){
	return trim(s.substr(12, 2));
}

double countMass(pdbqtFileStrings & ligand){
	double mass = 0;
	string element, name;
	for (int i = 0; i < ligand.size(); ++i){
		if (ligand.strings[i].substr(0, 4) == "ATOM"){
			element = getElementPdbqtStrings(ligand.strings[i]);
			mass += elementsRepo::Instanse()->getValue(element);
		}
	}
	return mass;
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
    ofstream mout("masses.txt");
    elementsRepo elements();
    pdbqtFileStrings file;
    multipdbqtFileStringsReader multipdbqtfile(inputFile, 0);
    multipdbqtfile.getNextPdbqt(file);
    while (file.size()){
		mout << file.name << "\t" << countMass(file) << "\n";
		multipdbqtfile.getNextPdbqt(file);
    }
	mout.close();
	return 0;
}

