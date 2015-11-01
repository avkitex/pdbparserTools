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
//ATOM      1  C   UNK A  0       58.634  43.041  74.974  0.00  0.00    +0.600 C
//ATOM     22 HD21 ASN A  1       61.505  44.535  78.193  0.00  0.00    +0.395 HD
	string el = trim(s.substr(77, 2)), atom = trim(s.substr(12, 2)), res = "";
	if (el.size() == 2 && isUpperLetter(el[1])){
		res += atom[0];
	}
	else{
		res += atom;
	}
	return res;
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
    mout << "Ligand\tMass\n";
    while (file.size()){
		mout << file.name << "\t" << countMass(file) << "\n";
		multipdbqtfile.getNextPdbqt(file);
    }
    multipdbqtfile.close();
	mout.close();
	return 0;
}

