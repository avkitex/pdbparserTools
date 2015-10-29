#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>


#include "../common/pdbqtFileStringsRW.h"
#include "../common/pdbqtFileStringsRW.cpp"

#include "../common/commonFuncs.h"
#include "../common/commonFuncs.cpp"

#define VERSION "0.1"

#define ELEMENTS_FILE "elements.dat"

using namespace std;

class elementsRepo{
	map<string, double> masses;

public:
	elementsRepo(){
		readDatFile();
	}
	void readDatFile(){
		ifstream din (ELEMENTS_FILE);

		string s;
		double val;
		vector splitedStr;

		while (getline(din, s)){
			splitedStr = split(s, '\t');
			if (masses.count(splitedStr[1])){
				cerr << "Double " << splitedStr[1] << " entry\n";
			} else {
				val = strtodoub(splitedStr[3]);
				masses[splitedStr[1]] = val;
			}
		}

		din.close();
	}

	double getValue(string s){
		if (masses.count(s)){
			return masses[s];
		} else{
			cerr << "No entry for " << s << "\n";
		}
		return 0;
	}
};

int main()
{
	elementsRepo elements();
}
