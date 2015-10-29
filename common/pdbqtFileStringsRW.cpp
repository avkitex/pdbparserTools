#include "pdbqtFileStringsRW.h"
#include "commonFuncs.h"

void pdbqtFileStrings::clear(){
	strings.clear();
	name = "";
}

bool pdbqtFileStrings::outpdbqtFileStrings(){
	if (!name.size()){
		return 1;
	}
	ofstream fout((name + ".pdbqt").c_str());
	if (!fout){
		return 1;
	}
	for (unsigned int i = 0; i < strings.size(); ++i){
		fout << strings[i] << "\n";
	}
	fout.close();
	return 0;
}

void pdbqtFileStrings::add(string s){
	strings.push_back(s);
}

int pdbqtFileStrings::size(){
	return strings.size();
}

multipdbqtFileStringsReader::multipdbqtFileStringsReader(string file, int skipStr){
	string s;
	fileHandle.open(file.c_str());
	modelsReturned = 0;
	while (skipStr > 0 && getline(fileHandle, s)){
		if (s.size() >5 && s.substr(0, 6) == "ENDMDL"){
			modelsReturned++;
			skipStr--;
		}
	}
	cout << "Skipped " << modelsReturned << "\n";
}

void multipdbqtFileStringsReader::getNextPdbqt(pdbqtFileStrings &file){
	string s;
	file.clear();
	while (getline(fileHandle, s)){
		if (s.size() >5 && s.substr(0, 6) == "ENDMDL"){
			if (file.size()){
				modelsReturned++;
				cout << "Returning file " << file.name << "\n";
				return;
			}
		}
		if (s.size() > 4 && s.substr(0, 5) == "MODEL"){
			continue;
		}
		if (s.size() > 15 && s.substr(0, 15) == "REMARK  Name = "){
			file.name = s.substr(15, s.size() - 15);
		}
		file.add(s);
	}
	cout << "Nothing more to read\n";
	return;
}

void multipdbqtFileStringsReader::close(){
	fileHandle.close();
}

multipdbqtFileStringsConstructor::multipdbqtFileStringsConstructor(string file){
	fileHandle.open(file.c_str());
	modelCounter = 1;
	fileName = file;
	fileHandle.close();
}

void multipdbqtFileStringsConstructor::add(pdbqtFileStrings & file){
	fileHandle.open(fileName.c_str(), ofstream::out | ofstream::app);
	fileHandle << "MODEL " << modelCounter << "\n";
	for (int i = 0; i < file.size(); ++i){
		fileHandle << file.strings[i] << "\n";
	}
	fileHandle << "ENDMDL\n";
	fileHandle.close();
}
