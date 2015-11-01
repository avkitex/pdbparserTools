#include <iostream>
#include <fstream>

#include "../commonFuncs.h"
#include "logs.h"


using namespace std;

void ParserToolsLogger::setLogDepth(int depth){
	logDepth = depth;
}

void ParserToolsLogger::setProgrammName(string programm){
	programmModule = programm;
	setOutputFile(programmModule + ".log");
}

void ParserToolsLogger::setOutputFile(string fileName){
	file = fileName;
	fileHandle.open(file.c_str(), ofstream::out | ofstream::app);
	fileHandle << "PROGRAMM: " << programmModule << "\t";
	fileHandle << "DATE: " << getDate() << "\n";
 	fileHandle.close();
}

void ParserToolsLogger::logMsg(int type, string where, string msg){
	if (logDepth >= type)
	{
		fileHandle.open(file.c_str(), ofstream::out | ofstream::app);
		if (type == ERROR_MSG){
			fileHandle << "E|";
		}
		else if (type == INFO_MSG){
			fileHandle << "I|";
		}
		else{
			fileHandle << "L|";
		}

		fileHandle << getTime()<< "\t";
		fileHandle << "In " << where << ": " << msg << "\n";

		fileHandle.close();
	}
}
ParserToolsLogger * ParserToolsLogger::inst = 0;
