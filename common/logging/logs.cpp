#include "logs.h"



void ParserToolsLogger::setLogDepth(int depth){
	logDepth = depth;
}
void ParserToolsLogger::setOutputFile(string fileName){
	file = fileName;
	fileHandle.open(fileName.c_str());
	fileHandle << "PROGRAMM MODULE " << programmModule << "\n";
	fileHandle << "DATE: "<< ltm->tm_mday << ":" << 1 + ltm->tm_mon << ":" << 1900 + ltm->tm_year << "\n";
 	fileHandle.close();
}

void ParserToolsLogger::logMsg(int type, string where, string msg){
	if (logDepth >= type)
	{
		fileHandle.open(fileName, ofstream::out | ofstream::app);
		if (type == ERROR){
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

