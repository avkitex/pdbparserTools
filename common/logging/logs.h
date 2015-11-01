#ifndef COMMON_FUNCS
#define COMMON_FUNCS

#include <iostream>
#include <vector>
#include <ctime>

#define DEFAULTLOGDEPTH 2


/*
0 - no log
1 - errors only
2 - add normal messages
3 - add verbose messages
4 - add all all all messages as many as you want
*/

using namespace std;



class ParserToolsLogger{
private:
	enum LOG_DEPTH{
		NO_LOG = 0,
		ERRORS_ONLY = 1,
		NORMAL_LOG = 2,
		VERBOSE_LOG = 3,
		VERY_VERBOSE_LOG = 4
	}
	enum LOG_MSG_TYPES{
		ERROR = 1,
		NORMAL_LOG = 2,
		VERBOSE_LOG = 3,
		INFO_MSG = 4
	}
	int logDepth;
	string file, programmModule;
	ofstream fileHandle;
public:
	logger(string programm){
		this->logDepth = DEFAULTLOGDEPTH;
		programmModule = programm;
		setOutputFile(programmModule + ".log");
	}

	void setLogDepth(int depth);
	void setOutputFile(string fileName);

	void logMsg(int type, string where, string msg);

};

#endif
