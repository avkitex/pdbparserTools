#ifndef LOGGER
#define LOGGER

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#define DEFAULTLOGDEPTH 2

#define log(type, func, msg) ParserToolsLogger::I()->logMsg(type, func, msg)


/*
0 - no log
1 - errors only
2 - add normal messages
3 - add verbose messages
4 - add all all all messages as many as you want
*/

using namespace std;


enum LOG_DEPTH
{
	NO_LOG = 0,
	ERRORS_ONLY = 1,
	NORMAL_LOG = 2,
	VERBOSE_LOG = 3,
	VERY_VERBOSE_LOG = 4
};

enum LOG_MSG_TYPES
{
	ERROR_MSG = 1,
	NORMAL_LOG_MSG = 2,
	VERBOSE_LOG_MSG = 3,
	INFO_MSG = 4
};

class ParserToolsLogger
{
public:
	static ParserToolsLogger * I() {
        if(!inst){
            inst = new ParserToolsLogger();
        }
        return inst;
    }
	void setProgrammName(string programm);

	void setLogDepth(int depth);

	void setOutputFile(string fileName);

	void logMsg(int type, string where, string msg);
private:
	int logDepth;
	string file;
	string programmModule;
	ofstream fileHandle;

	ParserToolsLogger(){
		logDepth = DEFAULTLOGDEPTH;
	}
    static ParserToolsLogger * inst;
};

#endif
