#include <cstring>
#include <fstream>
#include <vector>

#include "commonFuncs.h"
#include "commonFuncs.cpp"
#include "vinaLogFileParser.h"

using namespace std;

bool vinaResult::parseLogFile(string log_file)
{
    store.clear();
    string s;
    bool dataSection = false;
    vector <string> buf;
    dockingResultValues resVal;
    ifstream din(log_file.c_str());
    if (!din){
        return 1;
    }
    while (getline(din, s)){
        if (dataSection){
            buf = split(s, " ");
            if (buf.size() > 0){
                if (isdigit(buf[0][0])) {
                    if (buf.size() == 4) {
                        resVal.energy = strtodoub(buf[1]);
                        resVal.rmsdLB = strtodoub(buf[2]);
                        resVal.rmsdUB = strtodoub(buf[3]);
                        store.push_back(resVal);
                    }
                    else {
                        return 1;
                    }
                }
                else if (s.find("done") != -1){
                    finished = true;
                }
            }

        }
        else{
            if (s.substr(0, 4) == "mode"){
                dataSection = true;
            }
        }
    }
    din.close();
    return store.size() == 0;
}

