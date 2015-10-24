#ifndef VINALOGFILES
#define VINALOGFILES

#include <vector>

using namespace std;

class dockingResultValues{
public:
    double energy, rmsdLB, rmsdUB;

    dockingResultValues(){
        energy = 0;
        rmsdLB = 0;
        rmsdUB = 0;
    }
};

class vinaResult{
public:
    vector<dockingResultValues> store;
    bool finished;

    vinaResult(){
        store.clear();
        finished = false;
    }
    void add(dockingResultValues a){
        store.push_back(a);
    }
    int size(){
        return store.size();
    }
    bool parseLogFile(string log_file);
};

#endif
