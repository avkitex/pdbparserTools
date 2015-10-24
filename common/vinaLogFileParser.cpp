#include "../common/commonFuncs.cpp"


class dockingResultValues{
    double energy, rmsdLB, rmsdUB;

    dockingResultValues(){
        energy = 0;
        rmsdLB = 0;
        rmsdUB = 0;
    }
};

class vinaResult{
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
    bool parseLogFile(string log_file)
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
                            store.add(resVal);
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
                if (buf.substr(0, 4) == "mode"){
                    dataSection = true;
                }
            }
        }
        din.close();
        return store.size() == 0;
    }
}
