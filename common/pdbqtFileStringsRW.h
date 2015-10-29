#ifndef PDBQT_FILE_STRINGS_RW
#define PDBQT_FILE_STRINGS_RW

class pdbqtFileStrings{
public:
    vector <string> strings;

    string name;

    pdbqtFileStrings(){
		clear();
    }
    ~pdbqtFileStrings(){}

    void clear();

    bool outpdbqtFileStrings();

    void add(string s);

    int size();
};

class multipdbqtFileStringsReader{
public:
    ifstream fileHandle;
    int modelsReturned;

    multipdbqtFileStringsReader(string file, int skipStr);

    ~multipdbqtFileStringsReader(){}

    void getNextPdbqt(pdbqtFileStrings &file);

    void close();
};

class multipdbqtFileStringsConstructor{
public:
    ofstream fileHandle;
    int modelCounter;
    string fileName;

    multipdbqtFileStringsConstructor(string file);

    ~multipdbqtFileStringsConstructor(){}

    void add(pdbqtFileStrings & file);
};

#endif
