//Last update 23.10.15 17:20
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <queue>
#include <stdio.h>

#include <unistd.h>
#include <ctime>

#define VERSION "0.1"

using namespace std;

void help(string p_name)
{
    cout << "Parallel multiPdbqtLigandsFile program starter v" << VERSION << "\n";
    cout << "Usage:\n\n";
    cout << p_name << " theProgram multiPdbqtLigandsFile.pdbqt restParamsStr\n";
    cout << "programm - full path for the programm to be run\n";
    return;
}

string itostr(int n)
{
    string s = "0";
    while (n > 0)
    {
        s += char(n % 10 + '0');
        n /= 10;
    }
    return s;
}

class pdbqtFile{
    vector <string> strings;

    string name;

    pdbqtFile(){
        clear();
    }

    void clear(){
        strings.clear();
        name = "";
    }

    bool outPdbqtFile(){
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

    void addString(string s){
        strings.push_back(s);
    }

    int size(){
        return strings.size();
    }
};

class multiPdbqtFileReader{
    ifstream fileHandle;
    int modelsReturned;

    multiPdbqtFileReader(string file, int skipStr){
        fileHandle = new ifstream(file.c_str());
        modelsReturned = 0;
        while (getline(fileHandle, s) && skipStr > 0){
            if (s.size() >5 && s.substr(0, 6) == "ENDMDL"){
                modelsReturned++;
                skipStr--;
            }
        }
    }

    pdbqtFile getNextPdbqt(pdbqtFile file){
        string s;
        file.clear();
        while (getline(fileHandle, s)){
            if (s.size() >5 && s.substr(0, 6) == "ENDMDL"){
                if (file.size()){
                    modelsReturned++;
                    return file;
                }
            }
            if (s.size() > 4 && s.substr(0, 5) == "MODEL"){
                continue;
            }
            if (s.size() > 15 && s.substr(0, 15) == "REMARK  Name = "){
                file.name = s.substr(15, s.size() - 15);
            }
            file.addString(s);
        }
        return file;
    }

    void close(){
        fileHandle.close();
    }
};

class multiPdbqtFileConstructor{
    ofstream fileHandle;
    int modelCounter;

    multiPdbqtFileConstructor(string file){
        fileHandle = new ofstream(file.c_str());
        modelCounter = 1;
    }

    void add(pdbqtFile & file){
        fileHandle << "MODEL " << modelCounter << "\n";
        for (int i = 0; i < file.size(); ++i){
            fileHandle << file.strings[i] << "\n";
        }
        fileHandle << "ENDMDL\n";
    }
    void close(){
        fileHandle.close();
    }
}

class freeCores{
	queue <int> cores;

	freeCores(){
        cores.clear();
	}

	void add(int rank){
        cores.push(rank);
	}

	int get(){
        if (cores.size()){
            int a = cores.front();
            cores.pop();
            return a;
        }
        return 0;
	}
};

class logSummaryConstructor{
    ofstream fileHandle;
    logSummaryConstructor(string file){
        fileHandle = new ofstream(file.c_str());
    }

    void add(string name, vinaOutLog &vinaRes){
        fileHandle << name;
        for (int i = 0; i < vinaRes.store.size(); ++i){
            fileHandle << "\t" << vinaRes.store[i];
        }
        fileHandle << "\n";
    }
    void finalize(){
        fileHandle.close();
    }
};

void deleteFile(string name){
    return remove(name.c_str());
}


int main(int argc, char **argv)
{
	int rank, a, am_works = 0, cores;
	//time_t timer;
	string multiPdbqtFileName, programName, restParams;

	vector <string> arguments;
	vector <string> path;
	if (argc == 4)
	{
		programm = argv[1];
        multiPdbqtFileName = argv[2];
        restParams = argv[3];
	}
	else
	{
		cout << "Wrong amount of arguments\n";
		help(argv[0]);
		exit(0);
	}


    char receivedCharArr[100];
    string receivedStr;


	MPI_Init(&argc, &argv);
	MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &cores);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        multiPdbqtFileConstructor MAINOUTPUT_pdbqtFileConstructor("AllResults.pdbqt");
        logSummaryConstructor MAINOUTPUT_logFileConstructor("AllLogs.logs");
        multiPdbqtFileReader masterFile(multiPdbqtFileName, 0);
        multiPdbqtFileReader vinaOutPdbqt;
        vinaResult vinaOutLog;
        pdbqtFile file, vinaResFile;
        cout << rank << " I'm the BOSS\n";
        masterFile.getNextPdbqt(file);
        while (file.size())
        {
            if (file.outPdbqtFile()){
                cout << "Error with file " << file.name << "\n";
            }
            MPI_Recv(receivedCharArr, 100, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (strlen(receivedCharArr)){
                receivedStr = new string(receivedCharArr);

                vinaOutPdbqt = new multiPdbqtFileReader(receivedStr + ".pdbqt", 0);
                vinaOutPdbqt.getNextPdbqt(vinaResFile);
                vinaOutPdbqt.close();
                MAINOUTPUT_pdbqtFileConstructor.add(vinaResFile);
                deleteFile(receivedStr + ".pdbqt");

                vinaOutLog.parseLogFile(receivedStr + ".log");
                MAINOUTPUT_logFileConstructor.add(vinaOutLog)
                deleteFile(receivedStr + ".log");
            }


            cout << "New job " << file.name << " for process " << status.MPI_SOURCE << "\n";

            MPI_Send(file.name.c_str(), file.name.size(), MPI_CHAR, status.MPI_SOURCE, 11, MPI_COMM_WORLD);




            masterFile.getNextPdbqt(file);
        }
        int countCores = 0;
        cout << "No more jobs!\n";
        while (countCores < cores - 1)
        {
            MPI_Recv(receivedCharArr, 100, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            countCores++;
            receivedCharArr = "";
            MPI_Send(receivedCharArr, 0, MPI_CHAR, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD);
        }
    }
    else
    {
        //time(&timer);
        receivedCharArr = "";
        MPI_Rsend(receivedCharArr, 0, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
        MPI_Recv(receivedCharArr, 100, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        receivedStr = new string(receivedCharArr);
        while (receivedStr.size())
        {
            cout << rank << " received ligand " << receivedStr <<  "\n";

            s = programm + " " + restParams + " --ligand " + receivedStr + ".pdbqt --log " + receivedStr + ".log";

            system(s.c_str());

            MPI_Send(receivedStr.c_str(), receivedStr.size(), MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
            MPI_Recv(receivedCharArr, 100, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            receivedStr = new string(receivedCharArr);
        }
    }
	MPI_Finalize();
	return 0;
}

