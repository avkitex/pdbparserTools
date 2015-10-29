//Last update 25.10.15 00:30
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

#include "../common/vinaLogFileParser.h"
#include "../common/vinaLogFileParser.cpp"

#include "../common/pdbqtFileStringsRW.h"
#include "../common/pdbqtFileStringsRW.cpp"

#define VERSION "1.0"

using namespace std;

void help(string p_name)
{
    cout << "Parallel multiPdbqtLigandsFile program starter v" << VERSION << "\n";
    cout << "Usage:\n\n";
    cout << p_name << " theProgram multiPdbqtLigandsFile.pdbqt restParamsStr\n";
    cout << "programm - full path for the programm to be run\n";
    return;
}


class freeCores{
public:
	queue <int> cores;

	freeCores(){
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
public:
    ofstream fileHandle;
    int maxResults;
    string fileName;
    logSummaryConstructor(string file){
        fileName = file;
        fileHandle.open(fileName.c_str());
        maxResults = 3;
        fileHandle << "Name\t1_Energy";
        for (int i = 1; i < maxResults; ++i){
            fileHandle << "\t" << i + 1 << "_Energy\t" << i + 1 << "_RMSDLB\t" << i + 1 << "_RMSDRB";
        }
        fileHandle << "\n";
        fileHandle.close();
    }

    void add(string name, vinaResult &vinaRes){
        fileHandle.open(fileName.c_str(), ofstream::out | ofstream::app);
        fileHandle << name;
        if (vinaRes.store.size()){
            fileHandle << "\t" << vinaRes.store[0].energy;
        }
        int mm = min(maxResults, (int)vinaRes.store.size());
        for (int i = 1; i < mm; ++i){
            fileHandle << "\t" << vinaRes.store[i].energy;
            fileHandle << "\t" << vinaRes.store[i].rmsdLB;
            fileHandle << "\t" << vinaRes.store[i].rmsdUB;
        }
        fileHandle << "\n";
        fileHandle.close();
    }
};

void deleteFile(string name){
    remove(name.c_str());
}

string readParams(string file){
    ifstream fin(file.c_str());

    string s;

    getline(fin, s);
    fin.close();

    return s;
}

int main(int argc, char **argv)
{
	int rank, a, am_works = 0, cores;
	//time_t timer;
	string multiPdbqtFileName, programName, restParams, programm, name, dir, ext;

	vector <string> arguments;
	vector <string> path;
	if (argc == 4)
	{
		programm = argv[1];
        multiPdbqtFileName = argv[2];
        cout << argv[3] << "\n";
        restParams = argv[3];
        parseFileName(restParams, dir, name, ext);
        if (ext == "txt" || ext == "param")
        {
            restParams = readParams(argv[3]);
        }
        cout << restParams << "\n";

	}
	else
	{
		cout << "Wrong amount of arguments\n";
		help(argv[0]);
		exit(0);
	}

    char receivedCharArr[100];
    string receivedStr;
    int length;


	MPI_Init(&argc, &argv);
	MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &cores);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        multipdbqtFileStringsConstructor MAINOUTPUT_pdbqtFileConstructor("AllResults.pdbqt");
        logSummaryConstructor MAINOUTPUT_logFileConstructor("AllLogs.log");
        multipdbqtFileStringsReader masterFile(multiPdbqtFileName, 0);
        multipdbqtFileStringsReader * vinaOutPdbqt;
        vinaResult vinaOutLog;
        pdbqtFileStrings file, vinaResFile;
        cout << rank << " I'm the BOSS\n";
        masterFile.getNextPdbqt(file);
        while (file.size())
        {
            if (file.outpdbqtFileStrings()){
                cout << "Error with file " << file.name << "\n";
            }
            MPI_Recv(receivedCharArr, 100, MPI_CHAR, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_CHAR, &length); //Getting exact length of the string
            if (length){

                receivedStr = string(receivedCharArr, length); //Constructing the string
                deleteFile(receivedStr + ".pdbqt");

                vinaOutPdbqt = new multipdbqtFileStringsReader(receivedStr + "_out.pdbqt", 0);
                vinaOutPdbqt->getNextPdbqt(vinaResFile);
                vinaOutPdbqt->close();
                MAINOUTPUT_pdbqtFileConstructor.add(vinaResFile);

                deleteFile(receivedStr + "_out.pdbqt");

                if (vinaOutLog.parseLogFile(receivedStr + "_out.log")){
                    cout << "Error " << receivedStr << "\n";
                }
                MAINOUTPUT_logFileConstructor.add(receivedStr, vinaOutLog);

                deleteFile(receivedStr + "_out.log");
            }


            cout << "New job " << file.name << " for process " << status.MPI_SOURCE << "\n";

            MPI_Send((char *)file.name.c_str(), file.name.size(), MPI_CHAR, status.MPI_SOURCE, 11, MPI_COMM_WORLD);




            masterFile.getNextPdbqt(file);

        }
        masterFile.close();
        int countCores = 0;
        cout << "No more jobs!\n";
        receivedStr = string("");
        while (countCores < cores - 1) // Zero is me =)
        {
            MPI_Recv(receivedCharArr, 100, MPI_CHAR, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_CHAR, &length); //Getting exact length of the string
            if (length){

                receivedStr = string(receivedCharArr, length); //Constructing the string
                deleteFile(receivedStr + ".pdbqt");

                vinaOutPdbqt = new multipdbqtFileStringsReader(receivedStr + "_out.pdbqt", 0);
                vinaOutPdbqt->getNextPdbqt(vinaResFile);
                vinaOutPdbqt->close();
                MAINOUTPUT_pdbqtFileConstructor.add(vinaResFile);

                deleteFile(receivedStr + "_out.pdbqt");

                if (vinaOutLog.parseLogFile(receivedStr + "_out.log")){
                    cout << "Error " << receivedStr << "\n";
                }
                MAINOUTPUT_logFileConstructor.add(receivedStr, vinaOutLog);

                deleteFile(receivedStr + "_out.log");
            }


            countCores++;
            MPI_Send((char *)receivedStr.c_str(), 0, MPI_CHAR, status.MPI_SOURCE, 11, MPI_COMM_WORLD);
        }
    }
    else
    {
        receivedStr = "";

        MPI_Send((char *)receivedStr.c_str(), receivedStr.size(), MPI_CHAR, 0, 11, MPI_COMM_WORLD); //Sending an ask for a new job
        MPI_Recv(receivedCharArr, 100, MPI_CHAR, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, &status); //Receiving a new job new job
        MPI_Get_count(&status, MPI_CHAR, &length); //Getting exact length of the string
        receivedStr = string(receivedCharArr, length); //Constructing the string the string

        string executeStr;
        while (receivedStr.size()) //If length is zero - no more jobs
        {
            cout << rank << " received ligand " << receivedStr <<  "\n";

            executeStr = programm;
            executeStr += " ";
            executeStr += restParams;
            executeStr += " --ligand ";
            executeStr += receivedStr;
            executeStr += ".pdbqt --log ";
            executeStr += receivedStr;
            executeStr += "_out.log --out ";
            executeStr += receivedStr;
            executeStr += "_out.pdbqt";

            system(executeStr.c_str()); //Running the job
            MPI_Send((char *)receivedStr.c_str(), receivedStr.size(), MPI_CHAR, 0, 11, MPI_COMM_WORLD); //Sending an ask for a new job

            MPI_Recv(receivedCharArr, 100, MPI_CHAR, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, &status); //Receiving a new job new job
            MPI_Get_count(&status, MPI_CHAR, &length); //Getting exact length of the string
            receivedStr = string(receivedCharArr, length); //Constructing the string the string
        }
    }
	MPI_Finalize();
	return 0;
}

