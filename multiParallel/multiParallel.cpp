//Last update 23.10.15 17:20
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <queue>

#include <unistd.h>
#include <ctime>

#define VERSION "1.0.1"

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

    bool outPdbqtFile(string filePath){
        if (!name.size()){
            return 1;
        }
        string fullFileName = filePath + name;
        ofstream fout(fullFileName.c_str());
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

class multiPdbqtFile{
    ifstream fileHandle;
    int modelsReturned;

    multiPdbqtFile(string file, int skipStr){
        fileHandle = new ifstream(file.c_str());
        modelsReturned = 0;
        while (getline(fileHandle, s) && skipStr > 0){
            if (s.size() >5 && s.substr(0, 6) == "ENDMDL"){
                modelsReturned++;
                skipStr--;
            }
        }
    }

    pdbqtFile getNextPdbqt(string fileName){
        pdbqtFile file;
        string s;
        while (getline(fileHandle, s)){
            if (s.size() >5 && s.substr(0, 6) == "ENDMDL"){
                if (file.size()){
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
    }
};

class freeCores{
	queue <int> cores;

	freeCores(){
        cores.clear();
	}

};


int main(int argc, char **argv)
{
	int rank, a, am_works = 0, cores;
	int t[2];
	time_t timer;
	string multiPdbqtfile, programName, restParams, s;

	vector <string> arguments;
	vector <string> path;
	if (argc == 4)
	{
		programm = argv[1];
        multiPdbqtfile = argv[2];
        restParams = argv[3];
	}
	else
	{
		cout << "Wrong amount of arguments\n";
		help(argv[0]);
		exit(0);
	}

	multiPdbqtFile masterFile(multiPdbqtfile, 0);
	pdbqtFile file;


	MPI_Init(&argc, &argv);
	MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD,&cores);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        for (int i = 0; i < arguments.size(); ++i)
        {
            q.push(i);
        }
        cout << rank << " I'm the BOSS\n";
        while (q.size() > 0)
        {
            MPI_Recv(t, 2, MPI_INT, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, &status);
            cout << "New job " << q.front() << " for process " << status.MPI_SOURCE << "\n";
            t[0] = q.front();
            q.pop();
            MPI_Send(&t, 2, MPI_INT, status.MPI_SOURCE, 11, MPI_COMM_WORLD);
        }
        while (q2.size() < cores - 1)
        {
            MPI_Recv(t, 2, MPI_INT, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, &status);
            cout << "No more jobs! " << status.MPI_SOURCE << "\n";
            q2.push(t[0]);
            t[0] = -1;
            MPI_Send(&t, 2, MPI_INT, status.MPI_SOURCE, 11, MPI_COMM_WORLD);
        }

    }
    else
    {
        time(&timer);

        t[0] = rank;
        MPI_Rsend(&t, 2, MPI_INT, 0, 11, MPI_COMM_WORLD);
        MPI_Recv(t, 2, MPI_INT, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, &status);
        while (t[0] != -1)
        {
            cout << rank << " received job number " << t[0] <<  "\n";
            if (chdir(path[t[0]].c_str()) == -1)
            {
                cout << rank << " Unable to change directory. Job terminated\n";
                exit(1);
            }
            s = programm + " " + arguments[t[0]];// + " 2>&1 >> " + itostr(rank) + "log.txt";
            //system("date");
            system(s.c_str());
            //system("date");

            t[0] = rank;
            MPI_Send(&t, 2, MPI_INT, 0, 11, MPI_COMM_WORLD);
            MPI_Recv(t, 2, MPI_INT, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, &status);
        //execl(programm.c_str(), programm.c_str(), arguments[n].c_str(), (char *)0);
        }
    }
	MPI_Finalize();
	return 0;
}

