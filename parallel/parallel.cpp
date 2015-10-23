/*
 * Compile & run:
 *
 * $ mpicxx <file.cpp> -o <binary>
 * $ mpirun -np <cores amount> -q query  ./<binary> programm_name list_with_arguments
 *
 */
//Last update 12.12.13 20:20
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <queue>

#include <unistd.h>
#include <ctime>

#define VERSION "1.0.7"

using namespace std;

void help(string p_name)
{
    cout << "Parallel program starter v" << VERSION << "\n";
    cout << "Usage:\n\n";
    cout << p_name << " programm list\n";
    cout << "programm - full path for the programm to be run\n";
    cout << "list - file, containig few strings, each string contains all parameters for 1 run of the programm\n";
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

int main(int argc, char **argv) {
// Type your code here

	// Get my number. the first is '0'
	int rank, a, am_works = 0, cores;
	int t[2];
	time_t timer;
	queue <int> q, q2;
	string file, programm, s, num;

	vector <string> arguments;
	vector <string> path;
	if (argc == 3)
	{
		programm = argv[1];
        file = argv[2];
	}
	else
	{
		cout << "Wrong amount of arguments\n";
		help(argv[0]);
		exit(0);
	}
	ifstream pin (file.c_str());
	while (getline(pin, s))
    {
        a = s.find('\t');
        path.push_back(s.substr(0, a));
        arguments.push_back(s.substr(a + 1, s.size() - a));
    }
    pin.close();

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

