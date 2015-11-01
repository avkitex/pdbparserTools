#ifndef COMMON_FUNCS
#define COMMON_FUNCS

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;


bool in(char a, string s);

bool isdigit(char a);

double strtodoub(string s);

void parseFileName(string s, string &path, string &name, string &ext);

vector <string> split(string s, string chars);

int strtoint(string s);

string trim(string s);

bool isletter(char a);

bool isUpperLetter(char a);

bool isLowerLetter(char a);

string toLowerCase(string s);

string getDate();

string getTime();

#endif
