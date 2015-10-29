#include <iostream>

#include "commonFuncs.h"


using namespace std;

bool in(char a, string s)
{
    for (int i = 0; i < s.size(); ++i)
    {
        if (a == s[i])
        {
            return 1;
        }
    }
    return 0;
}

bool isdigit(char a) {
	return (a >= '0' && a <= '9');
}

double strtodoub(string s)
{
	int i = 0, mantissa = 1;
	double result = 0;
	bool minus = false;
	while(i < s.size() && !isdigit(s[i]) && s[i] != '-')
	{
		++i;
	}
	if (s[i] == '-')
	{
		minus = true;
		++i;
	}
	while(i < s.size() && isdigit(s[i]))
	{
		result *= 10;
		result += int(s[i] - '0');
		i++;
	}
	if (s[i] == '.' || s[i] == ',' )
	{
		i++;
	}
	while(i < s.size() && isdigit(s[i]))
	{
	    mantissa *= 10;
		result *= 10;
		result += int(s[i] - '0');
		i++;
	}
	return result * (minus ? -1 : 1) / mantissa;
}

vector <string> split(string s, string chars)
{
    string curs = "";
    vector <string> res;
    res.clear();
    for (int i = 0; i < s.size(); ++i)
    {
        if (in(s[i], chars))
        {
            if (curs.size() > 0)
            {
                res.push_back(curs);
                curs = "";
            }
        }
        else
        {
            curs += s[i];
        }
    }
    if (curs.size() > 0)
    {
        res.push_back(curs);
    }
    return res;
}

void parseFileName(string s, string &path, string &name, string &ext)
{
    int i = s.size() - 1;
    ext = "";
    name = "";
    path = "";
    while (i >= 0 && s[i] != '/' && s[i] != '.')
    {
        --i;
    }
    if (s[i] == '.')
    {
        ext = s.substr(i + 1, s.size() - i - 1);
        s = s.substr(0, i);
    }
    while (i >= 0 && s[i] != '/')
    {
        --i;
    }
    if (s[i] == '/')
    {
        name = s.substr(i + 1, s.size() - i - 1);
        path = s.substr(0, i + 1);
    }
    if (i < 0)
    {
        name = s;
        path = "";
    }
    return;
}

int strtoint(string s)
{
	int i = 0, ans = 0;
	bool minus = false;
	while(i < s.size() && !isdigit(s[i]) && s[i] != '-')
	{
		++i;
	}
	if (s[i] == '-')
	{
		minus = true;
		++i;
	}
	while(i < s.size() && isdigit(s[i]))
	{
		ans *= 10;
		ans += s[i] - '0';
		++i;
	}
	return ans * (minus ? -1 : 1);
}

