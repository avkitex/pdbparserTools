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
