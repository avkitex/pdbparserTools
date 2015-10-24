#include <iostream>

using namespace std;



bool isdigit(char a) {
	return (a >= '0' && a <= '9');
}

int strtodoub(string s)
{
	int ans = 0, i = 0;
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
		i++;
	}
	if (s[i] == '.' || s[i] == ',' )
	{
		i++;
	}
	while(i < s.size() && isdigit(s[i]))
	{
		ans *= 10;
		ans += s[i] - '0';
		i++;
	}
	return ans * (minus ? -1 : 1);
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


void warning(string function, string message) {
    if (!NOWARNINGS) {
        system("tput setf 3");
        cerr << "Warning!!! Function " << function << ". Msg: " << message;
        system("tput setf 7");
        cerr << "\n";
    }
	return;
}
