//***********************************************************
// Author: Kotlov Nikita, FBB, MSU
//
#define VERSION "1.3.5"
// Last modification date: 22.05.2015 23:00
//***********************************************************

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//Bugs fixed:
//No connects to hetatm residues in cut mode - FIXED 25.11.13
//Inserted residues wrote in wrong way       - FIXED 25.11.13
//           105
//          105A
//           106
//
//During modifying names to add chain label
//searched for first dot, not last
//adding .pdb if need                         - FIXED 12.12.13
//wrong charges output                       - FIXED 13.01.14
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


//Что умеем:
//Парсить раздел данных с координатами и конектами
//(то, что начинается с ATOM и HETATM и CONECT)
//и сохранять это в структуру данных с довольно удобным доступом
//Выводить все, что введи в pdb файл в нужном формате
//Вырезать лиганд в отдельный файл
//Удалять ненужные лиганды
//Перенумеровывать сериалы атомов
//Перенумеровывать резидью номера (при выводе структуры) притом с указанного номера
//Добавлять теры(при выводе структуры)
//Сортировать атомы внутри резидью по шаблону
//Все строки, кроме тех, которые начинаются с ATOM, HETATM, CONNECT и TER игнорируются и не выводятся в новый файл.
//Резать на цепи
//Предсказание нестандартных лигандов
//Словарь резьдью (сортировка)
//Чтение и разбор dlg файлов (результаты работы autodock)
//Компоновка результатов в xls файл
//Добавлено исправление некоторых ошибок OPEN BABEL
//Добавлено исправление некоторых ошибок VMD
//Удалиние альтернативных положений атомов
//Замена с хлоре CL -> CLA
//Обновление словаря резьдью
//И еще куча всего за полгода
//*********************************************************
//!!!! Что в процессе написания???

//!!!! Добавление резидью в структуру
//!!!! Проверка на близосьть < 0.5 ангстрем

//!!!! Сортировка в зависимости от наличия водородов
//*********************************************************


#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

//Jpeg compressor library
#include "jpge.h"
#include "jpge.cpp"
//Jpeg compressor library

#define MAXWARNLIGCOR 3
#define MINIMUM_RESIDUES_IN_CHAIN 100

#define M_PI           3.14159265358979323846  /* pi */

#define PICTUREWIDTH 800
#define PICTUREHEIGHT 600
#define PICTURESIZE PICTUREWIDTH*PICTUREHEIGHT

using namespace jpge;
using namespace std;

//Декларации всяких нужных функций
string inttostr(int val, int length, bool direction);
int strtodoub(string s);
int strtoint(string s);
double lig_correcter_dif(string sample, string coord, string output);
string doubtostr(int val, int length, int point);
bool SILENT = false, NOWARNINGS = false, NO_DLG_OUT_ERRORS = false, POSSIBLE_BABEL = false;
bool RREN_CHLOR = false, RREN_ZINK = false, LIG_GEOMETRIC_CENTER = false, ADD_CHARGE = false, CORRECT_HETATM_CHAINS = false, REPLASE_EX_FILE = false;
bool NO_OC = false, NO_TF = false, NO_TER = false;
int no_charge = -1, ALTERNATIVE_COORDINATES = 0, VMDCORRECT = 0, IGNORE_SEGMENT = 0;
string SUFF = "";


int MAGIC = 20037;//DO NOT TOUCH!!!! DON'T WORK WITHOUT IT!!!!!!
template<typename CharT>
class DecimalSeparator : public std::numpunct<CharT>
{
    public:
        DecimalSeparator(CharT Separator)
        : m_Separator(Separator)
        {}

    protected:
        CharT do_decimal_point()const
        {
            return m_Separator;
        }

    private:
        CharT m_Separator;
};

struct atom
{
    char chain;//цепь
    string atom_name, element, resname, resnumb;//хим, элемент, название резидью
    char alt;//Альтернативные положения
    int serial;//Порядковый номер атома
    int x, y, z, oq, tf, ch_val, ch_sign;//Кооринаты,
    //x, y, z - 3 digits after dot
    //oq, tf - 2 digits after dot
//==================================================
//integer values are more comfortable, because
//we don't care about fidelity
//сохранение в целочисленном типе удобно,
//потому что не нужно беспокоиться о точности,
//а количество знаков постоянное
//==================================================
};

struct residue{
    vector <atom> array;
    string number;
};

struct segment{
    vector<residue> atom;
    vector<residue> hetatm;
    string cname;
};

struct chain
{
    vector <segment> array;
    char letter;
};

struct pdb_file
{
    vector <chain> chains;
    map <int, multiset<int> > conect;
    vector < string > info;
};

struct mol2bond
{
    int a, b;
    string type;
};


struct mol2file
{
    string mol_name, ch_type, mol_type;
    vector <atom> atoms;
    vector <mol2bond> bonds;
};



pdb_file pdbstruct;
map<string, vector<string> > res_base;


//#####################HEADERS########################
string inttostr_2(int n);
void correct_h_chains(pdb_file &);
//#####################HEADERS########################

void help(string parsername)
{

	cout << "PDB Parser. Version " << VERSION <<"\n\n";
    cout << "Usage:\n";
    cout << parsername << " -s input.pdb [-bf basefile.dat -hf hetfile.dat -os out.pdb -og ligand.pdb -rr 1 -rs -st -ct]\n";
    cout << "Warning! It will erase all data in output files.\n";
    cout << "   -s  file - input PDB file\n";
    cout << "   -os file - output PDB file\n";
    cout << "   -og file - output ligand PDB file\n";
    cout << "   -bf file - residue base file for sorting.\n";
    cout << "   	  Base residues file format:\n";
    cout << "   		 Residue name N.CA.O.CB.CG.OG1. (dot must be after each atom name)\n";
    cout << "   -ub - update BASEFILE, using current receptor structure. Requires -bf parameter.\n";
    cout << "   -hf file - HETATM file file (special file, where you can write, which HETFILE need to delete and which cut to ligand output file)\n";
    cout << "   -fasta file - make a fasta file with protein sequence\n";
    cout << "   	  HETATM file format:\n";
    cout << "   		 # - comment\n";
    cout << "   		 LIG FLP A #Ligand FLP from chain A will be written to ligand.pdb\n";
    cout << "   		 DEL HOH A #Hetatm HOH from chain A will be deleted and not written to output.pdb\n";
    cout << "   		 DEL XZF ! #Hetatm XZF from each (!) chain will be deleted and not written to output.pdb\n";
    cout << "   		 LIG FIND  #Program will search for non standard ligand\n";
    cout << "   		 DEL FIND !#Founded non standard ligand will be deleted from each chain\n";
    cout << "   		 REN HEME HEM !#Rename HEME residue as HEM\n";
    cout << "   		 REN HSD HIS !#Rename HSD residue as HIS\n";
    cout << "   		 CHA ZN2 ZN 2 !#Add charge +2 to atom Zn in res ZN2\n";
    cout << "   	  If there is LIG beginning string (USING LAST) in hetatm file, the -ol parameter must be\n";
    cout << "   -renres - rename residues (hetfile)\n";
    cout << "   -rr [num] - renumber residues (default 1)\n";
    cout << "   -rs [num] - renumber atom serial numbers (default 1)\n";
    cout << "   -st - sort residues, -bf parameter required\n";
    cout << "   -ct - cut structure to chains. N - amount of chains in PDB. Creating N receptor files (each file with 1 chain) and N ligand files (each file with 1 ligand) (if -og parameter set)\n";
    cout << "   -nc [num] - no charge. 0 - Charge will not be written to atoms with no charge in structure. 1 - No charge will be written at all. Default write existing.\n";
    cout << "   -silent - Silent mode. No comments will be shown.\n";
    cout << "   -cif file ligname - make ligand from mmCif.\n";
    cout << "   -nw - no warnings. No warnings will be shown.\n";
    cout << "   -dp file1 file2 - dlg parse parameter. File1 must contain a number of full path's to dlg files with docking results. Second - output file for parsed information (.txt.xls extension is recommended).\n";
    cout << "   -dprls file1 file2 - dlg parse parameter. File1 must be runlist. Second - output file for parsed information (.txt.xls extension is recommended).\n";
    cout << "   -vina_dp file1 file2 - vina log parse parameter. File1 must contain a number of full path's to vina log files with docking results. Second - output file for parsed information (.txt.xls extension is recommended).\n";
    cout << "   -vina_dprls file1 file2 - vina log parse parameter. File1 must be runlist. Second - output file for parsed information (.txt.xls extension is recommended).\n";
    cout << "   -ndoe - no dlg out error strings (if docking was not finished there won't be any strings in output file like \"Error there are only 4 of 100 results\").\n";
    cout << "   -pb - Possible OPEN BABEL errors. Some bugs such as double serial number on new connect string fixed. All atom hydrogens with name H will be treated as hetatm entries.\n";
    cout << "   -igseg - Ignore segment entries.\n";
    cout << "   -nooc - no occupancy output.\n";
    cout << "   -notf - no temperature factor out.\n";
    cout << "   -setoctf1 - set zero tf and oq to 1 (to avoid autopsf errors).\n";
    cout << "   -delallh - delete all hydrogens.\n";
    cout << "   -noter - no TER strings output.\n";
    cout << "   -vmd - Correcting pdb in NAMD format. Restoring chains according to segments. All residues in other segments except P or O will be deleted!!!.\n";
    cout << "   -clr - Rename clor ion residue and atom type. CL -> CLA.\n";
    cout << "   -zn2 - Rename zn residue and atom type. ZN -> ZN2.\n";
    cout << "   -ac - alternative coordinates. Deletes all except 1-st and remove alternative position labels.\n";
    cout << "   -ligcor file1 file2 file3 - Correct ligand using sample. Requires 3 files. Sample ligand file1, coordinate file2, output file3.\n";
    cout << "   -dlglig file1 [file2] [number] - Create ligand pdb from dlf files. If 1 parameter, new dlg for 1st cluster will be created in the directory, that contains dlg file with '.pdb' ending.\n";
    cout << "         If 2 parameters, input dlg and output pdb, new dlg for 1st cluster will be created with the output name and '.pdb' ending.\n";
    cout << "         If 3 parameters, input dlg and output pdb, new dlg for each cluster will be created with the output name and '#.pdb' ending (# - cluster number).\n";
    cout << "   -suf str - Add suffix. Adds suffix just before dot in every output filename.\n";
    cout << "   -evp file1 file2 file3 - With help of accordance of yes/no inhibitor (file1) and docking result (file2) make ROC plot (file3), calculate AUC value points and square.\n";
    cout << "   -chc - correct hetatm chains. For each hetatm residue with lack of chain or with chain that contains no atom residues gives new chain - the one, with with atom there is the nearest measure.\n";
    cout << "   -sdfc file dir - cut sdf to many mol2 files.\n";
    cout << "   -frep - replace existing files.\n";
    cout << "   -rotate num [num] [num] - rotate pdb coordinates on angle x, y, z (in that direct way).\n";
    cout << "   \n";
    cout << "   [arg] - argument is not required. If not given, default value will be used";
    cout << "\n\n\nCopyright Kotlov Nikita, kit.iz.179@gmail.com,  FBB MSU 2014\n";
	return;
}

void warning(string function, string message)
{
    if (!NOWARNINGS)
    {
        system("tput setf 3");
        cerr << "Warning!!! Function " << function << ". Msg: " << message;
        system("tput setf 7");
        cerr << "\n";
    }
	return;
}
void error_msg(string function, string message)
{
    system("tput setf 4");
    cerr << "ERROR!!! Function " << function << ". Msg: " << message;
    system("tput setf 7");
    cerr << "\n";
	return;
}

string to_upper(string s)
{
    string res = "";
    for (int i = 0; i < s.size(); ++i)
    {
        if (s[i] < 'z' && s[i] > 'a')
        {
            res += (s[i] + 'A' - 'a');
        }
        else
        {
            res += s[i];
        }
    }
}


string del_digits(string s)
{
    string res = "";
    for (int i = 0; i < s.size(); ++i)
    {
        if (!isdigit(s[i]))
        {
            res += s[i];
        }
    }
    return res;
}
string del_letters(string s)
{
    string res = "";
    for (int i = 0; i < s.size(); ++i)
    {
        if (isdigit(s[i]))
        {
            res += s[i];
        }
    }
    return res;
}

bool isdigit(char a)
{
	return (a >= '0' && a <= '9');
}

bool isint(string s)
{
    int i = 0;
    if (s[i] == '-')
    {
        i++;
    }
    for (i; i < s.size(); ++i)
    {
        if (s[i] < '0' || s[i] > '9')
        {
            return 0;
        }
    }
	return 1;
}

int min(int a, int b)
{
   return a <= b ? a : b;
}

bool isupletter(char a)
{
	return a <= 'Z' && a >= 'A';
}
bool isuplettersword(string s)
{
    for (int i = 0; i < s.size(); ++i)
    {
        if (!(s[i] <= 'Z' && s[i] >= 'A' || isdigit(s[i])))
        {
            return 0;
        }
    }
	return 1;
}

bool isletter(char a)
{
	return a <= 'Z' && a >= 'A' || a <= 'z' && a >= 'a';
}
bool file_exists(string fname)
{
    return ifstream(fname.c_str()) != NULL;
}
void filename_parse(string s, string &path, string &name, string &ext)
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

string get_file_name(string path, string name, string ext)
{
    int n;
    if (!REPLASE_EX_FILE)
    {
        if (file_exists(path + name + "." + ext))
        {
            n = 0;
            while(file_exists(path + name + inttostr_2(n) + "." + ext))
            {
                n++;
            }
            return path + name + inttostr_2(n) + "." + ext;
        }
        else
        {
            return path + name + "." + ext;
        }
    }
    else
    {
        return path + name + "." + ext;
    }
}


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

string trim(string s)
{
	int i = 0, j;
	while (i < s.size() && !isletter(s[i]))
	{
		i++;
	}
	j = i + 1;
	while (j < s.size() && (isletter(s[j]) || isdigit(s[j]) || s[j] == '\''))
	{
		j++;
	}
	return s.substr(i, j - i);
}

string trimdel(string s)
{
	int i = 0, j;
	while (i < s.size() && s[i] == ' ')
	{
		++i;
	}
	j = s.size() - 1;
	while (j > i && s[j] == ' ')
	{
		j--;
	}
	return s.substr(i, j - i + 1);
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
//Formatted string length - len
//direction 0 - left, 1 - right
string formatstr(string s, int len, bool direction)
{
	string answer = "";
	string scurf = "formatstr", smsg = "";
	if (s.size() > len)
	{
	    //smsg = "String length is more than in format parameter, returning full string";
		//warning(scurf, smsg);
		return s;
	}
	if (direction)
	{
		for (int i = 0; i < len - s.size(); ++i)
		{
			answer += " ";
		}
	}
	answer += s;
	if (!direction)
	{
		while (answer.size() < len)
		{
			answer += " ";
		}
	}
	return answer;
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

int find_chain(pdb_file &pdbstruct, char c)
{
    for (int i = 0; i < pdbstruct.chains.size(); ++i)
    {
        if (pdbstruct.chains[i].letter == c)
        {
            return i;
        }
    }
    return -1;
}
int find_segment(pdb_file &pdbstruct, int cur_chain, string seg)
{
    for (int i = 0; i < pdbstruct.chains[cur_chain].array.size(); ++i)
    {
        if (pdbstruct.chains[cur_chain].array[i].cname == seg)
        {
            return i;
        }
    }
    return -1;
}
int find_residue(pdb_file &pdbstruct, int cur_chain, int cur_segment, bool atom, string nresn)
{
    if (atom)
    {
        for (int i = 0; i < pdbstruct.chains[cur_chain].array[cur_segment].atom.size(); ++i)
        {
            if (pdbstruct.chains[cur_chain].array[cur_segment].atom[i].number == nresn)
            {
                return i;
            }
        }
    }
    else
    {
        for (int i = 0; i < pdbstruct.chains[cur_chain].array[cur_segment].hetatm.size(); ++i)
        {
            if (pdbstruct.chains[cur_chain].array[cur_segment].hetatm[i].number == nresn)
            {
                return i;
            }
        }
    }
    return -1;
}

void parse_pdb(string file, pdb_file & pdbstruct)
{
    string scurf = "parse_pdb", smsg = "";
	if (!file_exists(file))
    {
        smsg = "There is no file ";
        smsg += file;
        error_msg(scurf, smsg);
        exit(1);
    }
	ifstream pin (file.c_str());
	string s;
    bool begin = true, empty, warn = false, empty_chains = false;
	//Vars for coordinate section
	residue newresidue;
	atom buf;//Буферныя структура для атома в текущей строке
	segment newsegment;
	string nresn, seg;//Номер резидью
	chain newchain;//Пустая цепь
	//Vars for conect section
	int atom_to;//Номер атома, к которому добаляются коннекты
	int cur_chain, cur_seg, cur_res;
	multiset<int> a_ser_no;//Пустой сет для добавления в мап
	vector<int> babel;
	//Other
	int strno = 1;
	while (getline(pin, s))
	{
		if (s.substr(0, 6) == "ATOM  " || s.substr(0, 6) == "HETATM")
		{
            while (s.size() < 80)
            {
                s += ' ';
                warn = true;
            }
		    begin = false;
			buf.serial = strtoint(s.substr(6, 5));//Serial no
			buf.atom_name = trimdel(s.substr(12, 4));//Название атома
			buf.alt = s[16];//Альтернативное положение
			buf.resname = trimdel(s.substr(17, 4));//Название резидью

			buf.chain = s[21];//Название цепи
			if (!isupletter(buf.chain))
            {
                //buf.chain = 'A';
                empty_chains = true;
            }

			nresn = trimdel(s.substr(22, 5));//Номер резидью

			buf.x = strtodoub(s.substr(30, 8));//Координата x
			buf.y = strtodoub(s.substr(38, 8));//Координата y
			buf.z = strtodoub(s.substr(46, 8));//Координата z
			buf.oq = strtodoub(s.substr(54, 6));//Occupancy
			buf.tf = strtodoub(s.substr(60, 6));//Температурный фактор
			if (!IGNORE_SEGMENT)
            {
                seg = trimdel(s.substr(72, 4));
            }
            else
            {
                seg = "";
            }
			if (seg == "")
			{
			    seg = "EMPTY";
			}
			buf.element = trimdel(s.substr(76, 2));//Элемент
			if (s[78] == ' ' || s[79] == ' ')//Заряды
            {
                buf.ch_sign = -1;
                buf.ch_val = 0;
                //smsg = "There is no charge in string ";
                //smsg += trno;
				//warning(scurf, smsg);
            }
            else if ((s[78] == '-' || s[78] == '+') && isdigit(s[79]))
			{
                buf.ch_sign = (s[78] == '+' ? 1 : 0);
                buf.ch_val = strtoint("" + s[79]);
			}
			else if ((s[79] == '-' || s[79] == '+') && isdigit(s[78]))
			{
                buf.ch_sign = (s[79] == '+' ? 1 : 0);
                buf.ch_val = strtoint("" + s[78]);
			}
			else
			{
			    smsg = "Something wrong with charge in string ";
			    smsg += inttostr(strno, 5, true);
			    smsg += ". Chage ignored";
			    warning(scurf, smsg);
			}
            cur_chain = find_chain(pdbstruct, buf.chain);

			if (cur_chain == -1)//Если еще нет такой цепи, то добавить
			{
			    newchain.letter = buf.chain;
				pdbstruct.chains.push_back(newchain);
				cur_chain = pdbstruct.chains.size() - 1;
			}
            cur_seg = find_segment(pdbstruct, cur_chain, seg);
			if (cur_seg == -1)
			{
			    newsegment.cname = seg;
				pdbstruct.chains[cur_chain].array.push_back(newsegment);
				cur_seg = pdbstruct.chains[cur_chain].array.size() - 1;
			}
			if (s.substr(0, 6) == "ATOM  ")//Пишем атом - резидью
			{
                cur_res = find_residue(pdbstruct, cur_chain, cur_seg, 0, nresn);//Пытаемся найти хетатм с таким номером
			    if (POSSIBLE_BABEL && buf.atom_name == "H" && cur_res != -1)//Если возможны ошибки бабеля, нашелся хетатм и это водород, пишем в хетатм
                {
                    pdbstruct.chains[cur_chain].array[cur_seg].hetatm[cur_res].array.push_back(buf);
                }
                else
                {
                    cur_res = find_residue(pdbstruct, cur_chain, cur_seg, 1, nresn);
                    if (cur_res == -1)
                    {
                        newresidue.number = nresn;
                        pdbstruct.chains[cur_chain].array[cur_seg].atom.push_back(newresidue);
                        cur_res = pdbstruct.chains[cur_chain].array[cur_seg].atom.size() - 1;
                    }
                    pdbstruct.chains[cur_chain].array[cur_seg].atom[cur_res].array.push_back(buf);
                }
			}
			else if (s.substr(0, 6) == "HETATM")//Пишем hetatm - резидью
			{
			    cur_res = find_residue(pdbstruct, cur_chain, cur_seg, 0, nresn);
				if (cur_res == -1)
				{
					newresidue.number = nresn;
					pdbstruct.chains[cur_chain].array[cur_seg].hetatm.push_back(newresidue);
                    cur_res = pdbstruct.chains[cur_chain].array[cur_seg].hetatm.size() - 1;
				}
                pdbstruct.chains[cur_chain].array[cur_seg].hetatm[cur_res].array.push_back(buf);
			}
		}
		else if (s.substr(0, 6) == "CONECT")
		{
		    begin = false;
			if (strtoint(s.substr(6, min(s.size() - 6, 5))) == 0)
			{
			    smsg = "Wrong connection in string ";
			    smsg += inttostr(strno, 4, true);
			    warning(scurf, smsg);
				continue;
			}
			atom_to = strtoint(s.substr(6, 5));
			empty = true;
			for (int i = 11; i < min(s.size(), 16); ++i)
			{
                if (isdigit(s[i]))
                {
                    empty = false;
                }
			}
			if (empty)
			{
			    smsg = "Wrong connection in string ";
			    smsg += inttostr(strno, 4, true);
			    smsg += ". Only first atom listed";
			    if (!POSSIBLE_BABEL)
                {
                    warning(scurf, smsg);
                }
				continue;
			}
			if (strtoint(s.substr(11, min(s.size() - 11, 5))) == 0)
			{
			    smsg = "Wrong connection in string ";
			    smsg += inttostr(strno, 4, true);
			    smsg += ". There is zero serial atom number";
			    warning(scurf, smsg);
				continue;
			}
			a_ser_no.clear();
			if (!pdbstruct.conect.count(atom_to))
			{
				pdbstruct.conect[atom_to] = a_ser_no;
                for (int i = 11; i < 27 && strtoint(s.substr(i, min(s.size() - i, 5))); i += 5)
                {
                    pdbstruct.conect[atom_to].insert(strtoint(s.substr(i, min(s.size() - i, 5))));
                    if (pdbstruct.conect[atom_to].size() > 6)
                    {
                        smsg = "More than 6 connections with atom ";
                        smsg += inttostr(atom_to, 5, true);
                        smsg += " (string ";
                        smsg += inttostr(strno, 4, true);
                        smsg += ")";
                        warning(scurf, smsg);
                    }
                }
			}
			else
            {
                if (POSSIBLE_BABEL)
                {
                    babel.clear();
                    for (int i = 11; i < 27 && strtoint(s.substr(i, min(s.size() - i, 5))); i += 5)
                    {
                        babel.push_back(strtoint(s.substr(i, min(s.size() - i, 5))));
                    }
                    if (babel.size() == 2 && babel[0] == babel[1])
                    {
                        pdbstruct.conect[atom_to].insert(babel[0]);
                    }
                    else
                    {
                        for (int i = 0; i < babel.size(); ++i)
                        {
                            pdbstruct.conect[atom_to].insert(babel[i]);
                        }
                    }
                    if (pdbstruct.conect[atom_to].size() > 6)
                    {
                        smsg = "More than 6 connections with atom ";
                        smsg += inttostr(atom_to, 5, true);
                        smsg += " (string ";
                        smsg += inttostr(strno, 4, true);
                        smsg += ")";
                        warning(scurf, smsg);
                    }
                }
                else
                {
                    for (int i = 11; i < 27 && strtoint(s.substr(i, min(s.size() - i, 5))); i += 5)
                    {
                        pdbstruct.conect[atom_to].insert(strtoint(s.substr(i, min(s.size() - i, 5))));
                    }
                    if (pdbstruct.conect[atom_to].size() > 6)
                    {
                        smsg = "More than 6 connections with atom ";
                        smsg += inttostr(atom_to, 5, true);
                        smsg += " (string ";
                        smsg += inttostr(strno, 4, true);
                        smsg += ")";
                        warning(scurf, smsg);
                    }
                }
            }
		}
		else
		{
		    if (begin)
		    {
		        pdbstruct.info.push_back(s);
		    }
		}
		++strno;
	}
	if (warn && VMDCORRECT != 0)
    {
        smsg = "There were coordinate strings with length less than 80. They were filled with spaces";
        //warning(scurf, smsg);
    }
	if (!CORRECT_HETATM_CHAINS && empty_chains)
    {
        smsg = "There were coordinate strings with empty or wrong character chain fields";
        warning(scurf, smsg);
        CORRECT_HETATM_CHAINS = true;
    }
	return;
}

void out_coords(ofstream &sout, pdb_file &pdbstruct, bool atom_out)
{
	string s;
	int ccon;
	int ter_ser, ter_resnum;
	string ter_resn;
	char ter_ch;
	if (atom_out)
	{
        for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
        {
            for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
            {
                for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
                {
                    for (int i = 0; i < pdbstruct.chains[ch].array[seg].atom[res].array.size(); ++i)
                    {
                        s = "ATOM  ";
                        s += inttostr(pdbstruct.chains[ch].array[seg].atom[res].array[i].serial, 5, true);

                        s += " ";
                        if (pdbstruct.chains[ch].array[seg].atom[res].array[i].element.size() == 2 || pdbstruct.chains[ch].array[seg].atom[res].array[i].atom_name.size() == 4)
                        {
                            s += formatstr(pdbstruct.chains[ch].array[seg].atom[res].array[i].atom_name, 4, false);
                        }
                        else
                        {
                            s += " ";
                            s += formatstr(pdbstruct.chains[ch].array[seg].atom[res].array[i].atom_name, 3, false);
                        }
                        s += pdbstruct.chains[ch].array[seg].atom[res].array[i].alt;
                        s += formatstr(pdbstruct.chains[ch].array[seg].atom[res].array[i].resname, 4, false);
                        //По умолчанию строка с названием резидью из 3х символов
                        //и 4й пустой, но может быть и 2 и 4

                        s += pdbstruct.chains[ch].letter;
                        if (!isdigit(pdbstruct.chains[ch].array[seg].atom[res].number[pdbstruct.chains[ch].array[seg].atom[res].number.size() - 1]))
                        {
                            s += formatstr(pdbstruct.chains[ch].array[seg].atom[res].number, 5, true);
                        }
                        else
                        {
                            s += formatstr(pdbstruct.chains[ch].array[seg].atom[res].number, 4, true);
                            s += " ";
                        }
                        s += "   ";
                        s += doubtostr(pdbstruct.chains[ch].array[seg].atom[res].array[i].x, 8, 3);
                        s += doubtostr(pdbstruct.chains[ch].array[seg].atom[res].array[i].y, 8, 3);
                        s += doubtostr(pdbstruct.chains[ch].array[seg].atom[res].array[i].z, 8, 3);
                        if (NO_OC)
                        {
                            s += "      ";
                        }
                        else
                        {
                            s += doubtostr(pdbstruct.chains[ch].array[seg].atom[res].array[i].oq, 6, 2);
                        }
                        if (NO_TF)
                        {
                            s += "      ";
                        }
                        else
                        {
                            s += doubtostr(pdbstruct.chains[ch].array[seg].atom[res].array[i].tf, 6, 2);
                        }
                        s += "      ";
                        if (pdbstruct.chains[ch].array[seg].cname == "EMPTY")
                        {
                            s += "    ";
                        }
                        else
                        {
                            s += formatstr(pdbstruct.chains[ch].array[seg].cname, 4, false);
                        }
                        s += formatstr(pdbstruct.chains[ch].array[seg].atom[res].array[i].element, 2, true);
                        if (no_charge == -1 || (no_charge == 0 &&  pdbstruct.chains[ch].array[seg].atom[res].array[i].ch_sign != -1))
                        {
                            s += ('0' + pdbstruct.chains[ch].array[seg].atom[res].array[i].ch_val);
                            s += (pdbstruct.chains[ch].array[seg].atom[res].array[i].ch_sign == 1 ? "+" : "-");
                        }
                        else
                        {
                            s += "  ";
                        }
                        sout << s << "\n";
                    }
                }
            }
            if (pdbstruct.chains[ch].array[pdbstruct.chains[ch].array.size() - 1].atom.size() != 0 && !NO_TER)
            {
                s = "";
                s = "TER   ";
                s += inttostr(pdbstruct.chains[ch].array[pdbstruct.chains[ch].array.size() - 1].atom[pdbstruct.chains[ch].array[pdbstruct.chains[ch].array.size() - 1].atom.size() - 1].array[pdbstruct.chains[ch].array[pdbstruct.chains[ch].array.size() - 1].atom[pdbstruct.chains[ch].array[pdbstruct.chains[ch].array.size() - 1].atom.size() - 1].array.size() - 1].serial + 1, 5, true);

                s += "      ";
                s += formatstr(pdbstruct.chains[ch].array[pdbstruct.chains[ch].array.size() - 1].atom[pdbstruct.chains[ch].array[pdbstruct.chains[ch].array.size() - 1].atom.size() - 1].array[0].resname, 4, false);

                s += pdbstruct.chains[ch].letter;
                s += formatstr(pdbstruct.chains[ch].array[pdbstruct.chains[ch].array.size() - 1].atom[pdbstruct.chains[ch].array[pdbstruct.chains[ch].array.size() - 1].atom.size() - 1].number, 4, true);

                for (int i = s.size(); i < 80; ++i)
                {
                    s += " ";
                }
                sout << s << "\n";
            }
        }
	}
	for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
	{
	    for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
	    {
	        for (int res = 0; res < pdbstruct.chains[ch].array[seg].hetatm.size(); ++res)
	        {
	            for (int i = 0; i < pdbstruct.chains[ch].array[seg].hetatm[res].array.size(); ++i)
	            {

	                s = "HETATM";
                    s += inttostr(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].serial, 5, true);

					s += " ";
					if (pdbstruct.chains[ch].array[seg].hetatm[res].array[i].element.size() == 2 || pdbstruct.chains[ch].array[seg].hetatm[res].array[i].atom_name.size() == 4)
                    {
                        s += formatstr(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].atom_name, 4, false);
                    }
                    else
                    {
                        s += " ";
                        s += formatstr(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].atom_name, 3, false);
                    }
					s += pdbstruct.chains[ch].array[seg].hetatm[res].array[i].alt;
					s += formatstr(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].resname, 4, false);
					//По умолчанию строка с названием резидью из 3х символов
					//и 4й пустой, но может быть и 2 и 4

					s += pdbstruct.chains[ch].letter;
					if (!isdigit(pdbstruct.chains[ch].array[seg].hetatm[res].number[pdbstruct.chains[ch].array[seg].hetatm[res].number.size() - 1]))
                    {
                        s += formatstr(pdbstruct.chains[ch].array[seg].hetatm[res].number, 5, true);
                    }
                    else
                    {
                        s += formatstr(pdbstruct.chains[ch].array[seg].hetatm[res].number, 4, true);
                        s += " ";
                    }

					s += "   ";
					s += doubtostr(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].x, 8, 3);
					s += doubtostr(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].y, 8, 3);
					s += doubtostr(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].z, 8, 3);
					s += doubtostr(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].oq, 6, 2);
					s += doubtostr(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].tf, 6, 2);
					s += "      ";
					if (pdbstruct.chains[ch].array[seg].cname == "EMPTY")
					{
					    s += "    ";
					}
					else
					{
					    s += formatstr(pdbstruct.chains[ch].array[seg].cname, 4, false);
					}
					s += formatstr(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].element, 2, true);
					if (no_charge == -1 || (no_charge == 0 &&  pdbstruct.chains[ch].array[seg].hetatm[res].array[i].ch_sign != -1))
                    {
                        s += ('0' + pdbstruct.chains[ch].array[seg].hetatm[res].array[i].ch_val);
                        s += (pdbstruct.chains[ch].array[seg].hetatm[res].array[i].ch_sign == 1 ? "+" : "-");
                    }
                    else
                    {
                        s += "  ";
                    }
					sout << s << "\n";
	            }
	        }
	    }
	}
	for (map<int, multiset<int> >::iterator ch = pdbstruct.conect.begin(); ch != pdbstruct.conect.end(); ++ch)
    {
        ccon = 0;
        s = "CONECT";
        s += inttostr(ch->first, 5, true);
        for (multiset<int>::iterator it = pdbstruct.conect[ch->first].begin(); it != pdbstruct.conect[ch->first].end(); ++it)
        {
            if (ccon == 4)//По стандарту коннектов может быть только 4 числа после номера атома, с которым пишется конект
            {//Сия конструкция помогает писать по 4 числа (5 с первым) в строке
                for (int i = s.size(); i < 80; ++i)
                {
                    s += " ";
                }
                sout << s << "\n";

                ccon = 0;
                s = "CONECT";
                s += inttostr(ch->first, 5, true);
            }
            s += inttostr(*it, 5, true);
            ++ccon;
        }
        for (int i = s.size(); i < 80; ++i)
        {
            s += " ";
        }
        sout << s << "\n";
    }
	return;
}
void set_oc_to_1(pdb_file &pdbstruct)
{
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
	{
	    for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
	    {
	        for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
	        {
	            for (int i = 0; i < pdbstruct.chains[ch].array[seg].atom[res].array.size(); ++i)
	            {
	                if (pdbstruct.chains[ch].array[seg].atom[res].array[i].oq == 0)
                    {
                        pdbstruct.chains[ch].array[seg].atom[res].array[i].oq = 100;
                    }
	                if (pdbstruct.chains[ch].array[seg].atom[res].array[i].tf == 0)
                    {
                        pdbstruct.chains[ch].array[seg].atom[res].array[i].tf = 100;
                    }
	            }
	        }
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].hetatm.size(); ++res)
	        {
	            for (int i = 0; i < pdbstruct.chains[ch].array[seg].hetatm[res].array.size(); ++i)
	            {
	                if (pdbstruct.chains[ch].array[seg].hetatm[res].array[i].oq == 0)
                    {
                        pdbstruct.chains[ch].array[seg].hetatm[res].array[i].oq = 100;
                    }
	                if (pdbstruct.chains[ch].array[seg].hetatm[res].array[i].tf == 0)
                    {
                        pdbstruct.chains[ch].array[seg].hetatm[res].array[i].tf = 100;
                    }
	            }
	        }
	    }
	}
    return;
}
void delete_all_hydrogens(pdb_file &pdbstruct)
{
    vector <int> fdel;
	set <int> atom_for_del;
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
	{
	    for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
	    {
	        for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
	        {
	            fdel.clear();
	            for (int i = 0; i < pdbstruct.chains[ch].array[seg].atom[res].array.size(); ++i)
	            {
	                if (pdbstruct.chains[ch].array[seg].atom[res].array[i].element == "H" || pdbstruct.chains[ch].array[seg].atom[res].array[i].element == "h" )
                    {
                        fdel.push_back(i);
                        atom_for_del.insert(pdbstruct.chains[ch].array[seg].atom[res].array[i].serial);
                    }
	            }
	            for (int i = fdel.size() - 1; i >= 0; --i)
                {
                    pdbstruct.chains[ch].array[seg].atom[res].array.erase(pdbstruct.chains[ch].array[seg].atom[res].array.begin() + fdel[i]);
                }
	        }
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].hetatm.size(); ++res)
	        {
	            fdel.clear();
	            for (int i = 0; i < pdbstruct.chains[ch].array[seg].hetatm[res].array.size(); ++i)
	            {
	                if (pdbstruct.chains[ch].array[seg].hetatm[res].array[i].element == "H" || pdbstruct.chains[ch].array[seg].hetatm[res].array[i].element == "h" )
                    {
                        fdel.push_back(i);
                        atom_for_del.insert(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].serial);
                    }
	            }
	            for (int i = fdel.size() - 1; i >= 0; --i)
                {
                    pdbstruct.chains[ch].array[seg].hetatm[res].array.erase(pdbstruct.chains[ch].array[seg].hetatm[res].array.begin() + fdel[i]);
                }
	        }
	    }
	}
	map<int, multiset<int> > newcon;
    multiset<int> bufmap;
    for (map<int, multiset<int> >::iterator it = pdbstruct.conect.begin(); it != pdbstruct.conect.end(); ++it)
    {//Как всегда заново перекопируем то, что нужно оставить и заменим исходное
        if (atom_for_del.count(it->first) == 0)
        {//Ежели в нашем списке нету этого атома, то нужно проверить то, с чем он связан и перекопировать нужное
            bufmap.clear();
            for (multiset<int>::iterator de = it->second.begin(); de != it->second.end(); ++de)
            {
                if (atom_for_del.count(*de) == 0)
                {
                    bufmap.insert(*de);
                }
            }
            if (bufmap.size() != 0)
            {//Если все убилось, то не стоит писать пустую запись
                newcon[it->first] = bufmap;
            }
        }
    }
    pdbstruct.conect = newcon;//Заменить
    return;
}

//Function for deleting alternative atoms except first in dictionary order
//105A -> 105
//105B -> delete
void remove_alternative(pdb_file &pdbstruct)
{
    string scurf = "remove_alternative", smsg = "";
    map <string, int> atoms;
    vector <int> for_del;
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
            {
                atoms.clear();
                for_del.clear();
                for (int i = 0; i < pdbstruct.chains[ch].array[seg].atom[res].array.size(); ++i)
                {
                    if (pdbstruct.chains[ch].array[seg].atom[res].array[i].alt != ' ')
                    {
                        if (atoms.count(pdbstruct.chains[ch].array[seg].atom[res].array[i].atom_name))
                        {
                            if (pdbstruct.chains[ch].array[seg].atom[res].array[i].alt < pdbstruct.chains[ch].array[seg].atom[res].array[atoms[pdbstruct.chains[ch].array[seg].atom[res].array[i].atom_name]].alt)
                            {
                                for_del.push_back(atoms[pdbstruct.chains[ch].array[seg].atom[res].array[i].atom_name]);
                                atoms[pdbstruct.chains[ch].array[seg].atom[res].array[i].atom_name] = i;
                            }
                            else if (pdbstruct.chains[ch].array[seg].atom[res].array[i].alt == pdbstruct.chains[ch].array[seg].atom[res].array[atoms[pdbstruct.chains[ch].array[seg].atom[res].array[i].atom_name]].alt)
                            {
                                smsg = "Same alternative letters, last will be removed";
                                warning(scurf, smsg);
                                for_del.push_back(i);
                            }
                            else
                            {
                                for_del.push_back(i);
                            }
                        }
                        else
                        {
                            atoms[pdbstruct.chains[ch].array[seg].atom[res].array[i].atom_name] = i;
                        }
                    }
                }
                if (for_del.size() > 0)
                {
                    sort(for_del.begin(), for_del.end());
                    reverse(for_del.begin(), for_del.end());
                    for (int i = 0; i < for_del.size(); ++i)
                    {
                        pdbstruct.chains[ch].array[seg].atom[res].array.erase(pdbstruct.chains[ch].array[seg].atom[res].array.begin() + for_del[i]);
                    }
                }
                for (map<string, int>::iterator mp = atoms.begin(); mp != atoms.end(); ++mp)
                {
                    pdbstruct.chains[ch].array[seg].atom[res].array[mp->second].alt = ' ';
                }
            }
        }
    }
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].hetatm.size(); ++res)
            {
                atoms.clear();
                for_del.clear();
                for (int i = 0; i < pdbstruct.chains[ch].array[seg].hetatm[res].array.size(); ++i)
                {
                    if (pdbstruct.chains[ch].array[seg].hetatm[res].array[i].alt != ' ')
                    {
                        if (atoms.count(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].atom_name))
                        {
                            if (pdbstruct.chains[ch].array[seg].hetatm[res].array[i].alt < pdbstruct.chains[ch].array[seg].hetatm[res].array[atoms[pdbstruct.chains[ch].array[seg].hetatm[res].array[i].atom_name]].alt)
                            {
                                for_del.push_back(atoms[pdbstruct.chains[ch].array[seg].hetatm[res].array[i].atom_name]);
                                atoms[pdbstruct.chains[ch].array[seg].hetatm[res].array[i].atom_name] = i;
                            }
                            else if (pdbstruct.chains[ch].array[seg].hetatm[res].array[i].alt == pdbstruct.chains[ch].array[seg].hetatm[res].array[atoms[pdbstruct.chains[ch].array[seg].hetatm[res].array[i].atom_name]].alt)
                            {
                                smsg = "Same alternative letters, last will be removed";
                                warning(scurf, smsg);
                                for_del.push_back(i);
                            }
                            else
                            {
                                for_del.push_back(i);
                            }
                        }
                        else
                        {
                            atoms[pdbstruct.chains[ch].array[seg].hetatm[res].array[i].atom_name] = i;
                        }
                    }
                }
                if (for_del.size() > 0)
                {
                    sort(for_del.begin(), for_del.end());
                    reverse(for_del.begin(), for_del.end());
                    for (int i = 0; i < for_del.size(); ++i)
                    {
                        pdbstruct.chains[ch].array[seg].hetatm[res].array.erase(pdbstruct.chains[ch].array[seg].hetatm[res].array.begin() + for_del[i]);
                    }
                }
                for (map<string, int>::iterator mp = atoms.begin(); mp != atoms.end(); ++mp)
                {
                    pdbstruct.chains[ch].array[seg].hetatm[res].array[mp->second].alt = ' ';
                }
            }
        }
    }
    return;
}


void rename_res_add_charge(pdb_file &pdbstruct, map <string, string> &re_names, map < pair <string, string >, int> &charges)
{
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            //cout << "!";
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
            {
                for (int i = 0; i < pdbstruct.chains[ch].array[seg].atom[res].array.size(); ++i)
                {
                    if (re_names.count(pdbstruct.chains[ch].array[seg].atom[res].array[i].resname))
                    {
                        pdbstruct.chains[ch].array[seg].atom[res].array[i].resname = re_names[pdbstruct.chains[ch].array[seg].atom[res].array[i].resname];
                    }
                    if (charges.count(make_pair(pdbstruct.chains[ch].array[seg].atom[res].array[i].resname, pdbstruct.chains[ch].array[seg].atom[res].array[i].atom_name)))
                    {
                        pdbstruct.chains[ch].array[seg].atom[res].array[i].ch_val = abs(charges[make_pair(pdbstruct.chains[ch].array[seg].atom[res].array[i].resname, pdbstruct.chains[ch].array[seg].atom[res].array[i].atom_name)]);
                        pdbstruct.chains[ch].array[seg].atom[res].array[i].ch_sign = (charges[make_pair(pdbstruct.chains[ch].array[seg].atom[res].array[i].resname, pdbstruct.chains[ch].array[seg].atom[res].array[i].atom_name)] >= 0 ? 1 : 0);
                    }
                }
            }
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].hetatm.size(); ++res)
            {
                for (int i = 0; i < pdbstruct.chains[ch].array[seg].hetatm[res].array.size(); ++i)
                {
                    if (re_names.count(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].resname))
                    {
                        pdbstruct.chains[ch].array[seg].hetatm[res].array[i].resname = re_names[pdbstruct.chains[ch].array[seg].hetatm[res].array[i].resname];
                    }
                    if (charges.count(make_pair(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].resname, pdbstruct.chains[ch].array[seg].hetatm[res].array[i].atom_name)))
                    {
                        pdbstruct.chains[ch].array[seg].hetatm[res].array[i].ch_val = abs(charges[make_pair(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].resname, pdbstruct.chains[ch].array[seg].hetatm[res].array[i].atom_name)]);
                        pdbstruct.chains[ch].array[seg].hetatm[res].array[i].ch_sign = (charges[make_pair(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].resname, pdbstruct.chains[ch].array[seg].hetatm[res].array[i].atom_name)] >= 0 ? 1 : 0);
                    }
                }
            }
        }
    }
    return;
}

void correct_VMD(pdb_file &pdbstruct)
{
    string scurf = "correct_VMD", smsg = "";
    pdb_file corps;
    chain newchain;
    segment newsegment;
    newsegment.cname = "EMPTY";
    int n;
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        if (pdbstruct.chains[ch].letter == 'P')
        {
            corps.chains.resize(pdbstruct.chains[ch].array.size());
            for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
            {
                corps.chains[seg].letter = char('A' + strtoint(del_letters(pdbstruct.chains[ch].array[seg].cname)) - 1);
                corps.chains[seg].array.push_back(newsegment);
                corps.chains[seg].array[0].atom = pdbstruct.chains[ch].array[seg].atom;
            }
        }
        else if (pdbstruct.chains[ch].letter == 'O')
        {
            newchain.letter = ' ';
            corps.chains.push_back(newchain);
            corps.chains[corps.chains.size() - 1].array.push_back(newsegment);
            for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
            {
                for (int i = 0; i < pdbstruct.chains[ch].array[seg].atom.size(); ++i)
                {
                    corps.chains[corps.chains.size() - 1].array[0].hetatm.push_back(pdbstruct.chains[ch].array[seg].atom[i]);
                }
            }
            correct_h_chains(corps);
        }
        else
        {
            //newchain.letter = ' ';
            //corps.chains.push_back(newchain);
            for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
            {
                if (pdbstruct.chains[ch].array[seg].cname[0] == 'P' || pdbstruct.chains[ch].array[seg].cname[0] == 'p')
                {
                    newchain.letter = char('A' + strtoint(del_letters(pdbstruct.chains[ch].array[seg].cname)) - 1);
                    corps.chains.push_back(newchain);
                    corps.chains[corps.chains.size() - 1].array.push_back(newsegment);
                    corps.chains[corps.chains.size() - 1].array[0].atom = pdbstruct.chains[ch].array[seg].atom;
                }
                else if (pdbstruct.chains[ch].array[seg].cname[0] == 'O' || pdbstruct.chains[ch].array[seg].cname[0] == 'o')
                {
                    corps.chains[strtoint(del_letters(pdbstruct.chains[ch].array[seg].cname)) - 1].array[0].hetatm = pdbstruct.chains[ch].array[seg].atom;
                }
            }
            correct_h_chains(corps);
            /*
            for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
            {
                n = find_chain(corps, char('A' + strtoint(del_letters(pdbstruct.chains[ch].array[seg].cname)) - 1));
                if (n == -1)
                {
                    smsg = "There is chain with hetatm and no atom ";
                    smsg += pdbstruct.chains[ch].array[seg].cname;
                    smsg += "\n";
                    warning(scurf, smsg);
                    corps.chains.push_back(newchain);
                    n = corps.chains.size() - 1;
                    corps.chains[n].letter = char('A' + strtoint(del_letters(pdbstruct.chains[ch].array[seg].cname)) - 1);
                    corps.chains[n].array.push_back(newsegment);
                }
                corps.chains[n].array[0].hetatm = pdbstruct.chains[ch].array[seg].atom;
            }*/
        }
    }
    pdbstruct = corps;
    return;
}

//formatting int VAL to string with current LENGTH in current DIRECTION
//direction 0 - left, 1 - right
string inttostr(int val, int length, bool direction)
{
	string s = "";
	int k = val;
	bool minus;
	if (k < 0)
	{
		minus = true;
		k *= -1;
	}
	else if (k == 0)
    {
        s = "0";
    }
	else
	{
		minus = false;
	}
	char r;
	while (k > 0)
	{
		r = k % 10 + '0';
		k /= 10;
		s = r + s;
	}
	if (minus)
	{
		s = "-" + s;
	}
	while (s.size() < length)
	{
		if (direction)
		{
			s = " " + s;
		}
		else
		{
			s += " ";
		}
	}
	return s;
}
string inttostr_2 (int val)
{
    string s = "";
	int k = val;
	bool minus;
	if (k < 0)
	{
		minus = true;
		k *= -1;
	}
	else if (k == 0)
    {
        s = "0";
    }
	else
	{
		minus = false;
	}
	char r;
	while (k > 0)
	{
		r = k % 10 + '0';
		k /= 10;
		s = r + s;
	}
	if (minus)
	{
		s = "-" + s;
	}
	return s;
}

//Функция переводящая число из инта в формат с
//длиной length и point знаками после запятой
string doubtostr(int val, int length, int point)
{
	string s = "";
	int k = val;
	bool minus;
	if (k < 0)
	{
		minus = true;
		k *= -1;
	}
	else
	{
		minus = false;
	}
	char r;
	for (int i = 0; i < point; ++i)
	{
		r = k % 10 + '0';
		s = r + s;
		k /= 10;
	}
	s = "." + s;
	if (k == 0)
	{
		s = "0" + s;
	}
	else
	{
		while (k > 0)
		{
			r = k % 10 + '0';
			s = r + s;
			k /= 10;
		}
	}
	if (minus)
	{
		s = "-" + s;
	}
	while (s.size() < length)
	{
		s = " " + s;
	}
	return s;
}
string doubtostr_3(double a, int point)
{
    string s = "";
    bool minus = false;
    if (a < 0)
    {
        minus = true;
        a *= -1;
    }
    double b = a;
    char r;
    for (int i = 0; i < point; ++i)
    {
        b *= 10;
        r = (int(b) % 10) + '0';
        s += r;
    }
    s = inttostr_2(int(a) * (minus ? -1 : 1)) + "." + s;
    return s;
}

string doubtostr_2(int val, int point, char coma = ',')
{
    string s = "";
	int k = val;
	bool minus = false;
	if (k < 0)
	{
		minus = true;
		k *= -1;
	}
	char r;
	for (int i = 0; i < point; ++i)
	{
		r = k % 10 + '0';
		s = r + s;
		k /= 10;
	}
	s = coma + s;
	if (k == 0)
	{
		s = "0" + s;
	}
	else
	{
		while (k > 0)
		{
			r = k % 10 + '0';
			s = r + s;
			k /= 10;
		}
	}
	if (minus)
	{
		s = "-" + s;
	}
	return s;
}


set <string> del_het(pdb_file &pdbstruct, map<char, set<string> > hets)
{

    string scurf = "del_het", smsg = "";
    bool flig = false;
    int cur_chain;
	set <int> atom_for_del;
	set <string> deleted_hets;
    for (map<char, set<string> >::iterator del_chain = hets.begin(); del_chain != hets.end(); ++del_chain)
    {
        if (del_chain->first == '!')
        {

            for (set<string>::iterator del_res = del_chain->second.begin(); del_res != del_chain->second.end(); ++del_res)
            {
                for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
                {
                    for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
                    {
                        for (int res = 0; res < pdbstruct.chains[ch].array[seg].hetatm.size(); ++res)
                        {
                            if (pdbstruct.chains[ch].array[seg].hetatm[res].array[0].resname == *del_res)
                            {
                                deleted_hets.insert(*del_res);
                                for (int i = 0; i < pdbstruct.chains[ch].array[seg].hetatm[res].array.size(); ++i)
                                {
                                    atom_for_del.insert(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].serial);
                                }
                                pdbstruct.chains[ch].array[seg].hetatm.erase(pdbstruct.chains[ch].array[seg].hetatm.begin() + res);
                                res--;
                            }
                        }
                    }
                }
            }
        }
        else
        {
            cur_chain = find_chain(pdbstruct, del_chain->first);
            if (cur_chain == -1)
            {
                smsg = "There is no chain ";
                smsg += del_chain->first;
                smsg += " in pdb\n";
                warning(scurf, smsg);
            }
            else
            {
                for (set<string>::iterator del_res = del_chain->second.begin(); del_res != del_chain->second.end(); ++del_res)
                {
                    flig = false;
                    for (int seg = 0; seg < pdbstruct.chains[cur_chain].array.size(); ++seg)
                    {
                        for (int res = 0; res < pdbstruct.chains[cur_chain].array[seg].hetatm.size(); ++res)
                        {
                            if (pdbstruct.chains[cur_chain].array[seg].hetatm[res].array[0].resname == *del_res)
                            {
                                flig = true;
                                deleted_hets.insert(*del_res);
                                for (int i = 0; i < pdbstruct.chains[cur_chain].array[seg].hetatm[res].array.size(); ++i)
                                {
                                    atom_for_del.insert(pdbstruct.chains[cur_chain].array[seg].hetatm[res].array[i].serial);
                                }
                                pdbstruct.chains[cur_chain].array[seg].hetatm.erase(pdbstruct.chains[cur_chain].array[seg].hetatm.begin() + res);
                                break;
                            }
                        }
                        if (flig)
                        {
                            break;
                        }
                    }
                    if (!flig)
                    {
                        smsg = "There is no ligand ";
                        smsg += *del_res;
                        smsg += " in chain ";
                        smsg += del_chain->first;
                        smsg += " in pdb\n";
                        warning(scurf, smsg);
                    }
                }
            }
        }
    }
	map<int, multiset<int> > newcon;
    multiset<int> bufmap;
    for (map<int, multiset<int> >::iterator it = pdbstruct.conect.begin(); it != pdbstruct.conect.end(); ++it)
    {//Как всегда заново перекопируем то, что нужно оставить и заменим исходное
        if (atom_for_del.count(it->first) == 0)
        {//Ежели в нашем списке нету этого атома, то нужно проверить то, с чем он связан и перекопировать нужное
            bufmap.clear();
            for (multiset<int>::iterator de = it->second.begin(); de != it->second.end(); ++de)
            {
                if (atom_for_del.count(*de) == 0)
                {
                    bufmap.insert(*de);
                }
            }
            if (bufmap.size() != 0)
            {//Если все убилось, то не стоит писать пустую запись
                newcon[it->first] = bufmap;
            }
        }
    }
    pdbstruct.conect = newcon;//Заменить
	return deleted_hets;
}
void renum_res(pdb_file &pdbstruct, int res_num)
{
    int sres = res_num;
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
	{
	    res_num = sres;
	    for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
	    {
	        for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
	        {
	            if (res_num == 0)
                {
                    if (isdigit(pdbstruct.chains[ch].array[seg].atom[res].number[pdbstruct.chains[ch].array[seg].atom[res].number.size() - 1]))
                    {
                        res_num = strtoint(pdbstruct.chains[ch].array[seg].atom[res].number);
                    }
                    else
                    {
                        res_num = strtoint(pdbstruct.chains[ch].array[seg].atom[res].number.substr(0, pdbstruct.chains[ch].array[seg].atom[res].number.size() - 1));
                    }
                }
	            pdbstruct.chains[ch].array[seg].atom[res].number = inttostr_2(res_num);
	            ++res_num;
	        }
	    }
	    for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
	    {
	        for (int res = 0; res < pdbstruct.chains[ch].array[seg].hetatm.size(); ++res)
	        {
	            if (res_num == 0)
                {
                    if (isdigit(pdbstruct.chains[ch].array[seg].hetatm[res].number[pdbstruct.chains[ch].array[seg].hetatm[res].number.size() - 1]))
                    {
                        res_num = strtoint(pdbstruct.chains[ch].array[seg].hetatm[res].number);
                    }
                    else
                    {
                        res_num = strtoint(pdbstruct.chains[ch].array[seg].hetatm[res].number.substr(0, pdbstruct.chains[ch].array[seg].hetatm[res].number.size() - 1));
                    }
                }
	            pdbstruct.chains[ch].array[seg].hetatm[res].number = inttostr_2(res_num);
	            ++res_num;
	        }
	    }
	}
	return;
}

void renum_ser(pdb_file &pdbstruct, int serial, bool atom_ren)
{
	map<int, int> reconect;//мап для перенумерации конектов
	//Atom
	if (atom_ren)
	{
        for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
        {
            for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
            {
                for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
                {
                    for (int i = 0; i < pdbstruct.chains[ch].array[seg].atom[res].array.size(); ++i)
                    {
                        reconect[pdbstruct.chains[ch].array[seg].atom[res].array[i].serial] = serial;
                        pdbstruct.chains[ch].array[seg].atom[res].array[i].serial = serial;
                        //здесь проиcходит перенумерация, ее результат запоминается, далее переделываются конекты
                        serial++;
                    }
                }
            }
            if (pdbstruct.chains[ch].array[pdbstruct.chains[ch].array.size() - 1].atom.size() != 0 && !NO_TER)
            {
                ++serial;//Пустой сериал для тера
            }
        }
	}
	//Hetatm
	for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
	{
	    for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
	    {
	        for (int res = 0; res < pdbstruct.chains[ch].array[seg].hetatm.size(); ++res)
	        {
	            for (int i = 0; i < pdbstruct.chains[ch].array[seg].hetatm[res].array.size(); ++i)
	            {
	                reconect[pdbstruct.chains[ch].array[seg].hetatm[res].array[i].serial] = serial;
	                pdbstruct.chains[ch].array[seg].hetatm[res].array[i].serial = serial;
				    //здесь проиcходит перенумерация, ее результат запоминается, далее переделываются конекты
	                serial++;
	            }
	        }
	    }
	}
	//Переделываем конекты
	map<int, multiset<int> > newmap;//новый для конектов
	multiset<int> newms;//новый пустой мап для чиселок
	//Из-за наших продвинутых структур данных,
	//прийдется сделать новую на основе старой и пересохранить.
	newmap.clear();
	for (map<int, multiset<int> >::iterator ch = pdbstruct.conect.begin(); ch != pdbstruct.conect.end(); ++ch)
	{
	    newms.clear();
		for (multiset<int>::iterator it = pdbstruct.conect[ch->first].begin(); it != pdbstruct.conect[ch->first].end(); ++it)
		{
		    if (reconect.count(*it))
		    {
		        newms.insert(reconect[*it]);
		    }
		    else
		    {
		        newms.insert(*it);
		    }
		}
		if (reconect.count(ch->first))
		{
		    newmap[reconect[ch->first]] = newms;
		}
		else
		{
		    newmap[ch->first] = newms;
		}
	}
	pdbstruct.conect = newmap;//Пересохраняем
	return;
}


void out_end(ofstream &cout)
{
	cout << "END   ";
	return;
}

void out_info(ofstream &cout, pdb_file &pdbstruct)
{
    for(int i = 0; i < pdbstruct.info.size(); ++i)
    {
        cout << pdbstruct.info[i] << "\n";
    }
    return;
}
//Вывод pbd файла
void out_pdb(pdb_file &pdbstruct, string out_file, bool atom_out)
{
    string path, name, ext;
    filename_parse(out_file, path, name, ext);
    out_file = get_file_name(path, name, ext);
	ofstream cout (out_file.c_str());
	//out_info(cout, pdbstruct);
	out_coords(cout, pdbstruct, atom_out);
	//out_master
	out_end(cout);
	cout.close();
	return;
}

string next_lex(string s, int &num, char stop)
{
	if (num >= s.size())
	{
		return "";
	}
	int j = num, i = num;
	while(s.size() > i && s[i] == stop)
	{
		i++;
	}
	j = i;
	while (s.size() > j && s[j] != stop)
	{
		++j;
	}
	num = j + 1;
	return s.substr(i, j - i);
}
//Считывание и обработка хетатм файла
bool parse_het_file(string hetfile, vector<string> &std_ligs, map<char, set<string> >  &delhetlist, string &ligname, char &ligchain, map <string, string> &rename, map < pair < string, string > , int> &achs)
{
    string scurf = "parse_het_file", smsg = "";
    if (!file_exists(hetfile))
    {
        smsg = "There is no file ";
        smsg += hetfile;
        error_msg(scurf, smsg);
        exit(1);
    }
	ifstream hin (hetfile.c_str());
	string s, lex, ligand, l2;
	int i = 0, str = 1;
	bool lig = (ligname.size() ? 1:0);
	vector<char> delstdchain;
	set<string> newstr;
	while (getline(hin, s))
	{
		i = 0;
		lex = next_lex(s, i, ' ');
		if (lex.size() == 0)
		{
		    smsg = "Empty or space-filled string, continue working";
		    //warning(scurf, smsg);
		}
		else if (lex[0] == '#')
        {
            //Comment
        }
		else if (lex.size() != 3)
		{
		    smsg = "Something wrong with string type (length != 3) in string ";
		    smsg += inttostr_2(str);
		    smsg += ". Ignoring this string";
		    warning(scurf, smsg);
		}
		else if (lex == "LIG")
		{
			if (!lig)
			{
				lig = true;
			}
			else
			{
                smsg = "More than 1 LIG. String ";
                smsg += inttostr_2(str);
                smsg += ". Using last";
                warning(scurf, smsg);
			}
			lex = next_lex(s, i, ' ');
			if (lex.size() > 4 || lex.size() < 2)
			{
                smsg = "Very strange LIG ligand. String ";
                smsg += inttostr_2(str);
                warning(scurf, smsg);
			}
			ligname = lex;
			lex = next_lex(s, i, ' ');
			if (lex.size() < 1 || lex[0] == '#')
			{
			    ligchain = '?';
			}
			else if (lex.size() != 1 || !isupletter(lex[0]))
			{
                smsg = "Wrong chain in string ";
                smsg += inttostr_2(str);
                warning(scurf, smsg);
				ligchain = '?';
			}
			else
			{
				ligchain = lex[0];
			}
		}
		else if (lex == "DEL")
		{
			lex = next_lex(s, i, ' ');
			if (lex.size() < 2 || lex.size() > 4)
			{
                smsg = "Very strange ligand for DEL in string ";
                smsg += inttostr_2(str);
                warning(scurf, smsg);
			}
			if (lex == "STD")
			{
				lex = next_lex(s, i, ' ');
				if (lex.size() == 1 && (isupletter(lex[0]) || lex[0] == '!'))
				{
					delstdchain.push_back(lex[0]);
				}
				else
                {
                    smsg = "There is no chain in DEL STD in string ";
                    smsg += inttostr_2(str);
                    smsg += ". Ignoring this string";
                    warning(scurf, smsg);
                }
				continue;
			}
			ligand = lex;
			lex = next_lex(s, i, ' ');
			if (lex.size() == 1 && (isupletter(lex[0]) || lex[0] == '!'))
			{
				if (!delhetlist.count(lex[0]))
				{
					delhetlist[lex[0]] = newstr;
				}
				delhetlist[lex[0]].insert(ligand);
			}
			else
			{
                smsg = "There is no chain in DEL LIG in string ";
                smsg += inttostr_2(str);
                smsg += ". Ignoring this string";
                warning(scurf, smsg);
			}
		}
		else if (lex == "STD")
		{
			lex = next_lex(s, i, ' ');
			if (lex.size() > 4 || lex.size() < 2)
			{
                smsg = "Very strange STD ligand in string ";
                smsg += inttostr_2(str);
                warning(scurf, smsg);
			}
			std_ligs.push_back(lex);
		}
		else if (lex == "REN")
        {
            lex = next_lex(s, i, ' ');
			if (lex.size() > 4 || lex.size() < 2)
			{
                smsg = "Very strange 1st REN ligand in string ";
                smsg += inttostr_2(str);
                warning(scurf, smsg);
			}
			ligand = lex;
            lex = next_lex(s, i, ' ');
			if (lex.size() > 4 || lex.size() < 2)
			{
                smsg = "Very strange 2nd REN ligand in string ";
                smsg += inttostr_2(str);
                warning(scurf, smsg);
			}
			rename[ligand] = lex;
        }
        else if (lex == "CHA")
        {
            lex = next_lex(s, i, ' ');
			if (lex.size() > 4 || lex.size() < 2)
			{
                smsg = "Very strange 1st CHA resname in string ";
                smsg += inttostr_2(str);
                warning(scurf, smsg);
			}
			ligand = lex;
            l2 = next_lex(s, i, ' ');
			if (lex.size() > 4)
			{
                smsg = "Atom type name must be 1-4 letters long ";
                smsg += inttostr_2(str);
                warning(scurf, smsg);
			}
			lex = next_lex(s, i, ' ');
			if (!isint(lex))
            {
                smsg = lex;
                smsg += " is not charge. Ignoring it.";
                warning(scurf, smsg);
            }
            else
            {
                achs[make_pair(ligand, l2)] = strtoint(lex);
            }

        }
		else
		{
            smsg = "Wrong beginning of string ";
            smsg += inttostr_2(str);
            smsg += ". Ignoring this string";
            warning(scurf, smsg);
			continue;
		}
		lex = next_lex(s, i, ' ');
		if (lex.size() != 0 && lex[0] != '#')
		{
            smsg = "Wrong end of string ";
            smsg += inttostr_2(str);
            smsg += ". Ignoring end of this string";
            warning(scurf, smsg);
		}
		str++;
	}
	for (int i = 0; i < delstdchain.size(); ++i)
	{
		for (int j = 0; j < std_ligs.size(); ++j)
		{
			delhetlist[delstdchain[i]].insert(std_ligs[j]);
		}
	}
	return lig;
}
//Сия функция копирует из pdbstruct нужный лиганд
//-1 - нет такой цепи
// 0 - нет такого лиганда в этой цепи
// 1 - все норм
int make_lig(pdb_file &ligand, string ligname, char ligchain, pdb_file &pdbstruct)
{
    pdb_file free;
    ligand = free;
    int cur_chain = find_chain(pdbstruct, ligchain);
    if (cur_chain == -1)
    {
        return -1;
    }
	chain newchain;
	segment newsegment;
    for (int seg = 0; seg < pdbstruct.chains[cur_chain].array.size(); ++seg)
    {
        for (int res = 0; res < pdbstruct.chains[cur_chain].array[seg].hetatm.size(); ++res)
        {
            if (pdbstruct.chains[cur_chain].array[seg].hetatm[res].array[0].resname == ligname)
            {
                newchain.letter = ligchain;
                ligand.chains.push_back(newchain);
                newsegment.cname = pdbstruct.chains[cur_chain].array[seg].cname;
                ligand.chains[0].array.push_back(newsegment);
                ligand.chains[0].array[0].hetatm.push_back(pdbstruct.chains[cur_chain].array[seg].hetatm[res]);
                for (int i = 0; i < ligand.chains[0].array[0].hetatm[0].array.size(); ++i)
                {
                    if (pdbstruct.conect.count(ligand.chains[0].array[0].hetatm[0].array[i].serial))
                    {
                        ligand.conect[ligand.chains[0].array[0].hetatm[0].array[i].serial] = pdbstruct.conect[ligand.chains[0].array[0].hetatm[0].array[i].serial];
                    }
                }
                return 1;
            }
        }
    }
	return 0;//такого лиганда не нашлось
}
//Вводится база резидью
void parse_res_base(map<string, vector<string> > &res_base, string filein)
{
    string scurf = "parse_res_base", smsg = "";
    if (!file_exists(filein))
    {
        smsg = "There is no file ";
        smsg += filein;
        error_msg(scurf, smsg);
        exit(1);
    }
	ifstream bin(filein.c_str());
	string s, cur_res, cur_atom;
	int i;
	res_base.clear();
	while(getline(bin, s))
    {
        i = s.find(' ');
        if (i != string::npos)
        {
            cur_res = s.substr(0, i);
            if (cur_res.size() > 4)
            {
                smsg = "Residue with name " ;
                smsg += cur_res;
                smsg += " has very long name size > 4";
                warning(scurf, smsg);
            }
            i++;
            cur_atom = "";
            while (i < s.size())
            {
                if (s[i] == '.')
                {
                    if (cur_atom.size() > 0 && cur_atom.size() <= 4)
                    {
                        res_base[cur_res].push_back(cur_atom);
                        cur_atom = "";
                    }
                    else if (cur_atom.size() > 4)
                    {
                        smsg = "Wrong atom type length ";
                        smsg += cur_atom;
                        warning(scurf, smsg);
                        cur_atom = "";
                    }
                    else
                    {
                        smsg = "Zero length atom name or unexpected dot";
                        warning(scurf, smsg);
                        cur_atom = "";
                    }
                }
                else
                {
                    cur_atom += s[i];
                }
                i++;
            }
            if (cur_atom.size() > 0 && cur_atom.size() <= 4)
            {
                res_base[cur_res].push_back(cur_atom);
                cur_atom = "";
            }
            else
            {
                smsg = "Wrong atom type length at the end of the string";
                warning(scurf, smsg);
            }
        }
        else
        {
            smsg = "Base file has wrong format. There is no space after residue name. Current string was skipped";
            warning(scurf, smsg);
        }
    }
    bin.close();
    return;
}
//Сортировка атомов внутри резидью
void sort_res(vector <atom> &res, vector <string> order, string res_name)
{
    string scurf = "sort_res", smsg = "";
    int arr[res.size()], k;
    bool p;
    for (int i = 0; i < res.size(); ++i)
    {
        for (k = 0; k < order.size(); ++k)
        {
            if (res[i].atom_name == order[k])
            {
                arr[i] = k;
                break;
            }
        }
        if (k == order.size())
        {
            smsg = "There is no atom ";
            smsg += res[i].atom_name;
            smsg += " for residue ";
            smsg += res_name;
            smsg += " in base file";
            warning(scurf, smsg);
            arr[i] = order.size();
        }
    }
    p = true;
    for (int i = 0; i < res.size() && p; ++i)
    {
        p = false;
        for (int j = 1; j < res.size(); ++j)
        {
            if (arr[j - 1] > arr[j])
            {
                p = true;
                swap(arr[j - 1], arr[j]);
                swap(res[j - 1], res[j]);
            }
        }
    }
    return;
}
//Сортировка атомов внутри структуры
void sort_pdb(pdb_file &pdbstruct, map <string, vector<string> > &res_base)
{
    string scurf = "sort_res", smsg = "";
	int j, count;
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
            {
                if (res_base.count(pdbstruct.chains[ch].array[seg].atom[res].array[0].resname))
                {
                    sort_res(pdbstruct.chains[ch].array[seg].atom[res].array, res_base[pdbstruct.chains[ch].array[seg].atom[res].array[0].resname], pdbstruct.chains[ch].array[seg].atom[res].array[0].resname);
                }
                else
                {
                    smsg = "There is no residue ";
                    smsg += pdbstruct.chains[ch].array[seg].atom[res].array[0].resname;
                    smsg += " in base file";
                    //warning(scurf, smsg);
                }
            }
        }
    }
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].hetatm.size(); ++res)
            {
                if (res_base.count(pdbstruct.chains[ch].array[seg].hetatm[res].array[0].resname))
                {
                    sort_res(pdbstruct.chains[ch].array[seg].hetatm[res].array, res_base[pdbstruct.chains[ch].array[seg].hetatm[res].array[0].resname], pdbstruct.chains[ch].array[seg].hetatm[res].array[0].resname);
                }
                else
                {
                    smsg = "There is no residue ";
                    smsg += pdbstruct.chains[ch].array[seg].hetatm[res].array[0].resname;
                    smsg += " in base file";
                    //warning(scurf, smsg);
                }
            }
        }
    }
	return;
}

vector <string> order_atom(vector <atom> &res)
{
    vector <string> order;
    for (int i = 0; i < res.size(); ++i)
    {
        order.push_back(res[i].atom_name);
    }
    return order;
}
//Обновление базы порядка атомов в резидью
void update_base(pdb_file &pdbstruct, map <string, vector<string> > &res_base)
{
    string scurf = "update_base", smsg = "";
    vector <string> order, buf;
    bool cb = false;
    int ressize;
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            ressize = pdbstruct.chains[ch].array[seg].atom.size() - 1;
            for (int res = 1; res < ressize; ++res)//Концевые атомы НЕ учитываются
            {
                if (!res_base.count(pdbstruct.chains[ch].array[seg].atom[res].array[0].resname))//Если такой резидью еще не в базе - записать
                {
                    res_base[pdbstruct.chains[ch].array[seg].atom[res].array[0].resname] = order_atom(pdbstruct.chains[ch].array[seg].atom[res].array);
                }
                else//Иначе понять, такой же он или нет
                {
                    order = order_atom(pdbstruct.chains[ch].array[seg].atom[res].array);
                    //cout << order.size() << " " << res_base[pdbstruct.chains[ch].array[seg].atom[res].array[0].resname].size() << "\n";
                    if (pdbstruct.chains[ch].array[seg].atom[res].array[0].resname == "CYS" && abs(int(order.size() - res_base[pdbstruct.chains[ch].array[seg].atom[res].array[0].resname].size())) == 1)
                    {//Если это цистеин, то из-за цистеиновых мостиков могут отличаться на 1 количества атомов
                        if (order.size() - res_base[pdbstruct.chains[ch].array[seg].atom[res].array[0].resname].size() == 1)
                        {
                            res_base[pdbstruct.chains[ch].array[seg].atom[res].array[0].resname] = order;
                            cb = true;
                        }
                    }
                    else if (order != res_base[pdbstruct.chains[ch].array[seg].atom[res].array[0].resname])
                    {
                        buf = res_base[pdbstruct.chains[ch].array[seg].atom[res].array[0].resname];
                        smsg = "Different atom order for atom residue ";
                        smsg += pdbstruct.chains[ch].array[seg].atom[res].array[0].resname;
                        smsg += " residue number ";
                        smsg += pdbstruct.chains[ch].array[seg].atom[res].number;
                        smsg += ", chain ";
                        smsg += pdbstruct.chains[ch].letter;
                        sort(order.begin(), order.end());
                        sort(buf.begin(), buf.end());
                        if (order.size() > res_base[pdbstruct.chains[ch].array[seg].atom[res].array[0].resname].size())//Резидью отличаются количеством атомов
                        {
                            smsg += ". They have different amout of atoms. The longest was saved";
                            res_base[pdbstruct.chains[ch].array[seg].atom[res].array[0].resname] = order;
                        }
                        else if (order != buf)//Резидью отличаются атомами
                        {
                            smsg += ". They have different atom names. The previous was saved";
                        }
                        else
                        {
                            //Резидью отличаются только порядком атомов
                        }
                        warning(scurf, smsg);
                    }
                }
            }
        }
    }//Аналагично для хетатмов, только здесь концевые атомы учитываются
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].hetatm.size(); ++res)
            {
                if (!res_base.count(pdbstruct.chains[ch].array[seg].hetatm[res].array[0].resname))
                {
                    res_base[pdbstruct.chains[ch].array[seg].hetatm[res].array[0].resname] = order_atom(pdbstruct.chains[ch].array[seg].hetatm[res].array);
                }
                else
                {
                    order = order_atom(pdbstruct.chains[ch].array[seg].hetatm[res].array);
                    if (order != res_base[pdbstruct.chains[ch].array[seg].hetatm[res].array[0].resname])
                    {
                        buf = res_base[pdbstruct.chains[ch].array[seg].hetatm[res].array[0].resname];
                        smsg = "Different atom order for hetatm residue ";
                        smsg += pdbstruct.chains[ch].array[seg].hetatm[res].array[0].resname;
                        smsg += " residue number ";
                        smsg += pdbstruct.chains[ch].array[seg].hetatm[res].number;
                        smsg += ", chain ";
                        smsg += pdbstruct.chains[ch].letter;
                        sort(order.begin(), order.end());
                        sort(buf.begin(), buf.end());
                        if (order.size() > res_base[pdbstruct.chains[ch].array[seg].hetatm[res].array[0].resname].size())
                        {
                            smsg += ". They have different amout of atoms. The longest was saved";
                            res_base[pdbstruct.chains[ch].array[seg].hetatm[res].array[0].resname] = order;
                        }
                        else if (order != buf)
                        {
                            smsg += ". They have different atom names. The previous was saved";
                        }
                        warning(scurf, smsg);
                    }
                }
            }
        }
    }
    if (cb)
    {
        smsg = "There are cysteins with different amount of atoms. Possible cystein bridges";
        warning(scurf, smsg);
    }
	return;
}
void out_res_base(map<string, vector<string> > &res_base, string fileout)
{
    string scurf = "out_res_base", smsg = "";
    string path, name, ext;
    filename_parse(fileout, path, name, ext);
    fileout = get_file_name(path, name, ext);
	ofstream bout(fileout.c_str());
	string s, msg;
	for (map<string, vector<string> >::iterator it = res_base.begin(); it != res_base.end(); ++it)
    {
        s = it->first;
        s += " ";
        if (it->second.size() > 0)
        {
            s += it->second[0];
            for (int i = 1; i < it->second.size(); ++i)
            {
                s += '.';
                s += it->second[i];
            }
            bout << s << "\n";
        }
        else
        {
            smsg = "There is no atoms in residue ";
            smsg += it->first;
            warning(scurf, smsg);
        }
    }
    bout.close();
	return;
}
int find_lig(pdb_file &pdbstruct, vector <string> stdligands, string &ligand, char &chain)
{
	string scurf = "find_lig", smsg = "";
	vector <string> all_ligs;
	int lig_find = 0;
	int cur_chain;
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
	{
        if (chain == '?')
        {
            cur_chain = ch;
        }
        else
        {
            cur_chain = find_chain(pdbstruct, chain);
        }

	    for (int seg = 0; seg < pdbstruct.chains[cur_chain].array.size(); ++seg)
	    {
	        for (int res = 0; res < pdbstruct.chains[cur_chain].array[seg].hetatm.size(); ++res)
	        {
	            if (find(stdligands.begin(), stdligands.end(), pdbstruct.chains[cur_chain].array[seg].hetatm[res].array[0].resname) == stdligands.end())
	            {
	                if (lig_find)
	                {
	                    smsg = "More than 1 non-standard ligand: ";
	                    smsg += ligand;
	                    smsg += " in chain ";
	                    smsg += chain;
	                    smsg += " and ";
	                    smsg += pdbstruct.chains[cur_chain].array[seg].hetatm[res].array[0].resname;
	                    smsg += " in chain ";
	                    smsg += pdbstruct.chains[cur_chain].letter;
	                    smsg += ". Using first as main non-standard ligand";
	                    warning(scurf, smsg);
	                    return 2;
	                }
	                else
	                {
	                    ligand = pdbstruct.chains[cur_chain].array[seg].hetatm[res].array[0].resname;
	                    chain = pdbstruct.chains[cur_chain].letter;
	                    lig_find = 1;
	                }
	            }
	        }
	    }
	    if (chain != '?')
	    {
	        break;
	    }
	}
    return lig_find;
}

void geometric_coords(pdb_file &ligand, string lgc_file, string resname)
{
	string scurf = "geometric_coords", smsg = "";
    string path, name, ext;
    filename_parse(lgc_file, path, name, ext);
    lgc_file = get_file_name(path, name, ext);
	ofstream lout(lgc_file.c_str());
	long x = 0, y = 0, z = 0, n = 0;
    for (int ch = 0; ch < ligand.chains.size(); ++ch)
    {
        for (int seg = 0; seg < ligand.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < ligand.chains[ch].array[seg].hetatm.size(); ++res)
            {
                for (int i = 0; i < ligand.chains[ch].array[seg].hetatm[res].array.size(); ++i)
                {
                    x += ligand.chains[ch].array[seg].hetatm[res].array[i].x;
                    y += ligand.chains[ch].array[seg].hetatm[res].array[i].y;
                    z += ligand.chains[ch].array[seg].hetatm[res].array[i].z;
                    n++;
                }
            }
        }
    }
    lout << "Residue name\t" << resname << "\n";
    //cout << "\n" << x / 1000. / n << " " << doubtostr_3(x /1000. / n, 3)  << " " << y / 1000. / n << " " << z / 1000. / n << " " << n << "\n";
    lout << "Center coordinates\t" << doubtostr_3(x /1000. / n, 3) << " " << doubtostr_3(y / 1000. / n, 3) << " " << doubtostr_3(z / 1000. / n, 3);

	lout.close();
	return;
}
void geometric_coords_mf(mol2file &ligand, string lgc_file, string resname)
{
	string scurf = "geometric_coords", smsg = "";
    string path, name, ext;
    filename_parse(lgc_file, path, name, ext);
    lgc_file = get_file_name(path, name, ext);
	ofstream lout(lgc_file.c_str());
	long x = 0, y = 0, z = 0, n = 0;
    for (int a = 0; a < ligand.atoms.size(); ++a)
    {

        x += ligand.atoms[a].x;
        y += ligand.atoms[a].y;
        z += ligand.atoms[a].z;
        n++;
    }
    lout << "Residue name\t" << resname << "\n";
    //cout << "\n" << x / 1000. / n << " " << doubtostr_3(x /1000. / n, 3)  << " " << y / 1000. / n << " " << z / 1000. / n << " " << n << "\n";
    lout << "Center coordinates\t" << doubtostr_3(x /10000. / n, 3) << " " << doubtostr_3(y / 10000. / n, 3) << " " << doubtostr_3(z / 10000. / n, 3);

	lout.close();
	return;
}

void cutting_lig(pdb_file &pdbstruct, string &ligname, string ligout, vector <string> stdligands, int renumb_ser, int renumb_res, bool del)
{
    string scurf = "cutting_lig", smsg = "";
	pdb_file ligand;
	string postf = "", cur_lig_name = ligname;
	map < char, set < string > > buf;
	int k = ligout.size() - 1, cur_res, call = 1;
	while (k > 0 && ligout[k] != '.')
	{
		k--;
	}
	if (k == 0)
    {
        postf = ".pdb";
        k = ligout.size();
    }
    else
    {
        postf = ligout.substr(k, ligout.size() - k);
    }
    if (!SILENT) cout << "";
	for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
	{
	    if (cur_lig_name == "FIND")
	    {
            if (!SILENT) cout << "  Trying to FIND non-standard ligand for chain " << pdbstruct.chains[ch].letter << "\n";
            call = find_lig(pdbstruct, stdligands, ligname, pdbstruct.chains[ch].letter);
            if (call == 0)
            {
                smsg = "There is no non-standard ligands in chain ";
                smsg += pdbstruct.chains[ch].letter;
                warning(scurf, smsg);
            }
            else if (call == 1)
            {
                if (!SILENT) cout << "  There is 1 non-standard ligand " << ligname << " in chain " << pdbstruct.chains[ch].letter << "\n";
            }
            else if (call == 2)
            {
                if (!SILENT) cout << "  There is more than 1 non-standard ligand in chain " << pdbstruct.chains[ch].letter << "\n";
            }
            else
            {
                smsg = "Wired callback from function find_lig";
                warning(scurf, smsg);
            }
	    }
	    if (call == 1 || call == 2)
        {
            if (!SILENT) cout << "  Making " << ligname << " structure from chain " << pdbstruct.chains[ch].letter <<  "\n";
            if (make_lig(ligand, ligname, pdbstruct.chains[ch].letter, pdbstruct) == 0)
            {
                smsg = "There is no ligand ";
                smsg += ligname;
                smsg += " in chain ";
                smsg += pdbstruct.chains[ch].letter;
                warning(scurf, smsg);
            }
            else
            {
                if (!SILENT) cout << "  " << ligname << " structure from chain " << pdbstruct.chains[ch].letter << " done successfully\n";
                if (LIG_GEOMETRIC_CENTER)
                {
                    if (!SILENT) cout << "  Writing geometric coordinates of the ligand\n";
                    geometric_coords(ligand, ligout.substr(0, k) + "_" +  pdbstruct.chains[ch].letter +"_center.txt", ligname);
                }
                if (!SILENT) cout << "  Writing ligand file with";
                if (renumb_ser > 0)
                {
                    renum_ser(ligand, 1, false);
                }
                else
                {
                    if (!SILENT) cout << "out";
                }
                if (!SILENT) cout << " renumbering serials and with";
                if (renumb_res >= 0)
                {
                    renum_res(ligand, 1);
                }
                else
                {
                    if (!SILENT) cout << "out";
                }
                if (!SILENT) cout << " renumbering residues\n";
                out_pdb(ligand, ligout.substr(0, k) + "_" + pdbstruct.chains[ch].letter + SUFF + postf, false);
            }
            buf[pdbstruct.chains[ch].letter].insert(ligname);
        }
    }

    if (del)
    {
        if (!SILENT) cout << "Deleting ligands from all chains ...\n";
        set <string> tmp;
        tmp = del_het(pdbstruct, buf);
    }
	return;
}

void cutting_prot(pdb_file &pdbstruct, string proout, int renumb_res, int renumb_ser)
{
	string postf = "";
	int k = proout.size();
	chain newchain;
	pdb_file to_chain;
	to_chain.chains.push_back(newchain);
	while (k > 0 && proout[k] != '.')
	{
		k--;
	}
	if (k == 0)
    {
        postf = ".pdb";
        k = proout.size();
    }
    else
    {
        postf = proout.substr(k, proout.size() - k);
    }

	for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
	{
	    to_chain.chains[0] = pdbstruct.chains[ch];
	    for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
            {
                for (int i = 0; i < pdbstruct.chains[ch].array[seg].atom[res].array.size(); ++i)
                {
                    if (pdbstruct.conect.count(pdbstruct.chains[ch].array[seg].atom[res].array[i].serial))
                    {
                        to_chain.conect[pdbstruct.chains[ch].array[seg].atom[res].array[i].serial] = pdbstruct.conect[pdbstruct.chains[ch].array[seg].atom[res].array[i].serial];
                    }
                }

            }
        }
        for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].hetatm.size(); ++res)
            {
                for (int i = 0; i < pdbstruct.chains[ch].array[seg].hetatm[res].array.size(); ++i)
                {
                    if (pdbstruct.conect.count(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].serial))
                    {
                        to_chain.conect[pdbstruct.chains[ch].array[seg].hetatm[res].array[i].serial] = pdbstruct.conect[pdbstruct.chains[ch].array[seg].hetatm[res].array[i].serial];
                    }
                }

            }
        }
        if (renumb_ser > 0)
        {
            renum_ser(to_chain, renumb_ser, true);
        }
        if (renumb_res >= 0)
        {
            renum_res(to_chain, renumb_res);
        }
		out_pdb(to_chain, proout.substr(0, k) + "_" + pdbstruct.chains[ch].letter + SUFF + postf, true);
	}
	return;
}

struct dock_clust
{
    int rmsd;//1.000, A
    int ki;//nanomolar
    int dg;//-8.00 kcal/mol
    int number;//structures in cluster
    residue res;
};
struct dlg_file
{
    vector<dock_clust> clust;

    string ligand_types;//DPF> ligand_types A C F OA
    string node;//Successful Completion on "----"
    string time;//Real= -----
    int ga_num_generations;//DPF> ga_num_generations 27000
    int ga_num_evals;//DPF> ga_num_generations 27000
    int ga_pop_size;//DPF> ga_pop_size 150
    int runs_done;//DPF> ga_pop_size 150
    int torsdof;//DPF> torsdof 3
    int ga_run;//DPF> ga_run 50
    int x, y, z;//DPF> about 26.6874 33.4806 198.7594
};


struct vina_log{
    vector <int> affin;
    vector <int> rmsdlb;
    vector <int> rmsdub;
};

dlg_file parse_dlg_file(string dlgfile)
{
    string scurf = "parse_dlg_file", smsg = "";
    dlg_file file;
    file.ga_run = 0;
    if (!file_exists(dlgfile))
    {
        return file;
    }
    ifstream din(dlgfile.c_str());
    string s, ss;
    bool cluster = false, dpf = false, endit = false;
    int result = 0;
    int details = -1, a, str = 1;
    dock_clust buf;
    atom batom;
    while (getline(din, s))
    {
        if (s.find("LOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER") != string::npos)
        {
            cluster = true;
            continue;
        }
        else if (s.find("SETTING UP DEFAULT PARAMETER LIBRARY") != string::npos)
        {
            dpf = true;
            continue;
        }
        else if (s.find("BEGINNING GENETIC ALGORITHM DOCKING") != string::npos)
        {
            dpf = false;
            continue;
        }
        else if (s.find("FINAL GENETIC ALGORITHM ALGORITHM DOCKED STATE") != string::npos)
        {
            result++;
            continue;
        }
        if (cluster)
        {
            ss = s.substr(0, 6);
            if (ss == "USER  ")
            {
                if (s.find("RMSD from reference structure") != string::npos)
                {
                    a = s.find("=") + 2;
                    file.clust[file.clust.size() - 1].rmsd = strtodoub(next_lex(s, a, ' '));
                }
                else if (s.find("Number of conformations in this cluster") != string::npos)
                {
                    a = s.find("=") + 2;
                    file.clust[file.clust.size() - 1].number = strtoint(next_lex(s, a, ' '));
                }
                else if (s.find("Estimated Free Energy of Binding") != string::npos)
                {
                    a = s.find("=") + 2;
                    file.clust[file.clust.size() - 1].dg = strtodoub(next_lex(s, a, ' '));
                }
                else if (s.find("Estimated Inhibition Constant, Ki") != string::npos)
                {
                    a = s.find("=") + 2;
                    file.clust[file.clust.size() - 1].ki = strtodoub(next_lex(s, a, ' '));
                    a += 8;
                    ss = next_lex(s, a, ' ');
                    if (ss[0] == 'u')
                    {
                        file.clust[file.clust.size() - 1].ki *= 1000;
                    }
                    if (ss[0] == 'm')
                    {
                        file.clust[file.clust.size() - 1].ki *= 1000000;
                    }
                }
            }
            else if (ss == "ATOM  ")
            {
                batom.serial = strtoint(s.substr(6, 5));//Serial no
                batom.atom_name = trimdel(s.substr(12, 4));//Название атома
                batom.alt = s[16];//Альтернативное положение
                batom.resname = trimdel(s.substr(17, 4));//Название резидью
                batom.resnumb = trimdel(s.substr(22, 5));//Номер резидью;
				if (batom.resnumb == "")
				{
					batom.resnumb = "1";
				}
                batom.chain = s[21];//Название цепи

                batom.x = strtodoub(s.substr(30, 8));//Координата x
                batom.y = strtodoub(s.substr(38, 8));//Координата y
                batom.z = strtodoub(s.substr(46, 8));//Координата z
                batom.oq = 0;//Occupancy
                batom.tf = 0;//Температурный фактор
                batom.ch_val = 0;
                batom.ch_sign = -1;
                batom.element = del_digits(batom.atom_name);
				/*if (to_upper(batom.element) == "FE" || to_upper(batom.element) == "CO" || to_upper(batom.element) == "MN")
				{
					smsg = "There is atom with name ";
					smsg += batom.element;
					smsg += ". It seems to be a metal";
					warning(scurf, smsg);
				}
                else
				{

                }*/
                batom.element = batom.element.substr(0, 1);
				file.clust[file.clust.size() - 1].res.array.push_back(batom);
                file.clust[file.clust.size() - 1].res.number = batom.resnumb;
            }
            else if (ss == "MODEL ")
            {
                file.clust.push_back(buf);
                details = 1;
            }
            else if (ss == "ENDMDL")
            {
                details = 0;
            }
            else if (ss == "TER")
            {
                details = 0;
            }
            else
            {
                if (details == 1)
                {
                    smsg = "There is an error in dlg file in string ";
                    smsg += inttostr_2(str);
                    warning(scurf, smsg);
                }
                else if (details == 0)
                {
                    cluster = false;
                    endit = true;
                }
            }
        }
        else if (dpf)
        {
            if (s.substr(0, 4) == "DPF>")
            {
                a = 4;
                ss = next_lex(s, a, ' ');
                if (ss == "ligand_types")
                {
                    file.ligand_types = trimdel(s.substr(18, s.find("#") - 18));
                }
                else if (ss == "ga_num_generations")
                {
                    ss = next_lex(s, a, ' ');
                    file.ga_num_generations = strtoint(ss);
                }
                else if (ss == "ga_pop_size")
                {
                    ss = next_lex(s, a, ' ');
                    file.ga_pop_size = strtoint(ss);
                }
                else if (ss == "torsdof")
                {
                    ss = next_lex(s, a, ' ');
                    file.torsdof = strtoint(ss);
                }
                else if (ss == "ga_run")
                {
                    ss = next_lex(s, a, ' ');
                    file.ga_run = strtoint(ss);
                }
                else if (ss == "about")
                {
                    file.x = strtodoub(next_lex(s, a, ' '));
                    file.y = strtodoub(next_lex(s, a, ' '));
                    file.z = strtodoub(next_lex(s, a, ' '));
                }
                else if (ss == "ga_num_evals")
                {
                    ss = next_lex(s, a, ' ');
                    file.ga_num_evals = strtoint(ss);
                }
            }
        }
        else if (endit)
        {
            if (s.find("Successful Completion on") != string::npos)
            {
                a = s.find("\"") + 1;
                file.node = next_lex(s, a, '\"');
            }
            else if (s.find("Real= ") != string::npos)
            {
                a = s.find("=") + 2;
                file.time = next_lex(s, a, ',');
            }
        }
        ++str;
	}
	if (result < file.ga_run)
	{
	    smsg = "There are less results in dlg file than need ";
	    smsg += inttostr_2(result);
	    smsg += " < ";
	    smsg += inttostr_2(file.ga_run);
	    smsg += ".\nNo results recieved ";
	    smsg += dlgfile;
        warning(scurf, smsg);
        file.runs_done = result;
    }
	din.close();
	return file;
}

map <string, string> get_rec_info(string infofile)
{
    string scurf = "get_rec_info", smsg = "";
    if (!file_exists(infofile))
    {
        smsg = "There is no file ";
        smsg += infofile;
        error_msg(scurf, smsg);
        exit(1);
    }
    ifstream infin(infofile.c_str());
    string s, pdbcode;
    int k;
    map <string, string> desc;
    while (getline(infin, s))
    {
        if (s[0] == '#')
        {
            continue;
        }
        k = s.find('\t');
        if (k != string::npos)
        {
            pdbcode = s.substr(0, k);
            if (desc.count(pdbcode))
            {
                smsg = "Double pdbcode ";
                smsg += pdbcode;
                smsg += " declaration in receptor info file. Second entry was skipped";
                warning(scurf, smsg);
            }
            s = s.substr(k + 1, s.size() - k - 1);
        }
        else
        {
            smsg = "Wrong receptor info file format. String was skipped";
            warning(scurf, smsg);
            continue;
        }
        k = s.find('\t');
        if (k != string::npos)
        {
            desc[pdbcode] = s.substr(0, k);
        }
        else
        {
            smsg = "Wrong receptor info file format 2. String was skipped";
            warning(scurf, smsg);
        }
    }

    infin.close();
    return(desc);
}

bool compf(pair < string, dlg_file > item1, pair < string, dlg_file > item2)
{
    return item1.second.clust[0].dg < item2.second.clust[0].dg;
}

bool compf_vina(pair < string, vina_log > item1, pair < string, vina_log > item2)
{
    return item1.second.affin[0] < item2.second.affin[0];
}

string strict(string s, string prog)
{
    string scurf = "strict", smsg;
    vector <string> splited = split(s, " \t");
    string sch, path, name, exp;
    int n;
    if (prog == "auto")
    {
        sch = "-l";
    }
    else if (prog == "vina")
    {
        sch = "--log";
    }
    for (int i = 0; i < splited.size(); ++i)
    {
        if (splited[i] == sch)
        {
            //cout << splited[i] << "*\n";
            if (i + 1 < splited.size())
            {
                filename_parse(splited[i + 1], path, name, exp);
                s = string(splited[0] + name + "." + exp);
                return s;
            }
            else
            {
                smsg = "There is no output file in runlist\n";
                warning(scurf, smsg);
                return s;
            }
        }
    }
    smsg = "There is no output parameter in runlist\n";
    warning(scurf, smsg);
    return s;
}

void make_redock_runlist(string dlg_list_file, vector <string> redock)
{
    string path, name, ext;
    int n = -1;
    filename_parse(dlg_list_file, path, name, ext);
    if (name.find("_redock") == -1)
    {
        name += "_redock";
    }
    dlg_list_file = get_file_name(path, name, ext);
    ofstream rout(dlg_list_file.c_str());
    for (int i = 0; i < redock.size(); ++i)
    {
        rout << redock[i] << "\n";
    }
    return;
}

vector <string> dlg_parser(string file_list, string file_out, bool roc_plot = false, bool runlist = false)
{
    string scurf = "dlg_parser", smsg = "";
    if (!file_exists(file_list))
    {
        smsg = "There is no file ";
        smsg += file_list;
        error_msg(scurf, smsg);
        exit(1);
    }
    ifstream pin(file_list.c_str());
    string path, sname, ext;
    filename_parse(file_out, path, sname, ext);
    file_out = get_file_name(path, sname, ext);
    ofstream pout(file_out.c_str());
    map <string, string> descript;
    dlg_file file;
    double square;
    vector <pair < string, dlg_file > > codefiles;
    string name, pdbcode, curcode = "";
    vector <string> redock;
    char curchain, ch;
    int k, p = 1, z, cc, x = 0, y = 0, amnot, amyes, sc;
    bool info = true, rif = false;
    getline(pin, name);
    pin.sync();
    //cout << name;
    if (name[0] == '#')
    {
        redock.push_back(name);
        k = 0;
        for (int i = 0; i < name.size(); ++i)
        {
            if (name[i] == '/')
            {
                k = i;
            }
        }
        if (!SILENT) cout << "Receptors info file: " << name.substr(k + 1, name.size() - k - 1) << "\n";
        descript = get_rec_info(name.substr(1, name.size() - 1));
        rif = true;
        getline(pin, name);
    }
    curchain = '-';
    curcode = "";
    sc = 0;
    do
    {
        pin.sync();
        if (runlist)
        {
            redock.push_back(name);
            name = strict(name, "auto");
        }
        if (sc == 500)
        {
            if (!SILENT) cout << "Processing " << name << "\n";
            sc = 1;
        }
        sc++;

        file = parse_dlg_file(name);
        filename_parse(name, path, name, ext);
        if (file.clust.size() > 0)
        {
            k = name.find('_');
            pdbcode = name.substr(0, k);
            name = name.substr(k + 1, name.size() - k - 1);
            z = 0;
            k = name.find('_');
            cc = -1;
            p = 0;
            while (k != string::npos)
            {
                if (k - z == 1)
                {
                    if (isupletter(name[z]))
                    {
                        if (cc != -1)
                        {
                            smsg = "Double 1-letter parameter in file name. Using 1st as chain";
                            warning(scurf, smsg);
                        }
                        else
                        {
                            ch = name[z];
                            cc = p;
                        }
                    }
                }
                z = k + 1;
                name.replace(k, 1, "\t");
                k = name.find('_');
                p++;
            }
            if (name.size() - z == 1)
            {
                if (isupletter(name[z]))
                {
                    if (cc != -1)
                    {
                        smsg = "Double 1-letter parameter in file name. Using 1st as chain";
                        warning(scurf, smsg);
                    }
                    else
                    {
                        ch = name[z];
                        cc = p;
                    }
                }
            }
            if (cc == -1)
            {
                smsg = "There is no 1-big-letter parameter in file name. Using chain A";
                warning(scurf, smsg);
                ch = 'A';
            }
            if (roc_plot)
            {
                //cout << curcode << "*" << pdbcode << "*" << curchain << "*" << ch << "*\n";
                if (curcode == "" || curchain == '-')
                {
                    curcode = pdbcode;
                    curchain = ch;
                    codefiles.push_back(make_pair(name, file));
                    continue;
                }
                if (curchain != ch || curcode != pdbcode)
                {
                    sort(codefiles.begin(), codefiles.end(), compf);
                    amnot = 0;
                    amyes = 0;
                    for (int i = 0; i < codefiles.size(); ++i)
                    {
                        if (codefiles[i].first.find("NOT") != string::npos || codefiles[i].first.find("not") != string::npos)
                        {
                            amnot++;
                        }
                        else
                        {
                            amyes++;
                        }
                    }
                    //cout << amnot << "&" << amyes << "&\n";
                    x = 0;
                    y = 0;
                    square = 0.;
                }
                else
                {
                    codefiles.push_back(make_pair(name, file));
                    continue;
                }
            }
            else
            {
                codefiles.clear();
                codefiles.push_back(make_pair(name, file));

            }
            for (int i = 0; i < codefiles.size(); ++i)
            {
                if (p > 0 && info)
                {
                    info = false;
                    pout << "LIGFILE\t";
                    pout << "PDBCode\t";
                    if (rif)
                    {
                        pout << "Description\t";
                    }
                    for (int j = 0; j <= p; ++j)
                    {
                        if (j == cc)
                        {
                            pout << "Chain\t";
                        }
                        else
                        {
                            pout << "par" << j + 1 << "\t";
                        }

                    }
                    pout << "ga_num_generations\t";
                    pout << "ga_num_evals\t";
                    pout << "ga_pop_size\t";
                    pout << "ga_run\t";
                    pout << "Amount of active torsions\t";
                    pout << "Calculation time\t";
                    pout << "Number of clusters\t";
                    pout << "Amount in 1st\t";
                    pout << "1st best pose rmsd from reference\t";
                    pout << "1st best delta G\t";
                    pout << "Amount in 2nd\t";
                    pout << "2nd best pose rmsd from reference\t";
                    pout << "2nd best delta G";
                    if (roc_plot)
                    {
                        pout << "\t";
                        pout << "xnplot\tynplot\tSquare";
                    }
                    pout << "\n";
                }
                pout << redock.back() << "\t";
                if (roc_plot)
                {
                    if (codefiles[i].first.find("NOT") != string::npos || codefiles[i].first.find("not") != string::npos)
                    {
                        x++;
                        square += 1. / amnot * ((amyes - y) * 1. / amyes);
                    }
                    else
                    {
                        y++;
                    }
                    pout << curcode << "\t";
                }
                else
                {
                    pout << pdbcode << "\t";
                }
                if (rif)
                {
                    if (descript.count(pdbcode))
                    {
                        pout << descript[pdbcode] << "\t";
                    }
                    else
                    {
                        pout << "none\t";
                    }
                }
                pout << codefiles[i].first << "\t";
                pout << codefiles[i].second.ga_num_generations << "\t";
                pout << codefiles[i].second.ga_num_evals << "\t";
                pout << codefiles[i].second.ga_pop_size << "\t";
                pout << codefiles[i].second.ga_run << "\t";
                pout << codefiles[i].second.torsdof << "\t";
                pout << codefiles[i].second.time << "\t";
                pout << codefiles[i].second.clust.size() << "\t" << codefiles[i].second.clust[0].number << "\t" << doubtostr_2(codefiles[i].second.clust[0].rmsd, 3) << "\t" <<  doubtostr_2(codefiles[i].second.clust[0].dg, 2);
                if (codefiles[i].second.clust.size() > 1)
                {
                    pout << "\t" << codefiles[i].second.clust[1].number << "\t" << doubtostr_2(codefiles[i].second.clust[1].rmsd, 3) << "\t" <<  doubtostr_2(codefiles[i].second.clust[1].dg, 2);
                }
                else
                {
                    pout << "\t-\t-\t-";
                }
                if (roc_plot)
                {
                    pout.imbue(locale(cout.getloc(), new DecimalSeparator<char>(',')));
                    pout << "\t" << x * 1. / amnot << "\t" << y * 1. / amyes;
                    if (i == codefiles.size() - 1)
                    {
                        pout << "\t" << 1 - square;
                    }
                    //cout << x * 1. / amnot << "\t" << y * 1. / amyes << "\n";
                }
                pout << "\n";
            }
            if (roc_plot)
            {
                codefiles.clear();
                codefiles.push_back(make_pair(name, file));
                curchain = ch;
                curcode = pdbcode;
            }
            redock.pop_back();
        }
        else
        {
            if (!NO_DLG_OUT_ERRORS)
            {
                pout << "File with no results\t" << file.runs_done << "\n";
            }
            smsg = name;
            smsg += " - file with no clusters";
            warning(scurf, smsg);
        }
    }
    while (getline(pin, name));
    if (roc_plot)
    {
        sort(codefiles.begin(), codefiles.end(), compf);
        amnot = 0;
        amyes = 0;
        for (int i = 0; i < codefiles.size(); ++i)
        {
            if (codefiles[i].first.find("ingNOT") != string::npos)
            {
                amnot++;
            }
            else
            {
                amyes++;
            }
        }
        x = 0;
        y = 0;
        square = 0.;
        for (int i = 0; i < codefiles.size(); ++i)
        {
            if (codefiles[i].first.find("ingNOT") != string::npos)
            {
                x++;
                square += 1. / amnot * ((amyes - y) * 1. / amyes);
            }
            else
            {
                y++;
            }
            pout << curcode << "\t";
            if (rif)
            {
                if (descript.count(pdbcode))
                {
                    pout << descript[pdbcode] << "\t";
                }
                else
                {
                    pout << "none\t";
                }
            }
            pout << codefiles[i].first << "\t";
            pout << codefiles[i].second.ga_num_generations << "\t";
            pout << codefiles[i].second.ga_num_evals << "\t";
            pout << codefiles[i].second.ga_pop_size << "\t";
            pout << codefiles[i].second.ga_run << "\t";
            pout << codefiles[i].second.torsdof << "\t";
            pout << codefiles[i].second.time << "\t";
            pout << codefiles[i].second.clust.size() << "\t" << codefiles[i].second.clust[0].number << "\t" << doubtostr_2(codefiles[i].second.clust[0].rmsd, 3) << "\t" <<  doubtostr_2(codefiles[i].second.clust[0].dg, 2);
            if (codefiles[i].second.clust.size() > 1)
            {
                pout << "\t" << codefiles[i].second.clust[1].number << "\t" << doubtostr_2(codefiles[i].second.clust[1].rmsd, 3) << "\t" <<  doubtostr_2(codefiles[i].second.clust[1].dg, 2);
            }
            else
            {
                pout << "\t-\t-\t-";
            }
            pout.imbue(locale(cout.getloc(), new DecimalSeparator<char>(',')));
            pout << "\t" << x * 1. / amnot << "\t" << y * 1. / amyes;
            if (i == codefiles.size() - 1)
            {
                pout << "\t" << 1 - square;
            }
            //cout << x * 1. / amnot << "\t" << y * 1. / amyes << "\n";
            pout << "\n";
        }
    }
    pin.close();
    pout.close();
    return redock;
}

vina_log parse_vina_log_file(string log_file)
{
    string scurf = "parse_vina_log_file", smsg = "";
    string s;
    bool idata = false;
    vector <string> buf;
    vina_log res;
    if (!file_exists(log_file))
    {
        return res;
    }
    ifstream din(log_file.c_str());
    while (getline(din, s))
    {
        buf = split(s, " ");
        if (buf.size() > 0)
        {
            if (buf[0] == "mode")
            {
                idata = true;
            }
            else
            {
                if (idata)
                {
                    if (isdigit(buf[0][0]))
                    {
                        if (buf.size() == 4)
                        {
                            res.affin.push_back(strtodoub(buf[1]));
                            res.rmsdlb.push_back(strtodoub(buf[2]));
                            res.rmsdub.push_back(strtodoub(buf[3]));
                        }
                        else
                        {
                            smsg = "Wrong VINA log file format. Skipped ";
                            warning(scurf, smsg);
                        }
                    }
                }
            }
        }
    }
    din.close();
    return res;
}

vector <string> vina_log_parser(string file_list, string file_out, bool roc_plot = false, bool runlist = false)
{

    string scurf = "vina_log_parser", smsg = "";
    if (!file_exists(file_list))
    {
        smsg = "There is no file ";
        smsg += file_list;
        error_msg(scurf, smsg);
        exit(1);
    }
    ifstream pin(file_list.c_str());
    string path, basename, ext;
    filename_parse(file_out, path, basename, ext);
    file_out = get_file_name(path, basename, ext);
    ofstream pout(file_out.c_str());
    map <string, string> descript;
    vina_log file;
    double square;
    vector <string> namecode;
    vector <pair < string, vina_log > > codefiles;
    string name, pdbcode = "", curcode = "", dir = "";
    vector <string> redock;
    char curchain, ch;
    int p = 1, z, cc, x = 0, y = 0, amnot, amyes;
    bool info = true, rif = false;
    getline(pin, name);
    pin.sync();
    if (name[0] == '#')
    {
        redock.push_back(name);
        name = name.substr(1, name.size() - 1);
        filename_parse(name, path, basename, ext);
        if (!SILENT) cout << "Receptors info file: " << basename << "\n";
        descript = get_rec_info(name);
        rif = true;
        getline(pin, name);
    }
    curchain = '-';
    curcode = "";

    do
    {
        pin.sync();
        if (runlist)
        {
            redock.push_back(name);
            name = strict(name, "vina");
        }
        filename_parse(name, path, basename, ext);
        //if (!SILENT) cout << "Processing " << basename << "\n";
        file = parse_vina_log_file(name);
        if (file.affin.size() > 0)
        {
            dir = "";
            if (path.size() > 0)
            {

                for (int i = path.size() - 1; i >= 0; --i)
                {
                    if (path[i] == '/' && path.size() - i > 1)
                    {
                        dir = path.substr(i + 1, path.size() - i - 2);
                        break;
                    }
                }
                if (dir == "")
                {
                    dir = path;
                }
            }
            if (dir == "")
            {
                smsg = "Path to the log file ";
                smsg += path;
                smsg += " must contain at list receptors directory";
                warning(scurf, smsg);
                namecode.clear();
                pdbcode = "UNCR";
                ch = 'A';
            }
            else
            {
                //cout << dir << "\n";
                namecode = split(dir, "_/");
                pdbcode = "";
                ch = '-';
                for (int i = 0; i < namecode.size(); ++i)
                {
                    if (namecode[i] == "vina" || namecode[i] == "VINA" || namecode[i] == "Vina" || namecode[i] == "out" || namecode[i] == "OUT" || namecode[i] == "ed")
                    {
                        continue;
                    }
                    if (pdbcode == "")
                    {
                        if (namecode[i].size() == 4 && isuplettersword(namecode[i]))
                        {
                            pdbcode = namecode[i];
                        }
                    }
                    if (ch == '-')
                    {
                        if (namecode[i].size() == 1 && isuplettersword(namecode[i]))
                        {
                            ch = namecode[i][0];
                        }
                    }
                }
                if (pdbcode == "")
                {
                    smsg = "There is no receptor name for file ";
                    smsg += name;
                    smsg += ". Skipped";
                    warning(scurf, smsg);
                    pdbcode = "UNCR";
                }
                if (ch == '-')
                {
                    smsg = "There is no chain big letter (_C_) in file name. Using chain A";
                    warning(scurf, smsg);
                    ch = 'A';
                }
            }

            if (roc_plot)
            {
                //cout << curcode << "*" << pdbcode << "*" << curchain << "*" << ch << "*\n";
                if (curcode == "" || curchain == '-')
                {
                    curcode = pdbcode;
                    curchain = ch;
                    codefiles.push_back(make_pair(basename, file));
                    continue;
                }
                if (curchain != ch || curcode != pdbcode)
                {
                    sort(codefiles.begin(), codefiles.end(), compf_vina);
                    amnot = 0;
                    amyes = 0;
                    for (int i = 0; i < codefiles.size(); ++i)
                    {
                        if (codefiles[i].first.find("NOT") != string::npos || codefiles[i].first.find("not") != string::npos)
                        {
                            amnot++;
                        }
                        else
                        {
                            amyes++;
                        }
                    }
                    //cout << amnot << "&" << amyes << "&\n";
                    x = 0;
                    y = 0;
                    square = 0.;
                }
                else
                {
                    codefiles.push_back(make_pair(basename, file));
                    continue;
                }
            }
            else
            {
                codefiles.clear();
                codefiles.push_back(make_pair(basename, file));

            }
            for (int i = 0; i < codefiles.size(); ++i)
            {
                if (p > 0 && info)
                {
                    info = false;
                    pout << "LIGFILE\t";
                    pout << "PDBCode\t";
                    pout << "Chain\t";
                    if (rif)
                    {
                        pout << "Descr\t";
                    }

                    pout << "Compound\t";

                    pout << "1 delta G\t";
                    pout << "2 delta G\t";
                    pout << "2 rmsd lb\t";
                    pout << "2 rmsd ub";
                    if (roc_plot)
                    {
                        pout << "\txnplot\tynplot\tSquare";
                    }
                    pout << "\n";
                }
                pout << redock.back() << "\t";
                if (roc_plot)
                {
                    if (codefiles[i].first.find("NOT") != string::npos || codefiles[i].first.find("not") != string::npos)
                    {
                        x++;
                        if  (amnot > 0)
                        {
                            square += 1. / amnot * ((amyes - y) * 1. / amyes);
                        }
                    }
                    else
                    {
                        y++;
                    }
                    pout << curcode << "\t";
                    pout << curchain << "\t";
                }
                else
                {
                    pout << pdbcode << "\t";
                    pout << ch << "\t";
                }
                if (rif)
                {
                    if (descript.count(pdbcode))
                    {
                        pout << descript[pdbcode] << "\t";
                    }
                    else
                    {
                        pout << "none\t";
                    }
                }
                pout << codefiles[i].first << "\t";
                pout << doubtostr_2(codefiles[i].second.affin[0], 1) << "0\t";
                if (codefiles[i].second.affin.size() > 1)
                {
                    pout << doubtostr_2(codefiles[i].second.affin[1], 1) << "0\t";
                    pout << doubtostr_2(codefiles[i].second.rmsdlb[1], 3) << "\t";
                    pout << doubtostr_2(codefiles[i].second.rmsdub[1], 3) << "\t";
                }
                else
                {
                    pout << "\t-\t-\t-";
                }
                if (roc_plot)
                {
                    pout.imbue(locale(cout.getloc(), new DecimalSeparator<char>(',')));
                    pout << "\t" << x * 1. / (amnot > 0 ? amnot : 1) << "\t" << y * 1. / amyes;
                    if (i == codefiles.size() - 1)
                    {
                        pout << "\t" << 1 - square;
                    }
                    //cout << x * 1. / amnot << "\t" << y * 1. / amyes << "\n";
                }
                pout << "\n";
            }
            if (roc_plot)
            {
                codefiles.clear();
                codefiles.push_back(make_pair(basename, file));
                curchain = ch;
                curcode = pdbcode;
            }
            redock.pop_back();
        }
        else
        {
            if (!NO_DLG_OUT_ERRORS)
            {
                pout << "File with no results\t" << name << "\n";
            }
            smsg = name;
            smsg += " - file with no results";
            warning(scurf, smsg);

        }
    }
    while (getline(pin, name));
    pin.close();
    pout.close();
    return redock;
}

void get_lig_dlg(string input_dlg, string out_file, int clust = 1)
{
    dlg_file file = parse_dlg_file(input_dlg);
    pdb_file ligand;
    chain newchain;
    segment newsegment;
    ligand.chains.push_back(newchain);
    ligand.chains[0].array.push_back(newsegment);
    ligand.chains[0].letter = 'A';
    ligand.chains[0].array[0].cname = "";
    int k = out_file.size();
    bool cl = true;
    if (clust == 1)
    {
        cl = false;
    }
    string postf = ".pdb";
    for (int i = out_file.size() - 1; i > 0; --i)
    {
        if (out_file[i] == '.')
        {
            k = i;
            break;
        }
    }
    if (file.clust.size() > 0)
    {
        for (int i = 0; i < min(clust, file.clust.size()); ++i)
        {
            ligand.chains[0].array[0].hetatm.clear();
            ligand.chains[0].array[0].hetatm.push_back(file.clust[i].res);
            out_pdb(ligand, out_file.substr(0, k) + (cl ? "_cl" + inttostr_2(i) : "") + postf, 0);
        }
    }
    return;
}

double distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

double lig_correcter(string sample, string coord, string output)
{
    string scurf = "lig_correcter", smsg = "", postf;
    pdb_file s_lig, c_lig;
    parse_pdb(sample, s_lig);
    parse_pdb(coord, c_lig);
    residue rs, rc;
    double loc_rmsd = 0, glob_rmsd = 0;
    bool reams = false, reamc = false;
    int a = 0, k = 0;
    for (int ch = 0; ch < s_lig.chains.size(); ++ch)
    {
        for (int seg = 0; seg < s_lig.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < s_lig.chains[ch].array[seg].atom.size(); ++res)
            {
                if (!reams)
                {
                    rs = s_lig.chains[ch].array[seg].atom[res];
                    reams = true;
                }
                else
                {
                    smsg = "More than 1 ligand in sample file, other function was run";
                    warning(scurf, smsg);
                    return lig_correcter_dif(sample, coord, output);
                }
            }
            for (int res = 0; res < s_lig.chains[ch].array[seg].hetatm.size(); ++res)
            {
                if (!reams)
                {
                    rs = s_lig.chains[ch].array[seg].hetatm[res];
                    reams = true;
                }
                else
                {
                    smsg = "More than 1 ligand in sample file, other function was run";
                    warning(scurf, smsg);
                    return lig_correcter_dif(sample, coord, output);
                }
            }
        }
    }
    for (int ch = 0; ch < c_lig.chains.size(); ++ch)
    {
        for (int seg = 0; seg < c_lig.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < c_lig.chains[ch].array[seg].atom.size(); ++res)
            {
                if (!reamc)
                {
                    rc = c_lig.chains[ch].array[seg].atom[res];
                    reamc = true;
                }
                else
                {
                    smsg = "More than 1 ligand in file with coordinates, other function was run";
                    warning(scurf, smsg);
                    return lig_correcter_dif(sample, coord, output);
                }
            }
            for (int res = 0; res < c_lig.chains[ch].array[seg].hetatm.size(); ++res)
            {
                if (!reamc)
                {
                    rc = c_lig.chains[ch].array[seg].hetatm[res];
                    reamc = true;
                }
                else
                {
                    smsg = "More than 1 ligand in file with coordinates, other function was run";
                    warning(scurf, smsg);
                    return lig_correcter_dif(sample, coord, output);
                }
            }
        }
    }
    if (rs.array.size() != rc.array.size())
    {
        smsg = "Sample ligand and ligand for correcting have different amounts of atoms, other function was run";
        warning(scurf, smsg);
        return lig_correcter_dif(sample, coord, output);
    }
    for (int i = 0; i < rs.array.size(); ++i)
    {
        if (rs.array[i].serial != rc.array[i].serial)
        {
            smsg = "Sample ligand and ligand for correcting have different serial numbers ";
            smsg += inttostr_2(rs.array[i].serial);
            smsg += " and ";
            smsg += inttostr_2(rc.array[i].serial);
            warning(scurf, smsg);
        }
        loc_rmsd = (distance(rs.array[i].x, rs.array[i].y, rs.array[i].z, rc.array[i].x, rc.array[i].y, rc.array[i].z) * 1.0) / 1000;
        if (loc_rmsd > MAXWARNLIGCOR)
        {
            smsg = "Single atom deviation ";
            smsg += doubtostr_3(loc_rmsd, 3);
            smsg += " more than max value. Atom serial ";
            smsg += inttostr_2(rs.array[i].serial);
            smsg += ". File ";
            for (int k = 0; k < sample.size() - 1; ++k)
            {
                if (sample[k] == '/')
                {
                    a = k;
                }
            }
            smsg += sample.substr(a + 1, sample.size() - a - 1);
            warning(scurf, smsg);
        }
        glob_rmsd += loc_rmsd;
        rs.array[i].x = rc.array[i].x;
        rs.array[i].y = rc.array[i].y;
        rs.array[i].z = rc.array[i].z;
    }
    reams = false;
    for (int ch = 0; ch < s_lig.chains.size(); ++ch)
    {
        for (int seg = 0; seg < s_lig.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < s_lig.chains[ch].array[seg].atom.size(); ++res)
            {
                if (!reams)
                {
                    for (int i = 0; i < s_lig.chains[ch].array[seg].atom[res].array.size(); ++i)
                    {
                        s_lig.chains[ch].array[seg].atom[res].array[i].x = rs.array[i].x;
                        s_lig.chains[ch].array[seg].atom[res].array[i].y = rs.array[i].y;
                        s_lig.chains[ch].array[seg].atom[res].array[i].z = rs.array[i].z;
                    }
                    reams = true;
                }
            }
            for (int res = 0; res < s_lig.chains[ch].array[seg].hetatm.size(); ++res)
            {
                if (!reams)
                {
                    for (int i = 0; i < s_lig.chains[ch].array[seg].hetatm[res].array.size(); ++i)
                    {
                        s_lig.chains[ch].array[seg].hetatm[res].array[i].x = rs.array[i].x;
                        s_lig.chains[ch].array[seg].hetatm[res].array[i].y = rs.array[i].y;
                        s_lig.chains[ch].array[seg].hetatm[res].array[i].z = rs.array[i].z;
                    }
                    reams = true;
                }
            }
        }
    }
    s_lig.conect = c_lig.conect;
    k = output.size() - 1;
    postf = "";
    while (k > 0 && output[k] != '.')
    {
        k--;
    }
    if (k == 0)
    {
        postf = ".pdb";
        k = output.size();
    }
    else
    {
        postf = output.substr(k, output.size() - k);
    }
    out_pdb(s_lig, output.substr(0, k) + SUFF + postf, 0);
    return glob_rmsd;
}

vector <string> get_coord_strings(string file)
{
    string scurf = "get_coord_strings", smsg = "";
    string s, ss;
    vector <string> ret;
    ret.clear();
    if (!file_exists(file))
    {
        smsg = "There is no file ";
        smsg += file;
        error_msg(scurf, smsg);
        exit(1);
    }
    ifstream fin (file.c_str());
    while (getline(fin, s))
    {
        ss = s.substr(0, 6);
        if (ss == "ATOM  " || ss == "HETATM" || ss == "CONECT")
        {
            ret.push_back(s);
        }
    }
    fin.close();
    return ret;
}

void out_lines(vector<string> s, string file)
{
    string path, name, ext;
    filename_parse(file, path, name, ext);
    file = get_file_name(path, name, ext);
    ofstream fout(file.c_str());
    for (int i = 0; i < s.size(); ++i)
    {
        fout << s[i] << "\n";
    }
    out_end(fout);
    fout.close();
    return;
}

double lig_correcter_dif(string sample, string coord, string output)
{
    string scurf = "lig_correcter_dif", smsg = "", postf;
    vector <string> s, c;
    int satoms = 0, catoms = 0, a = 0, k = 0;
    double loc_rmsd = 0, glob_rmsd = 0;
    atom sa, ca;
    s = get_coord_strings(sample);
    c = get_coord_strings(coord);
    for (int i = 0; i < s.size(); ++i)
    {
        if (s[i].substr(0, 6) == "ATOM  " || s[i].substr(0, 6) == "HETATM")
        {
            ++satoms;
        }
    }
    for (int i = 0; i < c.size(); ++i)
    {
        if (c[i].substr(0, 6) == "ATOM  " || c[i].substr(0, 6) == "HETATM")
        {
            ++catoms;
        }
    }
    if (satoms != catoms)
    {
        smsg = "Different amount of atoms ";
        smsg += inttostr_2(s.size());
        smsg += " coords in 1st ";
        smsg += inttostr_2(c.size());
        smsg += " coords in 2nd ";
        smsg += inttostr_2(satoms);
        smsg += " atom or hetatm in 1st ";
        smsg += inttostr_2(catoms);
        smsg += " atom or hetatm in 2nd";
        warning(scurf, smsg);
        return -1;
    }
    for (int i = 0; i < satoms; ++i)
    {
        if (strtoint(s[i].substr(6, 5)) != strtoint(c[i].substr(6, 5)))
        {
            smsg = "Different serial numbers ";
            smsg += inttostr_2(strtoint(s[i].substr(6, 5)));
            smsg += " and ";
            smsg += inttostr_2(strtoint(s[i].substr(6, 5)));
            warning(scurf, smsg);
        }
        sa.x = strtodoub(s[i].substr(30, 8));//Координата x
        sa.y = strtodoub(s[i].substr(38, 8));//Координата y
        sa.z = strtodoub(s[i].substr(46, 8));//Координата z
        ca.x = strtodoub(c[i].substr(30, 8));//Координата x
        ca.y = strtodoub(c[i].substr(38, 8));//Координата y
        ca.z = strtodoub(c[i].substr(46, 8));//Координата z
        loc_rmsd = (distance(sa.x, sa.y, sa.z, ca.x, ca.y, ca.z) * 1.0) / 1000;
        glob_rmsd += loc_rmsd;
        if (loc_rmsd > MAXWARNLIGCOR)
        {
            smsg = "Single atom deviation ";
            smsg += doubtostr_3(loc_rmsd, 3);
            smsg += " more than max value. Atom serial ";
            smsg += inttostr_2(strtoint(s[i].substr(6, 5)));
            smsg += ". File ";
            for (int k = 0; k < sample.size() - 1; ++k)
            {
                if (sample[k] == '/')
                {
                    a = i;
                }
            }
            smsg += sample.substr(a + 1, sample.size() - a);
            warning(scurf, smsg);
        }
        s[i] = s[i].substr(0, 30) + c[i].substr(30, 24) + s[i].substr(54, 80 - 54);
    }
    for (int i = 0; i < c.size(); ++i)
    {
        if (c[i].substr(0, 6) == "CONECT")
        {
            s.push_back(c[i]);
        }
    }
    k = output.size() - 1;
    postf = "";
    while (k > 0 && output[k] != '.')
    {
        k--;
    }
    if (k == 0)
    {
        postf = ".pdb";
        k = output.size();
    }
    else
    {
        postf = output.substr(k, output.size() - k);
    }
    out_lines(s, output.substr(0, k) + SUFF + postf);
    return glob_rmsd;
}

int make_cif_pdb(string ciffile, string lig, string outfile)
{
    string scurf = "make_cif_pdb", smsg = "";
    if (!file_exists(ciffile))
    {
        smsg = "There is no file ";
        smsg += ciffile;
        error_msg(scurf, smsg);
        exit(1);
    }
    ifstream fin (ciffile.c_str());
    string s, ss, curs1, curs2;
    bool coords = false, connect = false;
    residue ligres;
    atom resatom;
    int k = 0, at1, at2, amstereo = 1;
    resatom.resname = lig;
    resatom.resnumb = 1;
    resatom.alt = ' ';
    resatom.ch_sign = -1;
    resatom.ch_val = 0;
    resatom.oq = 0;
    resatom.tf = 0;
    ligres.number = "1";
    pdb_file lig_pdb;
    chain newchain;
    newchain.letter = 'A';
    lig_pdb.chains.push_back(newchain);
    segment newsegment;
    newsegment.cname = "";
    lig_pdb.chains[0].array.push_back(newsegment);
    vector<string> connects;
    while (getline(fin, s))
    {
        k = 0;
        ss = next_lex(s, k, ' ');
        if (ss == lig)
        {
            if (coords)
            {
                connects.push_back(next_lex(s, k, ' '));
                resatom.atom_name = next_lex(s, k, ' ');
                resatom.element = next_lex(s, k, ' ');
                for (int i = 0; i < 4; ++i)
                {
                    ss = next_lex(s, k, ' ');
                }
                ss = next_lex(s, k, ' ');
                if (ss != "N")
                {
                    amstereo *= 2;
                }
                resatom.x = strtodoub(next_lex(s, k, ' '));
                resatom.y = strtodoub(next_lex(s, k, ' '));
                resatom.z = strtodoub(next_lex(s, k, ' '));
                for (int i = 0; i < 5; ++i)
                {
                    ss = next_lex(s, k, ' ');
                }
                resatom.serial = strtoint(next_lex(s, k, ' '));
                ligres.array.push_back(resatom);
            }
            else if (connect)
            {
                curs1 = next_lex(s, k, ' ');
                curs2 = next_lex(s, k, ' ');
                ss = next_lex(s, k, ' ');
                for (int i = 0; i < connects.size(); ++i)
                {
                    if (curs1 == connects[i])
                    {
                        at1 = i + 1;
                        break;
                    }
                }
                for (int i = 0; i < connects.size(); ++i)
                {
                    if (curs2 == connects[i])
                    {
                        at2 = i + 1;
                        break;
                    }
                }
                lig_pdb.conect[at1].insert(at2);
                lig_pdb.conect[at2].insert(at1);
                if (ss == "DOUB")
                {
                    lig_pdb.conect[at1].insert(at2);
                    lig_pdb.conect[at2].insert(at1);
                }
                ss = next_lex(s, k, ' ');
                ss = next_lex(s, k, ' ');
                if (ss != "N")
                {
                    amstereo *= 2;
                }
            }
        }
        else if (ss == "#")
        {
            coords = false;
            connect = false;
        }
        else if (s.find("_chem_comp_atom") != string::npos)
        {
            coords = true;
        }
        else if (s.find("_chem_comp_bond") != string::npos)
        {
            connect = true;
        }
    }
    lig_pdb.chains[0].array[0].hetatm.push_back(ligres);
    out_pdb(lig_pdb, outfile, 0);
    return amstereo;
}

bool append_hetatm(pdb_file &pdbstruct, char chname, residue res)
{
    string scurf = "append_hetatm", smsg = "";
    bool p = true;
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        if (pdbstruct.chains[ch].letter == chname)
        {
            for (int i = 0; i < res.array.size(); ++i)
            {
                res.array[i].chain = pdbstruct.chains[ch].letter;
            }
            if (pdbstruct.chains[ch].array.size() == 0)
            {
                segment newsegment;
                pdbstruct.chains[ch].array.push_back(newsegment);
            }
            pdbstruct.chains[ch].array[pdbstruct.chains[ch].array.size() - 1].hetatm.push_back(res);
            p = false;
        }
    }
    if (p)
    {
        smsg = "There is no chain ";
        smsg += chname;
        smsg += ". Reside was not added";
        warning(scurf, smsg);
    }
    return p;
}

struct point_3d
{
    int x, y, z;
};

void correct_h_chains(pdb_file &pdbstruct)
{
    string scurf = "correct_h_chains", smsg = "";
    vector < vector < point_3d > > atom_points;
    vector < point_3d > chain_atoms;
    point_3d buf;
    vector <int> eraise_ch;
    double flen, minflen;
    int numl;
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        chain_atoms.clear();
        for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
            {
                for (int i = 0; i < pdbstruct.chains[ch].array[seg].atom[res].array.size(); ++i)
                {
                    buf.x = pdbstruct.chains[ch].array[seg].atom[res].array[i].x;
                    buf.y = pdbstruct.chains[ch].array[seg].atom[res].array[i].y;
                    buf.z = pdbstruct.chains[ch].array[seg].atom[res].array[i].z;
                    chain_atoms.push_back(buf);
                }
            }
        }
        atom_points.push_back(chain_atoms);
    }
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        if (atom_points[ch].size() == 0)
        {
            for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
            {
                for (int res = 0; res < pdbstruct.chains[ch].array[seg].hetatm.size(); ++res)
                {
                    numl = 0;
                    minflen = 10000000000;
                    for (int i = 0; i < pdbstruct.chains[ch].array[seg].hetatm[res].array.size(); ++i)
                    {
                        for (int j = 0; j < atom_points.size(); ++j)
                        {
                            for (int k = 0; k < atom_points[j].size(); ++k)
                            {
                                flen = distance(pdbstruct.chains[ch].array[seg].hetatm[res].array[i].x, pdbstruct.chains[ch].array[seg].hetatm[res].array[i].y, pdbstruct.chains[ch].array[seg].hetatm[res].array[i].z, atom_points[j][k].x, atom_points[j][k].y, atom_points[j][k].z);
                                if (flen < minflen)
                                {
                                    minflen = flen;
                                    numl = j;
                                }
                            }
                        }
                    }
                    if (append_hetatm(pdbstruct, pdbstruct.chains[numl].letter, pdbstruct.chains[ch].array[seg].hetatm[res]))
                    {
                        append_hetatm(pdbstruct, 'A', pdbstruct.chains[ch].array[seg].hetatm[res]);
                    }
                }
            }
            eraise_ch.push_back(ch);
        }
    }
    for (int i = eraise_ch.size() - 1; i >= 0; --i)
    {
        pdbstruct.chains.erase(pdbstruct.chains.begin() + eraise_ch[i]);
    }
    return;
}
bool check_chains(pdb_file &pdbstruct)
{
    string scurf = "check_chains", smsg = "";
    bool p = false;
    int amount;
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        p = false;
        amount = 0;
        for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
            {
                amount++;
                for (int i = 0; i < pdbstruct.chains[ch].array[seg].atom[res].array.size(); ++i)
                {
                    p = true;
                }
            }
        }
        if (amount < MINIMUM_RESIDUES_IN_CHAIN)
        {
            smsg = "There are less than ";
            smsg += inttostr_2(MINIMUM_RESIDUES_IN_CHAIN);
            smsg += " residues in chain ";
            smsg += pdbstruct.chains[ch].letter;
            warning(scurf, smsg);
        }
        if (!p)
        {
            return 1;
        }
    }
    return 0;
}

bool eq(double x, double p)
{
    return abs(p - x) < 0.0001;
}

void sline(uint8 *arr, int x1, int y1, int x2, int y2, uint8 color)
{
    int miny = min(y1, y2);
    int maxy = max(y1, y2);
    int minx = min(x1, x2);
    int maxx = max(x1, x2);
    double dx, dy;
    int x, y, curx, cury;
    if (x1 == x2)
    {
        for (int i = miny; i < maxy; i++)
        {
            arr[(i - 1) * PICTUREWIDTH + x1 - 1] = color;
        }
    }
    else if (y1 == y2)
    {
        for (int i = minx; i < maxx; i++)
        {
            arr[(y1 - 1) * PICTUREWIDTH + i - 1] = color;
        }
    }
    else
    {
        if (maxx - minx > maxy - miny)
        {
            dx = 1;
            dy = (maxx - minx) * 1. / (maxy - miny);
        }
        else if (maxx - minx == maxy - miny)
        {
            dy = 1;
            dx = (maxy - miny) * 1. / (maxx - minx);
        }
        x = 0;
        curx = x1;
        y = 0;
        cury = y1;
        arr[(cury - 1) * PICTUREWIDTH + curx - 1] = color;
        while (curx != x2 && cury != y2)
        {
            x++;
            y++;
            if (eq(x, dx))
            {
                x = 0;
                curx += 1 * (x1 < x2 ? 1 : -1);
            }
            else if (abs(dx - x) < 1)
            {
                x = -1;
                curx += 1 * (x1 < x2 ? 1 : -1);
            }
            else if (abs(dx - x) > 1)
            {
                x = 1;
                curx += 1 * (x1 < x2 ? 1 : -1);
            }
            if (eq(y, dy))
            {
                y = 0;
                cury += 1 * (y1 < y2 ? 1 : -1);
            }
            else if (abs(dy - y) < 1)
            {
                y = -1;
                cury += 1 * (y1 < y2 ? 1 : -1);
            }
            else if (abs(dy - y) > 1)
            {
                y = 1;
                cury += 1 * (y1 < y2 ? 1 : -1);
            }
            if (curx > maxx)
            {
                curx = maxx;
            }
            if (cury > maxy)
            {
                cury = maxy;
            }
            arr[(cury - 1) * PICTUREWIDTH + curx - 1] = color;
        }
    }
    return;
}

bool paintsymbol(uint8 *arr, int a, int x, int y, uint8 color, int mult = 1)
{
    bool digits[16][7][4] = {{{1,1,1,1},
    {1,0,0,1},
    {1,0,0,1},
    {1,0,0,1},
    {1,0,0,1},
    {1,0,0,1},
    {1,1,1,1}},
    {{0,0,1,0},
    {0,1,1,0},
    {1,0,1,0},
    {0,0,1,0},
    {0,0,1,0},
    {0,0,1,0},
    {1,1,1,1}},
    {{1,1,1,1},
    {0,0,0,1},
    {0,0,0,1},
    {0,0,1,0},
    {0,1,0,0},
    {1,0,0,0},
    {1,1,1,1}},
    {{1,1,1,1},
    {0,0,0,1},
    {0,0,0,1},
    {0,1,1,1},
    {0,0,0,1},
    {0,0,0,1},
    {1,1,1,1}},
    {{1,0,0,1},
    {1,0,0,1},
    {1,0,0,1},
    {1,1,1,1},
    {0,0,0,1},
    {0,0,0,1},
    {0,0,0,1}},
    {{1,1,1,1},
    {1,0,0,0},
    {1,0,0,0},
    {1,1,1,1},
    {0,0,0,1},
    {0,0,0,1},
    {1,1,1,1}},
    {{1,1,1,1},
    {1,0,0,0},
    {1,0,0,0},
    {1,1,1,1},
    {1,0,0,1},
    {1,0,0,1},
    {1,1,1,1}},
    {{1,1,1,1},
    {0,0,0,1},
    {0,0,0,1},
    {0,0,1,0},
    {0,0,1,0},
    {0,0,1,0},
    {0,0,1,0}},
    {{1,1,1,1},
    {1,0,0,1},
    {1,0,0,1},
    {1,1,1,1},
    {1,0,0,1},
    {1,0,0,1},
    {1,1,1,1}},
    {{1,1,1,1},
    {1,0,0,1},
    {1,0,0,1},
    {1,1,1,1},
    {0,0,0,1},
    {0,0,0,1},
    {1,1,1,1}},
    {{1,1,1,1},
    {1,0,0,0},
    {1,0,0,0},
    {1,1,1,1},
    {1,0,0,0},
    {1,0,0,0},
    {1,1,1,1}},
    {{1,0,0,1},
    {1,0,0,1},
    {1,0,0,1},
    {1,0,0,1},
    {0,1,1,0},
    {0,1,1,0},
    {0,1,1,0}},
    {{0,0,0,0},
    {0,0,0,0},
    {1,1,1,1},
    {0,0,0,0},
    {1,1,1,1},
    {0,0,0,0},
    {0,0,0,0}},
    {{1,1,1,1},
    {1,0,0,1},
    {0,0,0,1},
    {0,0,1,0},
    {0,0,1,0},
    {0,0,0,0},
    {0,0,1,0}},
    {{0,0,0,0},
    {0,0,0,0},
    {0,0,0,0},
    {1,1,1,1},
    {0,0,0,0},
    {0,0,0,0},
    {0,0,0,0}},
    {{0,0,0,1},
    {0,0,1,0},
    {0,0,1,0},
    {0,1,0,0},
    {0,1,0,0},
    {1,0,0,0},
    {1,0,0,0}}};
    if (a > 15 || a < 0)
    {
        return 1;
    }
    for (int i = 0; i < 7; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            if (digits[a][i][j])
            {
                for (int k = 0; k < mult; ++k)
                {
                    for (int p = 0; p < mult; ++p)
                    {
                        arr[(y + k - 1 + i * mult) * PICTUREWIDTH + x + p - 1 + j* mult] = color;
                    }
                }

            }

        }
    }
    return 0;
}
void paint_dot(uint8 *arr, int x, int y, int dsize, uint8 color)
{
    for (int j = y; j < y + dsize; ++j)
    {
        for (int i = x; i < x + dsize; ++i)
        {
            arr[(j - 1) * PICTUREWIDTH + i] = color;
        }
    }
    return;
}

void paintstr(uint8 *arr, string s, int x, int y, bool dir, uint8 color)
{
    //cout << s << "\n";
    string scurf = "paintstr", smsg = "";
    int curx = x;
    map <char, int> dist;
    dist['0'] = 0;
    dist['1'] = 1;
    dist['2'] = 2;
    dist['3'] = 3;
    dist['4'] = 4;
    dist['5'] = 5;
    dist['6'] = 6;
    dist['7'] = 7;
    dist['8'] = 8;
    dist['9'] = 9;
    dist['E'] = 10;
    dist['e'] = 10;
    dist['V'] = 11;
    dist['v'] = 11;
    dist['='] = 12;
    dist['?'] = 13;
    dist['-'] = 14;
    dist['/'] = 15;
    if (dir)
    {
        for (int i = 0; i < s.size(); ++i)
        {
            if (dist.count(s[i]))
            {
                paintsymbol(arr, dist[s[i]], curx, y - 14, color, 2);
                curx += 9 * (dir ? 1 : -1);
            }
            else if (s[i] == '.' || s[i] == ',')
            {
                paint_dot(arr, curx, y - 2, 2, color);
                curx += 4 * (dir ? 1 : -1);
            }
            else if (s[i] == ' ')
            {
                curx += 4 * (dir ? 1 : -1);
            }
            else
            {
                smsg = "There is no bitmap for symbol ";
                smsg += s[i];
                warning(scurf, smsg);
            }
        }
    }
    else
    {
        for (int i = s.size() - 1; i >= 0; --i)
        {
            if (dist.count(s[i]))
            {
                paintsymbol(arr, dist[s[i]], curx, y - 7, color, 2);
                curx += 9 * (dir ? 1 : -1);
            }
            else if (s[i] == '.' || s[i] == ',')
            {
                paint_dot(arr, curx, y - 2, 2, color);
                curx += 4 * (dir ? 1 : -1);
            }
            else if (s[i] == ' ')
            {
                curx += 4 * (dir ? 1 : -1);
            }
            else
            {
                smsg = "There is no bitmap for symbol ";
                smsg += s[i];
                warning(scurf, smsg);
            }
        }
    }
    return;
}

void paint_eplot(string file, vector < pair < double, double > > enpoints, int rgb = 1)
{
    int c_x = 0, c_y = 0, padding = 35;
    double square = 1;
    uint8 color = 0, backcolor = 255;
    uint8 carr[PICTURESIZE*rgb];
    for (int i = 0; i < PICTURESIZE*rgb; ++i)
    {
        carr[i] = backcolor;
    }
    //Chart axis y
    sline(carr, padding * 2, padding, padding * 2, PICTUREHEIGHT - padding * 2, color);
    sline(carr, padding * 2 + 10, padding + 10, padding * 2, padding, color);
    sline(carr, padding * 2 - 10, padding + 10, padding * 2, padding, color);
    //Chart axis x
    sline(carr, padding * 2, PICTUREHEIGHT - padding * 2, PICTUREWIDTH - padding, PICTUREHEIGHT - padding * 2, color);
    sline(carr, PICTUREWIDTH - padding, PICTUREHEIGHT - padding * 2, PICTUREWIDTH - padding - 10, PICTUREHEIGHT - padding * 2 - 10, color);
    sline(carr, PICTUREWIDTH - padding, PICTUREHEIGHT - padding * 2, PICTUREWIDTH - padding - 10, PICTUREHEIGHT - padding * 2 + 10, color);
    //Plot line
    if (eq(enpoints[0].first, 0))
    {
        c_x++;
    }
    if (eq(enpoints[0].second, 0))
    {
        c_y++;
    }
    sline(carr, padding * 2, PICTUREHEIGHT - padding * 2, padding * 2 + (PICTUREHEIGHT - padding * 4) * enpoints[0].first, PICTUREHEIGHT - padding * 2 - ((PICTUREHEIGHT - padding * 4) * enpoints[0].second), color);
    //Left stick
    sline(carr, padding * 2 - 5, PICTUREHEIGHT - padding * 2 - ((PICTUREHEIGHT - padding * 4) * enpoints[0].second), padding * 2 + 5, PICTUREHEIGHT - padding * 2 - ((PICTUREHEIGHT - padding * 4) * enpoints[0].second), color);
    //Bottom stick
    sline(carr, padding * 2 + (PICTUREWIDTH - padding * 4) * enpoints[0].first, PICTUREHEIGHT - padding * 2 - 5, padding * 2 + (PICTUREWIDTH - padding * 4) * enpoints[0].first, PICTUREHEIGHT - padding * 2 + 5, color);
    for (int i = 1; i < enpoints.size(); ++i)
    {
        if (eq(enpoints[i].first, enpoints[i - 1].first))
        {
            c_x++;
        }
        if (eq(enpoints[i].second, enpoints[i - 1].second))
        {
            c_y++;
        }
        square -= (enpoints[i].first - enpoints[i - 1].first) * (1 - enpoints[i].second);
        sline(carr, padding * 2 - 5, PICTUREHEIGHT - padding * 2 - ((PICTUREHEIGHT - padding * 4) * enpoints[i].second), padding * 2 + 5, PICTUREHEIGHT - padding * 2 - ((PICTUREHEIGHT - padding * 4) * enpoints[i].second), color);

        sline(carr, padding * 2 + (PICTUREWIDTH - padding * 4) * enpoints[i].first, PICTUREHEIGHT - padding * 2 - 5, padding * 2 + (PICTUREWIDTH - padding * 4) * enpoints[i].first, PICTUREHEIGHT - padding * 2 + 5, color);

        //Plot line
        sline(carr, padding * 2 + (PICTUREWIDTH - padding * 4) * enpoints[i - 1].first, PICTUREHEIGHT - padding * 2 - ((PICTUREHEIGHT - padding * 4) * enpoints[i - 1].second), padding * 2 + (PICTUREWIDTH - padding * 4) * enpoints[i].first, PICTUREHEIGHT - padding * 2 - ((PICTUREHEIGHT - padding * 4) * enpoints[i].second), color);
    }
    paintstr(carr, "EV = " + doubtostr_3(square, 3), PICTUREWIDTH - padding * 3, padding - 10, 1, color);
    paintstr(carr, inttostr_2(c_x) + " / " + inttostr_2(c_x), padding * 2 - 20, padding + 5, 0, color);
    paintstr(carr, inttostr_2(c_y) + " / " + inttostr_2(c_y), PICTUREWIDTH - padding * 2 - 20, PICTUREHEIGHT - padding * 2 + 25, 1, color);
    compress_image_to_jpeg_file(file.c_str(), PICTUREWIDTH, PICTUREHEIGHT, rgb, carr);
    return;
}

map<string, bool> parse_yn_ing(string nrp_parf)
{
    string scurf = "parse_yn_ing", smsg = "";
    map<string, bool> res;
    string s, ss;
    if (!file_exists(nrp_parf))
    {
        smsg = "There is no file ";
        smsg += nrp_parf;
        error_msg(scurf, smsg);
        exit(1);
    }
    ifstream yin(nrp_parf.c_str());
    while (getline(yin, s))
    {
        ss = s.substr(s.find('\t') + 1, s.size() - s.find('\t') - 1);
        s = s.substr(0, s.find('\t'));
        if (ss.find("YES") != string::npos || ss.find("yes") != string::npos || ss == "1")
        {
            if (res.count(s) == 0)
            {
                res[s] = true;
            }
            else
            {
                if (res[s] != true)
                {
                    smsg = "Different values for compound ";
                }
                else
                {
                    smsg = "Double values for compound ";
                }
                smsg += s;
                warning(scurf, smsg);
            }
        }
        else
        {
            if (res.count(s) == 0)
            {
                res[s] = false;
            }
            else
            {

                if (res[s] != true)
                {
                    smsg = "Different values for compound ";
                }
                else
                {
                    smsg = "Double values for compound ";
                }
                smsg += s;
                warning(scurf, smsg);
            }
        }
    }
    yin.close();
    return res;
}

bool cmp2(pair<string, int> a, pair<string, int> b)
{
    return a.second < b.second;
}

void rampce_dock_results(string nrp_parf, string nrp_resf, string plot_file)
{
    string scurf = "rampce_dock_results", smsg = "";
    map < string, bool> accord = parse_yn_ing(nrp_parf);
    if (!file_exists(nrp_resf))
    {
        smsg = "There is no file ";
        smsg += nrp_resf;
        error_msg(scurf, smsg);
        exit(1);
    }
    ifstream rin(nrp_resf.c_str());
    string s, ss, curpdb, curchain, comp;
    string path, name, ext;
    int k = 0, pc_num = -1, ch_num = -1, co_num = -1, eg_num = -1, z = 0;
    int n = -1;
    int amnot, amyes;
    double square = 1;
    double x, y;
    set <string> misscomp;
    vector <pair < string, int> > res_pdb;
    vector < pair < double, double > > xypoints;
    getline(rin, s);
    filename_parse(plot_file, path, name, ext);
    if (ext.size() == 0)
    {
        ext = "jpg";
    }
    ss = get_file_name(path, name + SUFF, "txt.xls");

    ofstream rout(ss.c_str());
    rout.imbue(locale(cout.getloc(), new DecimalSeparator<char>(',')));
    rout << "PDBCode\tChain\tCompound\tEnergy\tRank_x\tRank_Y\tEV\n";
    ss = next_lex(s, k, '\t');
    z = 0;

    while (ss.size() > 0)
    {
        if (ss == "PDBCode")
        {
            pc_num = z;
        }
        else if (ss == "Chain" || ss == "chain")
        {
            ch_num = z;
        }
        else if (ss == "Compound" || ss == "par8")
        {
            co_num = z;
        }
        else if (ss == "1st best delta G" || ss == "Energy" || ss == "delta G" || ss == "1 delta G")
        {
            eg_num = z;
        }
        ss = next_lex(s, k, '\t');
        z++;
    }
    if (pc_num == -1 || ch_num == -1 || co_num == -1 || eg_num == -1)
    {
        smsg = "Can't find one of required columns (PDBCode, Chain, Compound, Energy/delta G)\n";
        error_msg(scurf, smsg);
        return;
    }
    while (getline(rin, s))
    {
        if (s.find("File with no results") != string::npos)
        {
            continue;
        }
        else
        {
            break;
        }
    }
    k = 0;
    for (int i = 0; i < z; ++i)
    {
        if (i == pc_num)
        {
            curpdb = next_lex(s, k, '\t');
        }
        else if (i == ch_num)
        {
            curchain = next_lex(s, k, '\t');
        }
        else if (i == co_num)
        {
            comp = next_lex(s, k, '\t');
        }
        else if (i == eg_num)
        {
            ss = next_lex(s, k, '\t');
            res_pdb.push_back(make_pair(comp, strtodoub(ss)));
        }
        else
        {
            ss = next_lex(s, k, '\t');
        }
    }
    while(getline(rin, s))
    {
        if (s.find("File with no results") != string::npos)
        {
            continue;
        }
        //cout << s << "\n";
        k = 0;
        for (int i = 0; i < z; ++i)
        {
            if (i == pc_num)
            {
                ss = next_lex(s, k, '\t');
                //cout << ss << "\n";
                if (ss != curpdb)
                {
                    sort(res_pdb.begin(), res_pdb.end(), cmp2);
                    amnot = 0;
                    amyes = 0;
                    xypoints.clear();
                    //cout <<  res_pdb.size() << "\n";
                    for (int j = 0; j < res_pdb.size(); ++j)
                    {
                        //cout << res_pdb[j].first << "\n";
                        if (accord.count(res_pdb[j].first))
                        {
                            if (amyes + amnot > 0)
                            {
                                rout << "\n";
                            }
                            rout << curpdb << "\t" << curchain << "\t" << res_pdb[j].first << "\t" << res_pdb[j].second / 100. << "\t";
                            if (accord[res_pdb[j].first])
                            {
                                amyes++;
                                rout << "-" << "\t" << amyes;
                            }
                            else
                            {
                                amnot++;
                                rout << amnot << "\t" << "-";
                            }
                            xypoints.push_back(make_pair(amnot, amyes));
                        }
                        else
                        {
                            misscomp.insert(res_pdb[j].first);
                        }
                    }
                    //cout << amnot << " " << amyes << " " << xypoints.size() << "\n";
                    if (amnot + amyes == 0)
                    {
                        smsg = "There is no description for any compounds for pdb ";
                        smsg += curpdb;
                        smsg += " chain ";
                        smsg += curchain;
                        smsg += "\n";
                        warning(scurf, smsg);
                    }
                    else
                    {
                        if (amnot != 0)
                        {
                            xypoints[0].first /= (amnot * 1.);
                        }
                        if (amyes != 0)
                        {
                            xypoints[0].second /= (amyes * 1.);
                        }
                        square = 1.;
                        for (int j = 1; j < xypoints.size(); ++j)
                        {
                            if (amnot != 0)
                            {
                                xypoints[j].first /= (amnot * 1.);
                            }
                            if (amyes != 0)
                            {
                                xypoints[j].second /= (amyes * 1.);
                            }
                            square -= (xypoints[j].first - xypoints[j - 1].first) * (1 - xypoints[j].second);
                            //cout << xypoints[j].first << " " << xypoints[j].second << "\n";
                        }
                        rout << "\t" << doubtostr_3(square, 3) << "\n";
                        paint_eplot(get_file_name(path, string(name + "_" + curpdb + "_" + curchain + SUFF), ext), xypoints);
                    }
                    curpdb = ss;
                    curchain = "-";
                    res_pdb.clear();
                }
            }
            else if (i == ch_num)
            {
                ss = next_lex(s, k, '\t');
                //cout << ss << "\n";
                if (curchain == "-")
                {
                    curchain = ss;
                }
                else if (ss != curchain)
                {
                    sort(res_pdb.begin(), res_pdb.end(), cmp2);
                    amnot = 0;
                    amyes = 0;
                    xypoints.clear();
                    for (int j = 0; j < res_pdb.size(); ++j)
                    {
                        //cout << res_pdb[j].first << "\n";
                        if (accord.count(res_pdb[j].first))
                        {
                            if (amyes + amnot > 0)
                            {
                                rout << "\n";
                            }
                            rout << curpdb << "\t" << curchain << "\t" << res_pdb[j].first << "\t" << res_pdb[j].second / 100. << "\t";
                            if (accord[res_pdb[j].first])
                            {
                                amyes++;
                                rout << "-" << "\t" << amyes;
                            }
                            else
                            {
                                amnot++;
                                rout << amnot << "\t" << "-";
                            }
                            xypoints.push_back(make_pair(amnot, amyes));
                        }
                        else
                        {
                            misscomp.insert(res_pdb[j].first);
                        }
                    }
                    //cout << amnot << " " << amyes << " " << xypoints.size() << "\n";
                    if (amnot + amyes == 0)
                    {
                        smsg = "There is no description for any compounds for pdb ";
                        smsg += curpdb;
                        smsg += " chain ";
                        smsg += curchain;
                        smsg += "\n";
                        warning(scurf, smsg);
                    }
                    else
                    {
                        if (amnot != 0)
                        {
                            xypoints[0].first /= (amnot * 1.);
                        }
                        if (amyes != 0)
                        {
                            xypoints[0].second /= (amyes * 1.);
                        }
                        square = 1;
                        for (int j = 1; j < xypoints.size(); ++j)
                        {
                            if (amnot != 0)
                            {
                                xypoints[j].first /= (amnot * 1.);
                            }
                            if (amyes != 0)
                            {
                                xypoints[j].second /= (amyes * 1.);
                            }
                            square -= (xypoints[j].first - xypoints[j - 1].first) * (1 - xypoints[j].second);
                            //cout << xypoints[j].first << " " << xypoints[j].second << "\n";
                        }
                        rout << "\t" << doubtostr_3(square, 3) << "\n";
                        paint_eplot(get_file_name(path, string(name + "_" + curpdb + "_" + curchain + SUFF), ext), xypoints);
                    }
                    curchain = ss;
                    res_pdb.clear();
                }
            }
            else if (i == co_num)
            {
                ss = next_lex(s, k, '\t');
                comp = ss;
                //cout << ss << "\n";
            }
            else if (i == eg_num)
            {
                ss = next_lex(s, k, '\t');
                res_pdb.push_back(make_pair(comp, strtodoub(ss)));
            }
            else
            {
                ss = next_lex(s, k, '\t');
            }

        }
    }
    sort(res_pdb.begin(), res_pdb.end(), cmp2);
    amnot = 0;
    amyes = 0;
    xypoints.clear();
    for (int j = 0; j < res_pdb.size(); ++j)
    {
        if (accord.count(res_pdb[j].first))
        {
            if (amyes + amnot > 0)
            {
                rout << "\n";
            }
            rout << curpdb << "\t" << curchain << "\t" << res_pdb[j].first << "\t" << res_pdb[j].second / 100. << "\t";
            if (accord[res_pdb[j].first])
            {
                amyes++;
                rout << "-" << "\t" << amyes;
            }
            else
            {
                amnot++;
                rout << amnot << "\t" << "-";
            }
            xypoints.push_back(make_pair(amnot, amyes));
        }
        else
        {
            misscomp.insert(res_pdb[j].first);
        }
    }
    if (amnot + amyes == 0)
    {
        smsg = "There is no description for any compounds for pdb ";
        smsg += curpdb;
        smsg += " chain ";
        smsg += curchain;
        smsg += "\n";
        warning(scurf, smsg);
    }
    else
    {
        if (amnot != 0)
        {
            xypoints[0].first /= (amnot * 1.);
        }
        if (amyes != 0)
        {
            xypoints[0].second /= (amyes * 1.);
        }
        square = 1;
        for (int j = 1; j < xypoints.size(); ++j)
        {
            if (amnot != 0)
            {
                xypoints[j].first /= (amnot * 1.);
            }
            if (amyes != 0)
            {
                xypoints[j].second /= (amyes * 1.);
            }
            square -= (xypoints[j].first - xypoints[j - 1].first) * (1 - xypoints[j].second);
            //cout << xypoints[j].first << " " << xypoints[j].second << "\n";
        }
        rout << "\t" << doubtostr_3(square, 3) << "\n";
        paint_eplot(get_file_name(path, string(name + "_" + curpdb + "_" + curchain + SUFF), ext), xypoints);
    }
    rin.close();
    rout.close();
    if (misscomp.size() > 0)
    {

        smsg = "There are no ";
        for (set<string>::iterator it = misscomp.begin(); it != misscomp.end(); ++it)
        {
            smsg += *it;
            smsg += " ";
        }
        smsg += "compounds in parameters file. Skipping it";
        warning(scurf, smsg);
    }
    return;
}

double calculate_RMSD(residue res1, residue res2, int type = 0)
{
    string scurf = "calculate_RMSD", smsg = "";
    double result = 0;
    map<string, int> atmcrd1, atmcrd2;
    for (int i = 0; i < res1.array.size(); ++i)
    {
        if (res1.array[i].element != "h" && res1.array[i].element != "H")
        {
            if (atmcrd1.count(res1.array[i].atom_name))
            {
                smsg = "There are not unic atom names in 1st residue (";
                smsg += res1.array[i].atom_name;
                smsg += "). Terminated";
                warning(scurf, smsg);
                return -1;
            }
            atmcrd1[res1.array[i].atom_name] = i;
        }
    }
    for (int i = 0; i < res2.array.size(); ++i)
    {
        if (res2.array[i].element != "h" && res2.array[i].element != "H")
        {
            if (atmcrd2.count(res2.array[i].atom_name))
            {
                smsg = "There are not unic atom names in 2nd residue (";
                smsg += res2.array[i].atom_name;
                smsg += "). Terminated";
                warning(scurf, smsg);
                return -1;
            }
            atmcrd2[res2.array[i].atom_name] = i;
        }
    }
    if (atmcrd1.size() != atmcrd2.size())
    {
        smsg = "Two residues have different amount of non-H atoms ";
        smsg += inttostr_2(atmcrd1.size());
        smsg += " and ";
        smsg += inttostr_2(atmcrd2.size());
        smsg += ". Terminated";
        warning(scurf, smsg);
        return -1;
    }
    for (int i = 0; i < res2.array.size(); ++i)
    {
        if (res2.array[i].element != "h" && res2.array[i].element != "H")
        {
            if (atmcrd1.count(res2.array[i].atom_name))
            {
                result += distance(res1.array[atmcrd1[res2.array[i].atom_name]].x, res1.array[atmcrd1[res2.array[i].atom_name]].y, res1.array[atmcrd1[res2.array[i].atom_name]].z, res2.array[i].x, res2.array[i].y, res2.array[i].z);
            }
            else
            {
                smsg = "There is no atom with atom name ";
                smsg += res2.array[i].atom_name;
                smsg += " in 1st residue. Terminated";
                warning(scurf, smsg);
                return -1;
            }
        }
    }
    return result;
}

bool read_residue(string file, residue &nresi)
{
    string scurf = "read_residue", smsg = "";
    pdb_file buf;
    parse_pdb(file, buf);
    nresi.number = "@";
    for (int ch = 0; ch < buf.chains.size(); ++ch)
    {
        for (int seg = 0; seg < buf.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < buf.chains[ch].array[seg].atom.size(); ++res)
            {
                if (nresi.number == "@")
                {
                    nresi = buf.chains[ch].array[seg].atom[res];
                }
                else
                {
                    smsg = "There is more than 1 residue in input file ";
                    smsg += file;
                    smsg += ". Terminated";
                    warning(scurf, smsg);
                    return 1;
                }
            }
            for (int res = 0; res < buf.chains[ch].array[seg].hetatm.size(); ++res)
            {
                if (nresi.number == "@")
                {
                    nresi = buf.chains[ch].array[seg].hetatm[res];
                }
                else
                {
                    smsg = "There is more than 1 residue in input file ";
                    smsg += file;
                    smsg += ". Terminated";
                    warning(scurf, smsg);
                    return 1;
                }

            }
        }
    }
    return 0;
}

double lig_RMSD(string file1, string file2)
{
    residue res1, res2;
    if (read_residue(file1, res1) || read_residue(file2, res2))
    {
        return -1;
    }
    return calculate_RMSD(res1, res2)/1000.;
}


bool parse_mol2_file(string file, mol2file &mstruct)
{
    string scurf = "parse_mol2_file", smsg = "";
    ifstream mlin(file.c_str());
    string s, ss, k;
    atom buf;
    mol2bond bbuf;
    int natoms = 0, nbonds = 0;
    int mols = 0, atoms = 0, bonds = 0;
    while (getline(mlin, s))
    {
        if (s[0] == '@')
        {
            ss = s.substr(0, min(s.size(), 17));
            mols *= 2;
            atoms *= 2;
            bonds *= 2;
            if (ss == "@<TRIPOS>MOLECULE")
            {
                mols = 1;
            }
            else if (ss == "@<TRIPOS>ATOM")
            {
                atoms = 1;
            }
            else if (ss == "@<TRIPOS>BOND")
            {
                bonds = 1;
            }
        }
        else
        {
            if (mols == 1)
            {
                mstruct.mol_name = s;
                if (!getline(mlin, s))
                {
                    smsg = "Molecule section has ended unexpectably. Reading terminated";
                    warning(scurf, smsg);
                    return 1;
                }
                if (!getline(mlin, s))
                {
                    smsg = "Molecule section has ended unexpectably. Reading terminated";
                    warning(scurf, smsg);
                    return 1;
                }
                mstruct.mol_type = s;
                if (s != "SMALL")
                {
                    smsg = "Type of molecule is not SMALL";
                    warning(scurf, smsg);
                }
                if (!getline(mlin, s))
                {
                    smsg = "Molecule section has ended unexpectably. Reading terminated";
                    warning(scurf, smsg);
                    return 1;
                }
                mstruct.ch_type = s;
                if (s != "GASTEIGER" && s != "NO_CHARGES" && s != "GAST_HUCK")
                {
                    smsg = "Unknown charges type";
                    warning(scurf, smsg);
                }
                mols *= 2;
            }
            else if (atoms == 1)
            {
                buf.serial = strtoint(trimdel(s.substr(0, 7)));
                buf.atom_name = trimdel(s.substr(8, 4));
                buf.x = strtodoub(trimdel(s.substr(16, 10)));
                buf.y = strtodoub(trimdel(s.substr(26, 10)));
                buf.z = strtodoub(trimdel(s.substr(36, 10)));
                buf.element = trimdel(s.substr(47, 4));
                buf.resnumb = trimdel(s.substr(51, 5));
                buf.resname = trimdel(s.substr(58, 10));
                buf.ch_val = strtodoub(trimdel(s.substr(69, 7)));
                mstruct.atoms.push_back(buf);
            }
            else if (bonds == 1)
            {
                bbuf.a = strtoint(trimdel(s.substr(6, 6)));
                bbuf.b = strtoint(trimdel(s.substr(12, 6)));
                bbuf.type = trimdel(s.substr(18, 6));
                mstruct.bonds.push_back(bbuf);
            }
        }
    }
    if (mstruct.atoms.size() < 1)
    {
        smsg = "There is no atoms in file. Terminated";
        warning(scurf, smsg);
        return 1;
    }
    if (mstruct.bonds.size() < 1)
    {
        smsg = "There is no bonds in file. Something wrong";
        warning(scurf, smsg);
    }
    return 0;
}

void out_mol2_file(string file, mol2file &mstruct)
{
    ofstream mlout(file.c_str());
    string s;
    mlout << "@<TRIPOS>MOLECULE\n";
    //mlout << mstruct.mol_name << ". Generated by " << PARSERNAME << "\n";
    mlout << " " << mstruct.atoms.size() << " " << mstruct.bonds.size() << " 0 0 0\n";
    mlout << mstruct.mol_type << "\n";
    mlout << mstruct.ch_type << "\n";
    mlout << "\n";
    mlout << "@<TRIPOS>ATOM\n";
    for (int i = 0; i < mstruct.atoms.size(); ++i)
    {
        s = inttostr(mstruct.atoms[i].serial, 7, true);
        s += " ";
        if (mstruct.atoms[i].atom_name.size() < 4)
        {
            s += " ";
            s += formatstr(mstruct.atoms[i].atom_name, 3, false);
        }
        else
        {
            s += mstruct.atoms[i].atom_name;
        }
        s += "    ";
        s += doubtostr(mstruct.atoms[i].x, 10, 4);
        s += doubtostr(mstruct.atoms[i].y, 10, 4);
        s += doubtostr(mstruct.atoms[i].z, 10, 4);
        s += " ";
        s += formatstr(mstruct.atoms[i].element, 4, false);
        s += formatstr(mstruct.atoms[i].resnumb, 5, true);
        s += "  ";
        s += formatstr(mstruct.atoms[i].resname, 10, false);
        s += " ";
        s += doubtostr(mstruct.atoms[i].ch_val, 7, 4);
        mlout << s << "\n";
    }
    mlout << "@<TRIPOS>BOND\n";
    for (int i = 0; i < mstruct.bonds.size(); ++i)
    {
        s = inttostr(i + 1, 6, true);
        s += inttostr(mstruct.bonds[i].a, 6, true);
        s += inttostr(mstruct.bonds[i].b, 6, true);
        s += formatstr(mstruct.bonds[i].type, 5, true);
        mlout << s << "\n";
    }

    return;
}
/*
int dfs_pass(vector <bool> &visited, vector <vector <int> > &edges, vector <int> &ico, vector<vector<int> > &nedges, int parent = -1, int vertex = 0, int depth = 0)
{
    //cout << depth << "\n";
    int a, ret = vertex, group = vertex, branch = 0;
    vector <int> leaves, buf;
    leaves.clear();
    buf.clear();
    if (edges[vertex].size() == 1)
    {
        //cout << vertex + 1 << " -2\n";
        return -2;
    }
    for (int i = 0; i < edges[vertex].size(); ++i)
    {
        if (edges[vertex][i] != parent)
        {
            if (!visited[edges[vertex][i]])
            {
                visited[edges[vertex][i]] = true;
                a = dfs_pass(visited, edges, ico, nedges, vertex, edges[vertex][i], depth + 1);
                if (a == -2)
                {
                    leaves.push_back(edges[vertex][i]);
                }
                else if (a == vertex)
                {
                    ret = -1;
                    group = vertex;
                }
                else if (a == -1)
                {
                    //cout << ico[edges[vertex][i]] << "\n";
                    buf.push_back(ico[edges[vertex][i]]);
                    branch++;
                }
                else
                {
                    ret = a;
                    group = a;
                }
            }
            else
            {
                if (group == ico[edges[vertex][i]])
                {
                    ret = -1;
                }
                else
                {
                    group = edges[vertex][i];
                    ret = edges[vertex][i];
                }
            }
        }
    }
    for (int i = 0; i < buf.size(); ++i)
    {
        nedges[group].push_back(buf[i]);
    }
    ico[vertex] = group;
    for (int i = 0; i < leaves.size(); ++i)
    {
        ico[leaves[i]] = group;
    }
    if (edges[vertex].size() - 1 == branch + leaves.size())
    {
        //cout << vertex + 1 << " " << "-1\n";
        return -1;
    }
    //cout << vertex + 1<< " " << ret + 1 << "\n";
    return ret;
}*/

/*void compute_magic_rmsd(string file1, string file2)
{

    mol2file ml2f;
    parse_mol2_file(file1, ml2f);
    vector <bool> visited;
    vector <int> ico;
    vector <vector <int> > edges, nedges;
    vector <int> buf;
    for (int i = 0; i < ml2f.atoms.size(); ++i)
    {
        visited.push_back(false);
        ico.push_back(-1);
        edges.push_back(buf);
        nedges.push_back(buf);
    }
    for (int i = 0; i < ml2f.bonds.size(); ++i)
    {
        //cout << ml2f.bonds[i].a << " " << ml2f.bonds[i].b << "\n";
        edges[ml2f.bonds[i].a - 1].push_back(ml2f.bonds[i].b - 1);
        edges[ml2f.bonds[i].b - 1].push_back(ml2f.bonds[i].a - 1);
    }
    visited[0] = true;
    dfs_pass(visited, edges, ico, nedges, -1, 0, 0);
    int cur, k = 0;
    map<int, vector<int> > ssss;
    for (int i = 0; i < ml2f.atoms.size(); ++i)
    {
        if (!ssss.count(ico[i]))
        {
            ssss[ico[i]] = buf;
        }
        ssss[ico[i]].push_back(i);
    }
    for (map<int, vector<int> >::iterator it = ssss.begin(); it != ssss.end(); ++it)
    {
        cout << it->first << "  ";
        for (int i = 0; i < it->second.size(); ++i)
        {
            cout << ml2f.atoms[it->second[i]].atom_name << " ";
        }
        cout << "\n";
    }
    for (int i = 0; i < ml2f.atoms.size(); ++i)
    {
        if (nedges[i].size() != 0)
        {
            cout << i << " ";
            for (int j = 0; j < nedges[i].size(); ++j)
            {
                cout << nedges[i][j] << " ";
            }
            cout << "\n";
        }
    }
    out_mol2_file("2.mol2", ml2f);
}*/


void make_fasta(pdb_file pdbstruct, string fasta_file)
{

    string path, name, ext;
    filename_parse(fasta_file, path, name, ext);
    fasta_file = get_file_name(path, name, ext);
    ofstream fout(fasta_file.c_str());
    int counter = 0;
    map <string, char> amino_acids;
    amino_acids["ALA"]='A';
    amino_acids["ARG"]='R';
    amino_acids["ASN"]='N';
    amino_acids["ASP"]='D';
    amino_acids["CYS"]='C';
    amino_acids["GLU"]='E';
    amino_acids["GLN"]='Q';
    amino_acids["GLY"]='G';
    amino_acids["HIS"]='H';
    amino_acids["ILE"]='I';
    amino_acids["LEU"]='L';
    amino_acids["LYS"]='K';
    amino_acids["MET"]='M';
    amino_acids["PHE"]='F';
    amino_acids["PRO"]='P';
    amino_acids["SER"]='S';
    amino_acids["THR"]='T';
    amino_acids["TRP"]='W';
    amino_acids["TYR"]='Y';
    amino_acids["VAL"]='V';
    amino_acids["SEC"]='U';
    amino_acids["PYL"]='O';
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        fout << ">" << fasta_file << "_chain_" << pdbstruct.chains[ch].letter << " generated by PDBParser. Non-standard residues are framed in ().";
        for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
            {
                if (counter == 0)
                {
                    fout << " 1st amino acid has residue number " << pdbstruct.chains[ch].array[seg].atom[res].number << "\n";
                    counter++;
                }
                if (amino_acids.count(pdbstruct.chains[ch].array[seg].atom[res].array[0].resname))
                {
                    fout << amino_acids[pdbstruct.chains[ch].array[seg].atom[res].array[0].resname];
                }
                else
                {
                    fout << "(" << pdbstruct.chains[ch].array[seg].atom[res].array[0].resname << ")";
                }

                counter++;
                if (counter > 70)
                {
                    fout << "\n";
                    counter = 1;
                }
            }
        }
        fout << "\n";
        counter = 0;
    }
    fout.close();
}

void rotate_pdb_coords(pdb_file & pdbstruct, int x_a, int y_a, int z_a)
{
    int x, y, z;
    double angle;
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
            {
                for (int i = 0; i < pdbstruct.chains[ch].array[seg].atom[res].array.size(); ++i)
                {
                    if (x_a != 0)
                    {
                        x = pdbstruct.chains[ch].array[seg].atom[res].array[i].x;
                        y = pdbstruct.chains[ch].array[seg].atom[res].array[i].y;
                        z = pdbstruct.chains[ch].array[seg].atom[res].array[i].z;
                        angle = x_a * M_PI / 180;
                        pdbstruct.chains[ch].array[seg].atom[res].array[i].y = y * cos(angle) - z * sin(angle);
                        pdbstruct.chains[ch].array[seg].atom[res].array[i].z = y * sin(angle) + z * cos(angle);
                    }
                    if (y_a != 0)
                    {
                        x = pdbstruct.chains[ch].array[seg].atom[res].array[i].x;
                        y = pdbstruct.chains[ch].array[seg].atom[res].array[i].y;
                        z = pdbstruct.chains[ch].array[seg].atom[res].array[i].z;
                        angle = y_a * M_PI / 180;
                        pdbstruct.chains[ch].array[seg].atom[res].array[i].x = x * cos(angle) - z * sin(angle);
                        pdbstruct.chains[ch].array[seg].atom[res].array[i].z = x * sin(angle) + z * cos(angle);
                    }
                    if (z_a != 0)
                    {
                        x = pdbstruct.chains[ch].array[seg].atom[res].array[i].x;
                        y = pdbstruct.chains[ch].array[seg].atom[res].array[i].y;
                        z = pdbstruct.chains[ch].array[seg].atom[res].array[i].z;
                        angle = z_a * M_PI / 180;
                        pdbstruct.chains[ch].array[seg].atom[res].array[i].x = x * cos(angle) - y * sin(angle);
                        pdbstruct.chains[ch].array[seg].atom[res].array[i].y = x * sin(angle) + y * cos(angle);
                    }
                }
            }
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].hetatm.size(); ++res)
            {
                for (int i = 0; i < pdbstruct.chains[ch].array[seg].hetatm[res].array.size(); ++i)
                {
                    if (x_a != 0)
                    {
                        x = pdbstruct.chains[ch].array[seg].hetatm[res].array[i].x;
                        y = pdbstruct.chains[ch].array[seg].hetatm[res].array[i].y;
                        z = pdbstruct.chains[ch].array[seg].hetatm[res].array[i].z;
                        angle = x_a * M_PI / 180;
                        pdbstruct.chains[ch].array[seg].hetatm[res].array[i].y = y * cos(angle) - z * sin(angle);
                        pdbstruct.chains[ch].array[seg].hetatm[res].array[i].z = y * sin(angle) + z * cos(angle);
                    }
                    if (y_a != 0)
                    {
                        x = pdbstruct.chains[ch].array[seg].hetatm[res].array[i].x;
                        y = pdbstruct.chains[ch].array[seg].hetatm[res].array[i].y;
                        z = pdbstruct.chains[ch].array[seg].hetatm[res].array[i].z;
                        angle = y_a * M_PI / 180;
                        pdbstruct.chains[ch].array[seg].hetatm[res].array[i].x = x * cos(angle) - z * sin(angle);
                        pdbstruct.chains[ch].array[seg].hetatm[res].array[i].z = x * sin(angle) + z * cos(angle);
                    }
                    if (z_a != 0)
                    {
                        x = pdbstruct.chains[ch].array[seg].hetatm[res].array[i].x;
                        y = pdbstruct.chains[ch].array[seg].hetatm[res].array[i].y;
                        z = pdbstruct.chains[ch].array[seg].hetatm[res].array[i].z;
                        angle = z_a * M_PI / 180;
                        pdbstruct.chains[ch].array[seg].hetatm[res].array[i].x = x * cos(angle) - y * sin(angle);
                        pdbstruct.chains[ch].array[seg].hetatm[res].array[i].y = x * sin(angle) + y * cos(angle);
                    }
                }
            }
        }
    }
    return;
}

void cut_sdf_to_mol2(string sdf_m_file, string sdf_out_dir)
{
    string scurf = "cut_sdf_to_mol2", smsg = "";
    string path, name, ext, s;
    filename_parse(sdf_out_dir, path, name, ext);
    if (!file_exists(sdf_m_file))
    {
        smsg = "There is no file ";
        smsg += sdf_m_file;
        warning(scurf, smsg);
        return;
    }
    ifstream sin(sdf_m_file.c_str());
    vector <string> buf;
    buf.clear();
    bool p = 0;
    while (getline(sin, s))
    {
        if (s.substr(0, 17) == "@<TRIPOS>MOLECULE")
        {
            p = 1;
            if (buf.size())
            {
                out_lines(buf, path + name + ".mol2");
                buf.clear();
            }
        }
        else if (p)
        {
            name = s;
            p = 0;

        }
        buf.push_back(s);
    }
    if (buf.size())
    {
        out_lines(buf, path + name + ".mol2");
    }
    return;
}


void rename_atoms_pat(pdb_file & pdbstruct, pdb_file & pdbstruct_ren, bool sort_res)
{
    /*string scurf = "rename_atoms_pat", smsg;
    vector <int> minrmsd;
    double curmin, curm;
    int num;
    for (int ch = 0; ch < pdbstruct.chains.size(); ++ch)
    {
        for (int seg = 0; seg < pdbstruct.chains[ch].array.size(); ++seg)
        {
            if (pdbstruct.chains[ch].array[seg].atom.size() != pdbstruct_ren.chains[ch].array[seg].atom.size())
            {
                smsg = "Different residue amount in segment ";
                smsg += pdbstruct.chains[ch].array[i].cname;
                smsg += ". Terminated\n";
                error_msg(scurf, smsg);
                exit(0);
            }
            for (int res = 0; res < pdbstruct.chains[ch].array[seg].atom.size(); ++res)
            {
                if (pdbstruct.chains[ch].array[seg].atom[res].number != pdbstruct_ren.chains[ch].array[seg].atom[res].number)
                {
                    smsg = "Different resnumbers\n";
                    warning(scurf, smsg);
                }
                for (int i = 0; i < pdbstruct.chains[ch].array[seg].atom[res].array.size(); ++i)
                {
                    curmin = 10000000;
                    curmin = 10000000;
                    num = -1;
                    for (int j = 0; j < pdbstruct_ren.chains[ch].array[seg].atom[res].array.size(); ++j)
                    {
                        if (pdbstruct.chains[ch].array[seg].atom[res].array[i].element == pdbstruct_ren.chains[ch].array[seg].atom[res].array[j].element)
                        {
                            curm = distance(pdbstruct.chains[ch].array[seg].atom[res].array[i].x, pdbstruct.chains[ch].array[seg].atom[res].array[i].y, pdbstruct.chains[ch].array[seg].atom[res].array[i].z, pdbstruct_ren.chains[ch].array[seg].atom[res].array[j].x, pdbstruct_ren.chains[ch].array[seg].atom[res].array[j].y, pdbstruct_ren.chains[ch].array[seg].atom[res].array[j].z);
                            if (curmin > curm)
                            {
                                curmin = curm;
                                num = j;
                            }
                        }
                    }
                    minrmsd.push_back(num);
                }
            }
        }
    }*/
}

int main(int argc, char ** argv)
{
	string scurf = "main", smsg = "", PROGRAMMNAME, PROGRAMMPATH;
	string proin = "", proren, ligout, proout, hetfile, ss, basefile, dlg_list_file, dlg_out_file, lig_cor_sample, lig_cor_coord, fasta_file, postf, dlg_input, dlg_output;
	string ciffile, ciflig, nrp_parf, nrp_resf, plot_file, crmsdfile1, crmsdfile2, proout_ren, ligname = "", lmol2_file = "", sdf_m_file = "", sdf_out_dir = "";

    char ligch;
    int call;
	vector <string> pars_is, redock;
	bool pin = 0, pout = 0, lout = 0, plig = 0, hf = 0, bf = 0, sort_res = 0, cut = 0, FASTA = 0, update_base_file = 0, vina_log_parse = 0, SET_OC_TF_1 = 0, DEL_ALL_HYDROGENS = 0, sdfc_mode = 0;
	bool dlg_parse = 0, dlg_parse_runlist = 0, lig_cor = 0, DCLUST = 0, cif = 0, chc = 0, roc_plot = 0, crmsd = 0, renatpdb = 0, vina_log_parse_runlist = 0, RES_RENAME = 0;
	bool rotate_pdb = 0;
	int renumb_ser = -1, renumb_res = -1, k = 0, amclust = 1, rotate_x = 0, rotate_y = 0, rotate_z = 0;
	ss = argv[0];
	pars_is.clear();
	for (int i = ss.size() - 1; i >= 0; --i)
    {
        if (ss[i] == '/' || ss[i] == '\\')
        {
            PROGRAMMNAME = ss.substr(i + 1, ss.size() - i - 1);
            PROGRAMMPATH = ss.substr(0, i + 1);
            i = 0;
        }
    }

	if (argc > 1)
	{
	    for (int i = 0; i < argc; ++i)
	    {
	        ss = argv[i];
            ss = ss.substr(1, ss.size() - 1);
	        if (argv[i][0] == '-' && ss == "silent")
	        {
	            SILENT = 1;
	        }
	    }
        if (!SILENT) cout << "Running PDBParser version \"" << VERSION << "\", from path \"" << PROGRAMMPATH << "\"\n";
        for (int ar = 1; ar < argc; ++ar)
		{
			if (argv[ar][0] == '-')
			{
				ss = argv[ar];
				ss = ss.substr(1, ss.size() - 1);
				if (ss == "h" || ss == "H")//Help
				{
					help(argv[0]);
					exit(0);
				}
				else if (ss == "s" || ss == "S")//Input pdb
				{
					if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
						proin = argv[ar + 1];
						pin = true;
					}
					else
					{
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "hf")//Input hetfile
				{
					if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
						hetfile = argv[ar + 1];
						hf = true;
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg = " hf";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "suf")//Suffix
				{
					if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
						SUFF = argv[ar + 1];
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg = " suf";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "os")//Output protein file
				{
					if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
						proout = argv[ar + 1];
						pout = true;
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg = " os";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "og")//Output ligand file
				{
					if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
						ligout = argv[ar + 1];
						lout = true;
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg += " og";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "bf")//Residues base file
				{
                    bf = true;
					if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
						basefile = argv[ar + 1];
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg = " bf";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "ub")
				{
				    update_base_file = true;
				}
				else if (ss == "nc")
				{
				    if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
						no_charge = strtoint(argv[ar + 1]);
						if (no_charge != 1 && no_charge != 0)
						{
						    smsg = "Wrong charge parameter value, using default (0), Charge will not be written to atoms with no charge in structure";
						    warning(scurf, smsg);
						    no_charge = 0;
						}
					}
					else
					{
						no_charge = 0;
						pars_is.push_back("There is no charge parameter value, using default (0). Charges will not be written to atoms with no charge in structure");
					}
				}
				else if (ss == "rs")//Renum serial
				{
					if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
						renumb_ser = strtoint(argv[ar + 1]);
						if (renumb_ser < 1)
						{
						    smsg = "Wrong renumbering serials starting value, using default (1)";
						    warning(scurf, smsg);
						    renumb_ser = 1;
						}
					}
					else
					{
						renumb_ser = 1;
						pars_is.push_back("There is no serial renumbering starting value, using default (1)");
					}
				}
				else if (ss == "rr")//Renum residue
				{
					if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
						renumb_res = strtoint(argv[ar + 1]);
						if (renumb_res < 0)
						{
						    smsg = "Wrong renumbering residue starting value, using default (1)";
						    warning(scurf, smsg);
						    renumb_res = 1;
						}
					}
					else
					{
						renumb_res = 1;
						pars_is.push_back("There is no renumbering residue starting value, using default (1)");
					}
				}
				else if (ss == "st")//Sort residues
				{
					sort_res = true;
				}
				else if (ss == "ct")//Cut to chains
				{
					cut = true;
				}
				else if (ss == "nw")
				{
				    NOWARNINGS = true;
				}
				else if (ss == "silent")
				{
				    SILENT = true;
				}
				else if (ss == "cif")
                {
                    if (ar + 2 < argc && argv[ar + 1][0] != '-' && argv[ar + 1][0] != '-')
					{
					    cif = true;
					    ciffile = argv[ar + 1];
					    ciflig = argv[ar + 2];
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg += " cif";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
                }
				else if (ss == "rmsd")
                {
                    if (ar + 2 < argc && argv[ar + 1][0] != '-' && argv[ar + 1][0] != '-')
					{
					    crmsd = true;
					    crmsdfile1 = argv[ar + 1];
					    crmsdfile2 = argv[ar + 2];
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg += " crmsd";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
                }
				else if (ss == "dp")
				{
				    if (ar + 2 < argc && argv[ar + 1][0] != '-' && argv[ar + 2][0] != '-')
					{
						dlg_list_file = argv[ar + 1];
						dlg_out_file = argv[ar + 2];
                        dlg_parse = true;
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg += " dp";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "dprls")
				{
				    if (ar + 2 < argc && argv[ar + 1][0] != '-' && argv[ar + 2][0] != '-')
					{
						dlg_list_file = argv[ar + 1];
						dlg_out_file = argv[ar + 2];
                        dlg_parse_runlist = true;
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg += " dprls";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "vina_dp")
				{
				    if (ar + 2 < argc && argv[ar + 1][0] != '-' && argv[ar + 2][0] != '-')
					{
						dlg_list_file = argv[ar + 1];
						dlg_out_file = argv[ar + 2];
                        vina_log_parse = true;
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg += " vina_dp";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "sdfc")
				{
				    if (ar + 2 < argc && argv[ar + 1][0] != '-' && argv[ar + 2][0] != '-')
					{
						sdf_m_file = argv[ar + 1];
						sdf_out_dir = argv[ar + 2];
                        sdfc_mode = true;
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg += " sdfc";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "vina_dprls")
				{
				    if (ar + 2 < argc && argv[ar + 1][0] != '-' && argv[ar + 2][0] != '-')
					{
						dlg_list_file = argv[ar + 1];
						dlg_out_file = argv[ar + 2];
                        vina_log_parse_runlist = true;
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg += " vina_dprls";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "ligcor")
				{
				    if (ar + 3 < argc && argv[ar + 1][0] != '-' && argv[ar + 2][0] != '-' && argv[ar + 3][0] != '-')
					{
						lig_cor_sample = argv[ar + 1];
						lig_cor_coord = argv[ar + 2];
						ligout = argv[ar + 3];
                        lig_cor = true;
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg += " ligcor";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "dlglig")
				{
				    if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{

						DCLUST = true;
						dlg_input = argv[ar + 1];
						if (ar + 3 < argc && argv[ar + 2][0] != '-' && argv[ar + 3][0] != '-')
                        {
                            dlg_output = argv[ar + 2];
                            amclust = strtoint(argv[ar + 3]);
                        }
                        else if (ar + 2 < argc && argv[ar + 2][0] != '-')
                        {
                            if (del_digits(argv[ar + 2]) == "")
                            {
                                dlg_output = argv[ar + 1];
                                amclust = strtoint(argv[ar + 2]);
                            }
                            else
                            {
                                dlg_output = argv[ar + 2];
                                amclust = 1;
                            }
                        }
                        if (amclust < 1)
                        {
                            smsg = "Wrong amclust parameter value, using default (1)";
                            warning(scurf, smsg);
                            amclust = 1;
                        }
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg += " dlglig";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "ndoe")
				{
				    NO_DLG_OUT_ERRORS = true;
				}
				else if (ss == "pb")
				{
				    POSSIBLE_BABEL = true;
				}
				else if (ss == "vmd")
				{
				    if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
					    VMDCORRECT = strtoint(argv[ar + 1]);
					    if (VMDCORRECT > 2 || VMDCORRECT < 0)
                        {
                            smsg = "Wrong vmd parameter value, using default (0)";
                            warning(scurf, smsg);
                            VMDCORRECT = 0;
                        }
					}
					else
					{
						VMDCORRECT = 1;
					}
				}
				else if (ss == "renat")
				{
				    if (ar + 2 < argc && argv[ar + 1][0] != '-' && argv[ar + 2][0] != '-')
					{
					    renatpdb = 1;
					    proren = argv[ar + 1];
					    proout_ren = argv[ar + 2];

					}
					else
					{
						smsg = "Missing or wrong argument";
					    smsg += " renat";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "rotate")
                {
                    if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
					    rotate_pdb = 1;
					    rotate_x = strtoint(argv[ar + 1]);
					    if (rotate_x > 180 || rotate_x < 180)
                        {
                            smsg = "Rotate x param is not in 2PI diapason";
                            warning(scurf, smsg);
                        }
                        if (ar + 2 < argc && argv[ar + 2][0] != '-')
                        {
                            rotate_y = strtoint(argv[ar + 2]);
                            if (rotate_y > 180 || rotate_y < 180)
                            {
                                smsg = "Rotate y param is not in 2PI diapason";
                                warning(scurf, smsg);
                            }
                        }
                        if (ar + 3 < argc && argv[ar + 3][0] != '-')
                        {
                            rotate_z = strtoint(argv[ar + 3]);
                            if (rotate_z > 180 || rotate_z < 180)
                            {
                                smsg = "Rotate z param is not in 2PI diapason";
                                warning(scurf, smsg);
                            }
                        }
					}
					else
					{
						smsg = "Missing or wrong argument";
					    smsg += " rotate";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
                }
				else if (ss == "clr")
				{
				    RREN_CHLOR = true;
				}
				else if (ss == "nooc")
				{
				    NO_OC = true;
				}
				else if (ss == "setoctf1")
				{
				    SET_OC_TF_1 = true;
				}
				else if (ss == "delallh")
				{
				    DEL_ALL_HYDROGENS = true;
				}
				else if (ss == "notf")
				{
				    NO_TF = true;
				}
				else if (ss == "noter")
				{
				    NO_TER = true;
				}
				else if (ss == "zn2")
				{
				    RREN_ZINK = true;
				}
				else if (ss == "renres")
                {
                    RES_RENAME = true;
                }
				else if (ss == "lgc")
				{
				    LIG_GEOMETRIC_CENTER = true;
				}
				else if (ss == "frep")
				{
				    REPLASE_EX_FILE = true;
				}
				else if (ss == "addch")
				{
				    ADD_CHARGE = true;
				}
				else if (ss == "chc")
				{
				    CORRECT_HETATM_CHAINS = true;
				}
				else if (ss == "fasta")
				{
				    if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
						fasta_file = argv[ar + 1];
						FASTA = true;
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg += " fasta";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "lmol2")
				{
				    if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
						lmol2_file = argv[ar + 1];
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg += " lmol2";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "ligf")
				{
				    if (ar + 1 < argc && argv[ar + 1][0] != '-')
					{
						ligname = argv[ar + 1];
						if (ar + 2 < argc && argv[ar + 2][0] != '-')
                        {
                            ligch = argv[ar + 2][0];
                        }
                        else
                        {
                            ligch = '?';
                        }
					}
					else
					{
					    smsg = "Missing or wrong argument";
					    smsg += " ligf";
					    error_msg(scurf, smsg);
						help(argv[0]);
						exit(0);
					}
				}
				else if (ss == "evp")
				{
				    if (ar + 3 < argc && argv[ar + 1][0] != '-' && argv[ar + 2][0] != '-' && argv[ar + 3][0] != '-')
					{
                        roc_plot = true;
					    nrp_parf = argv[ar + 1];
                        nrp_resf = argv[ar + 2];
                        plot_file = argv[ar + 3];
					}
					else
					{
					    smsg = "There are no files for ROC plot calculating";
					    smsg += " evp";
					    warning(scurf, smsg);
					}
				}
				else if (ss == "ac")
				{
                    ALTERNATIVE_COORDINATES = 1;
				}
				else if (ss == "igseg")
				{
                    IGNORE_SEGMENT = 1;
				}
				else
				{
					smsg = "Unknown parameter";
					error_msg(scurf, smsg);
					help(argv[0]);
					exit(0);
				}
			}
		}
	}
	else
	{
	    //cout << "Running PDBParser version \"" << VERSION << "\" from path \"" << PROGRAMMPATH << "\"\n";
        smsg = "No parameters. Run terminated";
        error_msg(scurf, smsg);
		help(argv[0]);
		exit(0);
	}
	if (pars_is.size() > 0)
    {
        if (!SILENT) cout << "Parsing input parameters\n";
        for (int i = 0; i < pars_is.size(); ++i)
        {
            if (!SILENT) cout << "  " << pars_is[i] << "\n";
        }
    }
	if (bf)
	{
		if (!SILENT) cout << "Parsing base file ...\n";
		parse_res_base(res_base, basefile);
	}
    if (cif)
    {
        int num;
        if (!SILENT) cout << "Making pdb from cif file ...\n";
        if (lout)
        {
            num = make_cif_pdb(ciffile, ciflig, ligout);
            if (!SILENT) cout << "There are about " << num << " stereo tautomers\n";
        }
        else
        {
            smsg = "There is no ligand out file parameter specified. Terminated ...";
            error_msg(scurf, smsg);
        }
        return 0;
    }
    else if (crmsd)
    {
        double rmsd = lig_RMSD(crmsdfile1, crmsdfile2);
        if (rmsd < 0)
        {
            smsg = "There were fatal errors, terminated";
            warning(scurf, smsg);
        }
        else
        {
            cout.precision(4);
            cout << "RMSD = " << rmsd<< " A\n";
        }
		return 0;
    }
    else if (sdfc_mode)
    {
        cut_sdf_to_mol2(sdf_m_file, sdf_out_dir);
        return 0;
    }
	else if (lig_cor)
	{
        double globe_rmsd;
        if (!SILENT) cout << "Correcting ligand ...\n";
        CORRECT_HETATM_CHAINS = true;
        globe_rmsd = lig_correcter(lig_cor_sample, lig_cor_coord, ligout);
        if (globe_rmsd < 0)
        {
            smsg =  "There were fatal errors, terminated";
            warning(scurf, smsg);
        }
        else
        {
            cout.precision(3);
            if (!SILENT) cout << "  Summ deviation is " << globe_rmsd << "A" << "\n";
        }
		return 0;
	}
	else if (dlg_parse || dlg_parse_runlist)
    {
        if (!SILENT) cout << "Parsing DLG files ...\n";
        redock = dlg_parser(dlg_list_file, dlg_out_file, roc_plot, dlg_parse_runlist);
        if (redock.size() > 0)
        {
            if (!SILENT) cout << "Redocking for few runs is required ...\n";
            make_redock_runlist(dlg_list_file, redock);
        }
		return 0;
    }
    else if (vina_log_parse || vina_log_parse_runlist)
    {
        if (!SILENT) cout << "Parsing VINA log files ...\n";
        redock = vina_log_parser(dlg_list_file, dlg_out_file, roc_plot, vina_log_parse_runlist);
        //if (!SILENT) cout << "!!!\n";
        if (redock.size() > 0)
        {
            if (!SILENT) cout << "Redocking for few runs is required ...\n";
            make_redock_runlist(dlg_list_file, redock);
        }
		return 0;
    }
    else if (DCLUST)
    {
        if (!SILENT) cout << "Preparing LIG files from DLG file ...\n";
        get_lig_dlg(dlg_input, dlg_output, amclust);
		return 0;
    }
    else if (roc_plot)
    {
        if (!SILENT) cout << "Analyzing docking results, making ROC plots and calculating AUC values ...\n";
        rampce_dock_results(nrp_parf, nrp_resf, plot_file);
        return 0;
    }
    else if (lmol2_file != "")
    {
        if (!SILENT) cout << "Parsing mol2 ligand file ...\n";
        mol2file ml2f;
        parse_mol2_file(lmol2_file, ml2f);
        if (LIG_GEOMETRIC_CENTER)
        {
            k = lmol2_file.size() - 1;
            postf = "";
            while (k > 0 && lmol2_file[k] != '.')
            {
                k--;
            }
            if (k == 0)
            {
                postf = ".pdb";
                k = lmol2_file.size();
            }
            else
            {
                postf = lmol2_file.substr(k, lmol2_file.size() - k);
            }
            if (!SILENT) cout << "  Writing geometric coordinates of the ligand\n";
            geometric_coords_mf(ml2f, lmol2_file.substr(0, k) + "_center.txt", "");
        }
        return 0;
    }
	else if (pin)
	{
		if (!SILENT) cout << "Parsing PDB file ...\n";
		parse_pdb(proin, pdbstruct);
		if (SET_OC_TF_1)
        {
            set_oc_to_1(pdbstruct);
        }
        if (DEL_ALL_HYDROGENS)
        {
            delete_all_hydrogens(pdbstruct);
        }
        if (VMDCORRECT != 0)
        {
            if (!SILENT) cout << "Correcting VMD output ...\n";
            correct_VMD(pdbstruct);
        }
		if (!CORRECT_HETATM_CHAINS)
        {
            if (check_chains(pdbstruct))
            {
                smsg = "There were founded chains with no atoms in ATOM field";
                warning(scurf, smsg);
                CORRECT_HETATM_CHAINS = true;
            }
        }
		if (CORRECT_HETATM_CHAINS)
        {
            if (!SILENT) cout << "Fixing HETATM chain bugs...\n";
            correct_h_chains(pdbstruct);
        }
        if (FASTA)
        {
            if (!SILENT) cout << "Making fasta from pdb structure ...\n";
            make_fasta(pdbstruct, fasta_file);
        }
        if (ALTERNATIVE_COORDINATES == 1)
        {
            if (!SILENT) cout << "Removing alternative coordinates ...\n";
            remove_alternative(pdbstruct);
        }
		if (update_base_file)
		{
		    if (!SILENT) cout << "Updating base using current receptor structure ...\n";
		    update_base(pdbstruct, res_base);
		    if (!SILENT) cout << "Writing base to file ...\n";
		    out_res_base(res_base, basefile);
		    return 0;
		}
		if (rotate_pdb)
        {
            rotate_pdb_coords(pdbstruct, rotate_x, rotate_y, rotate_z);
        }
		if (sort_res)
		{
			if (bf)
			{
				if (!SILENT) cout << "Sorting atoms using BASE file and renumbering serials ...\n";
				sort_pdb(pdbstruct, res_base);
			}
			else
			{
				smsg = "There is no base file parameter. Can't sort";
				warning(scurf, smsg);
			}
		}
        if (renatpdb)
        {
            smsg = "Renaming atom function is not implemented yet\n";
            error_msg(scurf, smsg);
            exit(1);
            pdb_file pdbstruct_ren;
            if (!SILENT) cout << "Reading 2nd pdb structure ...\n";
            parse_pdb(proren, pdbstruct_ren);
            if (!SILENT) cout << "Renaming atoms in 2nd pdb according to the names in the 1st structure with";
            if (!sort_res)
            {
                if (!SILENT) cout << "out";
            }
            if (!SILENT) cout << " sorting atoms basing on 1st structure ...\n";
            rename_atoms_pat(pdbstruct, pdbstruct_ren, sort_res);
            k = proout_ren.size() - 1;
            postf = "";
            while (k > 0 && proout_ren[k] != '.')
            {
                k--;
            }
            if (k == 0)
            {
                postf = ".pdb";
                k = proout_ren.size();
            }
            else
            {
                postf = proout_ren.substr(k, proout_ren.size() - k);
            }
            if (!SILENT) cout << "Writing 2nd PDB file with";
            if (renumb_ser > 0)
            {
                renum_ser(pdbstruct, renumb_ser, true);
            }
            else
            {
                if (!SILENT) cout << "out";
            }
            if (!SILENT) cout << " renumbering serials and with";
            if (renumb_res >= 0)
            {
                renum_res(pdbstruct, renumb_res);
            }
            else
            {
                if (!SILENT) cout << "out";
            }
            if (!SILENT) cout << " renumbering residues ...\n";
            out_pdb(pdbstruct, proout_ren.substr(0, k) + SUFF + postf, true);
            return 0;
        }
		if (hf)
		{
			map<char, set<string> > delhetlist;
			set<string> deleted_hets;
			if (!SILENT) cout << "Parsing HETATM file ...\n";
			vector <string> stdligands;
			map <string, string> renresi;
			map < pair < string, string > , int> achs;
			plig = parse_het_file(hetfile, stdligands, delhetlist, ligname, ligch, renresi, achs);

			if (RES_RENAME || RREN_ZINK || RREN_CHLOR)
			{
                if (RREN_ZINK)
                {
                    renresi["ZN"] = "ZN2";
                    achs[make_pair("ZN2", "ZN")] = 2;
                    achs[make_pair("ZN2", "zn")] = 2;
                    achs[make_pair("ZN2", "Zn")] = 2;
                }
			    if (RREN_CHLOR)
                {
                    renresi["CL"] = "CLA";
                    achs[make_pair("CLA", "CL")] = -1;
                    achs[make_pair("CLA", "cl")] = -1;
                    achs[make_pair("CLA", "Cl")] = -1;
                }

			    if (!SILENT) cout << "Renaming residues\n";
			    rename_res_add_charge(pdbstruct, renresi, achs);
			}
			if (cut)
			{
			    if (!SILENT) cout << "Analyzing ligands from pdb in CUT mode\n";
			    if (plig)
			    {
			        if (lout)
			        {
                        if (!SILENT) cout << "  There is LIG beginning string in HETFILE\n";
			            cutting_lig(pdbstruct, ligname, ligout, stdligands, renumb_ser, renumb_res, false);//Удаляя лиганды из структуры
					    delhetlist['!'].insert(ligname);
			        }
			        else
                    {
                        if (!SILENT) cout << "  There is LIG beginning string in HETFILE, but no ligand output file\n";
                    }
			    }
			    else
                {
                    if (lout)
                    {
                        smsg += "There is ligand output file (og), but no LIG beginning string in HETFILE";
                        warning(scurf, smsg);
                    }
                    else
                    {
                        if (!SILENT) cout << "  There is no LIG beginning string in HETFILE\n";
                    }
                }
			}
			else
			{
			    if (!SILENT) cout << "Analysing ligands in structure ...\n";
                if (plig)
                {
                    if (!SILENT) cout << "  There is LIG beginning string in HETFILE\n";
                    if (ligname == "FIND")
                    {
                        if (!SILENT) cout << "  Trying to FIND non-standard ligand in chain " << ligch << "\n";
                        call = find_lig(pdbstruct, stdligands, ligname, ligch);
                        if (call == 0)
                        {
                            if (!SILENT) cout << "  There is no non-standard ligands in this structure\n";
                        }
                        else if (call == 1)
                        {
                            if (!SILENT) cout << "  There is 1 non-standard ligand " << ligname << " in chain " << ligch << "\n";
                        }
                        else if (call == 2)
                        {
                            if (!SILENT) cout << "  There is more than 1 non-standard ligand in this structure, using " << ligname << " from chain " << ligch << "\n";
                        }
                        else
                        {
                            smsg = "Wierd callback from function find_lig";
                            warning(scurf, smsg);
                        }
                        for (map < char, set < string > >::iterator it = delhetlist.begin(); it != delhetlist.end(); ++it)
                        {
                            if (it->second.count("FIND"))
                            {
                                it->second.erase(it->second.find("FIND"));
                                it->second.insert(ligname);
                            }
                        }
                    }
                    if (lout)
                    {
				        pdb_file ligand;
				        k = ligout.size() - 1;
                        postf = "";
                        while (k > 0 && ligout[k] != '.')
                        {
                            k--;
                        }
                        if (k == 0)
                        {
                            postf = ".pdb";
                            k = ligout.size();
                        }
                        else
                        {
                            postf = ligout.substr(k, ligout.size() - k);
                        }
                        if (!SILENT) cout << "Making ligand structure ...\n";
                        if (!make_lig(ligand, ligname, ligch, pdbstruct))
                        {
                            smsg += "There is no ligand ";
                            smsg += ligname;
                            smsg += " in chain ";
                            smsg += ligch;
                        }
                        else
                        {
                            if (LIG_GEOMETRIC_CENTER)
                            {
                                if (!SILENT) cout << "  Writing geometric coordinates of the ligand\n";
                                geometric_coords(ligand, ligout.substr(0, k) + "_center.txt", ligname);
                            }
                            if (!SILENT) cout << "  Writing ligand file with";
                            if (renumb_ser > 0)
                            {
                                renum_ser(ligand, 1, false);
                            }
                            else
                            {
                                if (!SILENT) cout << "out";
                            }
                            if (!SILENT) cout << " renumbering serials and with";
                            if (renumb_res >= 0)
                            {
                                renum_res(ligand, 1);
                            }
                            else
                            {
                                if (!SILENT) cout << "out";
                            }
                            if (!SILENT) cout << " renumbering residues ...\n";
                            out_pdb(ligand, ligout.substr(0, k) + SUFF + postf, false);
                        }
				    }
                }
			    else
                {
                    if (!SILENT) cout << "  There is no LIG beginning string in HETFILE\n";
                }
			}
			if (!SILENT) cout << "Deleting HETATMs ...\n";
			deleted_hets = del_het(pdbstruct, delhetlist);
			if (!SILENT) cout << "  Deleted HETATMs:";
			for (set<string>::iterator it = deleted_hets.begin(); it != deleted_hets.end(); ++it)
            {
                if (!SILENT) cout << " " << *it;
            }
            if (!SILENT) cout << "\n";
		}
		else
		{
			if (!SILENT) cout << "There is no HETFILE parameter. Do nothing about ligand.\n";
		}
		if (pout)
		{
			if (cut)
			{
				if (!SILENT) cout << "Writing PDB files in CUT mode with";
				if (!SILENT) cout << (renumb_res != -1 ? "" : "out");
				if (!SILENT) cout << " renumbering residues and with";
				if (!SILENT) cout << (renumb_ser != -1 ? "" : "out");
				if (!SILENT) cout << " renumbering serials ...\n";
				cutting_prot(pdbstruct, proout, renumb_res, renumb_ser);
			}
			else
			{
			    k = proout.size() - 1;
                postf = "";
                while (k > 0 && proout[k] != '.')
                {
                    k--;
                }
                if (k == 0)
                {
                    postf = ".pdb";
                    k = proout.size();
                }
                else
                {
                    postf = proout.substr(k, proout.size() - k);
                }
				if (!SILENT) cout << "Writing PDB file with";
                if (renumb_ser > 0)
                {
                    renum_ser(pdbstruct, renumb_ser, true);
                }
                else
                {
                    if (!SILENT) cout << "out";
                }
                if (!SILENT) cout << " renumbering serials and with";
                if (renumb_res >= 0)
                {
                    renum_res(pdbstruct, renumb_res);
                }
                else
                {
                    if (!SILENT) cout << "out";
                }
                if (!SILENT) cout << " renumbering residues ...\n";
				out_pdb(pdbstruct, proout.substr(0, k) + SUFF + postf, true);
			}
		}
	}
	else
	{
		if (!SILENT) cout << "There is no PDB file input\n";
	}
	if (!SILENT) cout << "\n";
	return 0;
}
