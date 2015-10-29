/* 
 * File:   elementsRepo.cpp
 * Author: Никита
 * 
 * Created on 29 октября 2015 г., 15:35
 */

#include <fstream>
#include <vector>

#include "elementsRepo.h"
#include "commonFuncs.h"

using namespace std;



elementsRepo::elementsRepo(){
    readDatFile();
}
void elementsRepo::readDatFile(){
    ifstream din (ELEMENTS_FILE);

    string s;
    double val;
    vector <string> splitedStr;

    while (getline(din, s)){
        splitedStr = split(s, "\t");
        if (masses.count(splitedStr[1])){
            cerr << "Double " << splitedStr[1] << " entry\n";
        } else {
                val = strtodoub(splitedStr[3]);
            masses[splitedStr[1]] = val;
        }
    }

    din.close();
}

double elementsRepo::getValue(string s){
    if (masses.count(s)){
        return masses[s];
    } else{
        cerr << "No entry for " << s << "\n";
    }
    return 0;
}




elementsRepo* elementsRepo::inst = 0;