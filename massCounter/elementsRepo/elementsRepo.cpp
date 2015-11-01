/*
 * File:   elementsRepo.cpp
 * Author: Никита
 *
 * Created on 29 октября 2015 г., 15:35
 */

#include "elementsRepo.h"
#include "../../common/commonFuncs.h"
#include "../../common/logger/logs.h"

using namespace std;



elementsRepo::elementsRepo(){
    readDatFile();
}
void elementsRepo::readDatFile(){
    ifstream din (ELEMENTS_FILE);

    string s, element;
    double val;
    vector <string> splitedStr;
	getline(din, s); //header
    while (getline(din, s)){
        splitedStr = split(s, "\t");
        element = toLowerCase(splitedStr[1]);
        log(INFO_MSG, "readDatFile", element);
        //cout << element << " ";
        if (masses.count(element)){
			log(ERROR_MSG, "readDatFile", "Double " + element + " entry");
        } else {
			val = strtodoub(splitedStr[3]);
            masses[element] = val;
        }
    }

    din.close();
}

double elementsRepo::getValue(string s){
	s = toLowerCase(s);
    if (masses.count(s)){
        return masses[s];
    } else{

		log(ERROR_MSG, "getValue", "No entry for |" + s + "|\n");
    }
    return 0;
}




elementsRepo* elementsRepo::inst = 0;
