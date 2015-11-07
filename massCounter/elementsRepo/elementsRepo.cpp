/*
 * File:   elementsRepo.cpp
 * Author: Никита
 *
 * Created on 29 октября 2015 г., 15:35
 */

#include "elementsRepo.h"
#include "../../common/commonFuncs.h"
#include "../../common/logger/logs.h"
#include "elements.h"

using namespace std;



elementsRepo::elementsRepo(){
    fillElementsRepo();
}
void elementsRepo::fillElementsRepo(){
    vector<string> elementsDat;
    fillElementsDat(elementsDat);

    string s, element;
    double val;
    vector <string> splitedStr;
    for (int i = 0; i < elementsDat.size(); ++i)
    {
        splitedStr = split(elementsDat[i], "\t");
        element = toLowerCase(splitedStr[1]);
        log(INFO_MSG, "fillElementsRepo", element);
        //cout << element << " ";
        if (masses.count(element)){
			log(ERROR_MSG, "fillElementsRepo", "Double " + element + " entry");
        } else {
			val = strtodoub(splitedStr[3]);
            masses[element] = val;
        }
    }
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
