/* 
 * File:   elementsRepo.h
 * Author: Никита
 *
 * Created on 29 октября 2015 г., 15:35
 */

#ifndef ELEMENTSREPO_H
#define	ELEMENTSREPO_H

#include <iostream>
#include <map>
using namespace std;

#define ELEMENTS_FILE "elements.dat"


class elementsRepo {
public:
    void readDatFile();
    double getValue(string s);
    
    static elementsRepo * Instanse() {
        if(!inst){           
            inst = new elementsRepo();
        }
        return inst;
    }
private:
    static elementsRepo * inst;
    elementsRepo();
    map<string, double> masses;


};

#endif	/* ELEMENTSREPO_H */