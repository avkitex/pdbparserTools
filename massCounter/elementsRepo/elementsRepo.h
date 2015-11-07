/*
 * File:   elementsRepo.h
 * Author: Никита
 *
 * Created on 29 октября 2015 г., 15:35
 */

#ifndef ELEMENTSREPO_H
#define	ELEMENTSREPO_H

#include <iostream>
#include <vector>
#include <map>

using namespace std;


class elementsRepo {
public:
    void fillElementsRepo();
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
