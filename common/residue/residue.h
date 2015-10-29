/* 
 * File:   residue.h
 * Author: Никита
 *
 * Created on 29 октября 2015 г., 14:53
 */

#ifndef RESIDUE_H
#define	RESIDUE_H

#include <vector>
#include "atom.h"

using namespace std;

class residue: private atom {
public:
    residue();
    
    double getMass();
    void countMass();
private:
    char chain;
    int resnumber;
    string ligName;
    vector <atom> atoms;
    double mass;
            
};

#endif	/* RESIDUE_H */

