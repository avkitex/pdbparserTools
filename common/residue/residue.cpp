/* 
 * File:   residue.cpp
 * Author: Никита
 * 
 * Created on 29 октября 2015 г., 14:53
 */

#include "residue.h"
#include "elementsRepo.h"

residue::residue() {
    atoms.clear();
    ligName = "";
    resnumber = 0;
    chain = '?';
    mass = -1;
}

void residue::countMass(){
    mass = 0;
    for (int i = 0; i < atoms.size(); ++i){
        mass += elementsRepo::Instanse()->getValue(atoms[i].getElement());
    }
}

double residue::getMass() {
    if (mass < 0,05){//zero
        countMass();
    }
    return mass;
}