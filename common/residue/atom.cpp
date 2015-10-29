/* 
 * File:   atom.cpp
 * Author: Никита
 * 
 * Created on 29 октября 2015 г., 14:59
 */

#include "atom.h"

atom::atom() {
    charge = 0;
    element = "";
    serial = 0;
}

atom::atom(const atom& orig) {
    charge = orig.charge;
    element = orig.element;
    serial = orig.serial;
    x = orig.x;
    y = orig.y;
    z = orig.z;
}

atom::~atom() {
}

