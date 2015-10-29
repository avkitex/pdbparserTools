/* 
 * File:   atom.h
 * Author: Никита
 *
 * Created on 29 октября 2015 г., 14:59
 */

#ifndef ATOM_H
#define	ATOM_H

#include <iostream>
#include "point3D.h"

using namespace std;

class atom: private point3D {
public:
    atom();
    atom(const atom& orig);
    virtual ~atom();
    
    string getElement(){
        return element;
    }
protected:
    int serial;
    string element;
    int charge;
};

#endif	/* ATOM_H */

