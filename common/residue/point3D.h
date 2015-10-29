/* 
 * File:   point3D.h
 * Author: Никита
 *
 * Created on 29 октября 2015 г., 14:54
 */

#ifndef POINT3D_H
#define	POINT3D_H

class point3D {
public:
    point3D();
    point3D(const point3D& orig);
    virtual ~point3D();
protected:
    double x, y, z;
};

#endif	/* POINT3D_H */

